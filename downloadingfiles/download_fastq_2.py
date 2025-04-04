#!/usr/bin/env python3
import os
import sys
import time
import argparse
import subprocess
import shutil
import glob
import gzip
import fcntl
import psutil
import atexit 
from concurrent.futures import ProcessPoolExecutor, as_completed
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

##########################
# Locking & Logging Helpers
##########################

def flock_exclusive(fd):
    """Acquire an exclusive (write) lock on an open file descriptor."""
    fcntl.flock(fd, fcntl.LOCK_EX)

def flock_release(fd):
    """Release a lock on an open file descriptor."""
    fcntl.flock(fd, fcntl.LOCK_UN)

def append_line_with_lock(line, out_file, lock_file):
    """
    Append one line to 'out_file' using an exclusive lock,
    so multiple processes do not overwrite each other.
    """
    try:
        with open(lock_file, "a+") as lk:
            flock_exclusive(lk)
            with open(out_file, "a") as f:
                f.write(line + "\n")
            flock_release(lk)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Could not append '{line}' to {out_file}: {e}\n")

def log_debug_message(msg, debug_file, debug_lock):
    """
    Write a debug message to 'debug_file' with an exclusive lock.
    """
    try:
        with open(debug_lock, "a+") as lk:
            flock_exclusive(lk)
            with open(debug_file, "a") as dbg:
                dbg.write(f"{time.ctime()} {msg}\n")
            flock_release(lk)
    except Exception as e:
        sys.stderr.write(f"[WARN] Could not log debug message '{msg}': {e}\n")

##########################
# File Validation Helpers
##########################

def is_valid_fastq(path, min_size_bytes=1024):
    """Check if file exists and is at least min_size_bytes in size."""
    return os.path.isfile(path) and os.path.getsize(path) >= min_size_bytes

def is_valid_gzip(path):
    """Check if gzip file can be read (at least 1 byte)."""
    try:
        with gzip.open(path, "rb") as f:
            f.read(1)
        return True
    except:
        return False

def remove_file_safely(path, debug_file, debug_lock):
    """Remove file if it exists, log any error."""
    if os.path.exists(path):
        try:
            os.remove(path)
        except Exception as e:
            log_debug_message(f"[ERROR] remove_file_safely({path}): {e}", debug_file, debug_lock)

##########################
# Command Runner with Retry
##########################

@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=2, max=10),
    retry=retry_if_exception_type(RuntimeError)
)
def run_command(cmd_list, err_msg, debug_file, debug_lock):
    """
    Run a shell command with check=True. On failure, log stderr and raise RuntimeError.
    We retry up to 3 times, with exponential backoff, using tenacity.
    """
    try:
        result = subprocess.run(cmd_list, check=True, text=True, stderr=subprocess.PIPE)
        return result.stdout
    except subprocess.CalledProcessError as e:
        msg = f"{err_msg}:\n{e.stderr}"
        log_debug_message(msg, debug_file, debug_lock)
        # Raise a RuntimeError so tenacity sees it and retries.
        raise RuntimeError(msg) from e

##########################
# FASTQ Steps
##########################

def compress_fastqs(accession, fastq_dir, debug_file, debug_lock):
    """Gzip or pigz all .fastq files for this accession."""
    pattern = os.path.join(fastq_dir, f"{accession}*.fastq")
    fastq_files = glob.glob(pattern)
    compressor = shutil.which("pigz") or shutil.which("gzip")
    if not compressor:
        raise RuntimeError("No gzip or pigz found on system!")

    for fq in fastq_files:
        if not fq.endswith(".gz"):
            cmd = [compressor]
            if "pigz" in compressor:
                cmd.extend(["-p", "4"])  # e.g. 4 threads
            cmd.append(fq)
            run_command(cmd, f"[compress_fastqs] {accession} compression failed on {fq}",
                        debug_file, debug_lock)

def cleanup_invalid_fastqs(accession, fastq_dir, debug_file, debug_lock):
    """Remove any .fastq.gz that appear invalid or too small."""
    gz_files = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz"))
    for gz in gz_files:
        if not is_valid_fastq(gz) or not is_valid_gzip(gz):
            log_debug_message(f"[WARN] Removing invalid FASTQ: {gz}",
                              debug_file, debug_lock)
            try:
                os.remove(gz)
            except Exception as e:
                log_debug_message(f"[ERROR] while removing {gz}: {e}",
                                  debug_file, debug_lock)

def classify_fastq_by_read_type(accession, fastq_dir, debug_file, debug_lock):
    """
    Runs 'vdb-dump <accession> -C READ_TYPE' to figure out BIOLOGICAL vs
    TECHNICAL reads, then moves them to subfolders.
    """
    vdb_dump = shutil.which("vdb-dump")
    if not vdb_dump:
        log_debug_message(f"[classify] vdb-dump not found; skipping for {accession}.",
                          debug_file, debug_lock)
        return

    cmd_read_type = f"vdb-dump {accession} -C READ_TYPE 2>/dev/null | head -n 1"
    try:
        result = subprocess.run(cmd_read_type, shell=True, capture_output=True, text=True, check=True)
        line = result.stdout.strip()
        if "READ_TYPE:" not in line:
            log_debug_message(f"[classify] No READ_TYPE info for {accession}.", debug_file, debug_lock)
            return

        # e.g. "READ_TYPE: BIOLOGICAL, BIOLOGICAL" => ["BIOLOGICAL", "BIOLOGICAL"]
        types_part = line.split(":", 1)[1].strip()
        read_types = [t.strip() for t in types_part.split(",")]

        bio_dir = os.path.join(fastq_dir, "fastq_biologicaldata")
        tech_dir = os.path.join(fastq_dir, "fastq_technicaldata")
        os.makedirs(bio_dir, exist_ok=True)
        os.makedirs(tech_dir, exist_ok=True)

        for i, read_type in enumerate(read_types):
            gz_suffixed = f"{accession}_{i+1}.fastq.gz"
            src_suffixed = os.path.join(fastq_dir, gz_suffixed)

            src_path = None
            if os.path.exists(src_suffixed):
                src_path = src_suffixed
                gz_filename = gz_suffixed
            elif i == 0:
                # Single-end might be unsuffixed
                gz_unsuffixed = f"{accession}.fastq.gz"
                src_unsuffixed = os.path.join(fastq_dir, gz_unsuffixed)
                if os.path.exists(src_unsuffixed):
                    src_path = src_unsuffixed
                    gz_filename = gz_unsuffixed

            if not src_path:
                raise FileNotFoundError(f"[classify] Missing FASTQ for read #{i+1} in {accession}")

            dest_dir = bio_dir if "BIOLOGICAL" in read_type else tech_dir
            dest_path = os.path.join(dest_dir, gz_filename)
            log_debug_message(f"[classify] Move {src_path} -> {dest_path}",
                              debug_file, debug_lock)
            shutil.move(src_path, dest_path)

    except subprocess.CalledProcessError as e:
        log_debug_message(f"[classify] vdb-dump failed for {accession}:\n{e.stderr}",
                          debug_file, debug_lock)
    except Exception as e:
        log_debug_message(f"[classify] Error for {accession}: {e}",
                          debug_file, debug_lock)

##########################
# Metadata Fetch
##########################

def fetch_sra_metadata(accession, metadata_dir, combined_meta, debug_file, debug_lock):
    """
    Uses esearch/efetch to get run-info CSV, merges it into a combined TSV file
    with exclusive lock on the debug lock (for simplicity).
    """
    csv_path = os.path.join(metadata_dir, f"{accession}-run-info.csv")
    # If user has set EMAIL or NCBI_API_KEY in environment, esearch will pick that up
    # Or we can pass -email explicitly. For EDirect, environment "EMAIL" is typically enough.
    cmd = f'esearch -db sra -query "{accession}" | efetch -format runinfo > "{csv_path}"'
    ret = subprocess.call(cmd, shell=True)
    if ret != 0 or not os.path.isfile(csv_path) or os.path.getsize(csv_path) == 0:
        # Not fatal, but log it
        log_debug_message(f"[metadata] Failed to fetch runinfo for {accession} (ret={ret})",
                          debug_file, debug_lock)
        return

    # Convert CSV -> TSV
    tsv_path = csv_path.replace(".csv", ".tsv")
    try:
        with open(csv_path, "r") as inf, open(tsv_path, "w") as outf:
            for line in inf:
                outf.write(line.replace(",", "\t"))

        # Append to combined file
        os.makedirs(os.path.dirname(combined_meta), exist_ok=True)
        with open(debug_lock, "a+") as lk:
            flock_exclusive(lk)
            # If combined file doesn't exist, copy the header
            if not os.path.isfile(combined_meta):
                with open(tsv_path, "r") as t_in:
                    header_line = t_in.readline()
                with open(combined_meta, "w") as cm:
                    cm.write(header_line)
            # Now append rows except header
            with open(tsv_path, "r") as t_in, open(combined_meta, "a") as cm:
                _ = t_in.readline()  # skip header
                for row in t_in:
                    cm.write(row)

            flock_release(lk)

        # Cleanup CSV/TSV
        remove_file_safely(csv_path, debug_file, debug_lock)
        remove_file_safely(tsv_path, debug_file, debug_lock)

    except Exception as e:
        log_debug_message(f"[metadata] Error merging runinfo for {accession}: {e}",
                          debug_file, debug_lock)

##########################
# SRA Download Route
##########################

def sra_download_route(accession, workdir, debug_file, debug_lock, combined_meta):
    """
    1) prefetch .sra
    2) vdb-validate
    3) fasterq-dump (with adaptive memory)
    4) fetch SRA metadata
    5) ensure we got .fastq
    Returns True if success, False if error
    """
    fastq_dir = os.path.join(workdir, "fastq_data")
    metadata_dir = os.path.join(workdir, "metadata")

    # Step 0: check memory at runtime
    mem_bytes = psutil.virtual_memory().available
    mem_gb = int(mem_bytes / (1024**3))
    mem_str = f"{mem_gb}G"

    sra_file = os.path.join(fastq_dir, f"{accession}.sra")

    # 1) prefetch
    try:
        run_command(["prefetch", accession, "--max-size", "200G", "--output-file", sra_file],
                    f"[sra_download_route] prefetch failed {accession}",
                    debug_file, debug_lock)
    except RuntimeError:
        return False

    # 2) vdb-validate
    try:
        subprocess.run(["vdb-validate", sra_file], check=True,
                       capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        log_debug_message(f"[sra_download_route] vdb-validate failed {accession}:\n{e.stderr}",
                          debug_file, debug_lock)
        return False

    # 3) fasterq-dump
    try:
        run_command([
            "fasterq-dump", sra_file,
            "--outdir", fastq_dir,
            "--threads", str(os.cpu_count()),
            "--mem", mem_str,
            "--split-files",
            "--include-technical"
        ], f"[sra_download_route] fasterq-dump failed {accession}",
           debug_file, debug_lock)
    except RuntimeError:
        return False

    # 4) fetch metadata
    fetch_sra_metadata(accession, metadata_dir, combined_meta, debug_file, debug_lock)

    # 5) check if we got .fastq
    fastqs = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq"))
    if not fastqs:
        log_debug_message(f"[sra_download_route] No FASTQs for {accession}", debug_file, debug_lock)
        return False
    # minimal size check
    for fq in fastqs:
        if not is_valid_fastq(fq, 1024):
            log_debug_message(f"[sra_download_route] Invalid or tiny FASTQ: {fq}",
                              debug_file, debug_lock)
            return False

    return True

##########################
# Worker Function
##########################

def process_accession(accession, args):
    """
    Process one accession:
      - run sra_download_route
      - compress and classify fastqs
      - remove .sra
      - if success => append to completed
      - if fail => append to failed
    """
    debug_lock = args.debug_lock
    debug_file = args.debug_file
    completed_file = args.completed_file
    completed_lock = args.completed_lock
    failed_file = args.failed_file
    failed_lock = args.failed_lock
    workdir = args.workdir

    try:
        # 1) Download from SRA
        ok = sra_download_route(accession, workdir, debug_file, debug_lock, args.combined_metadata)
        if not ok:
            raise RuntimeError(f"[process_accession] SRA route failed for {accession}")

        # 2) Compress
        compress_fastqs(accession, os.path.join(workdir, "fastq_data"), debug_file, debug_lock)

        # 3) Validate GZ
        gz_files = glob.glob(os.path.join(workdir, "fastq_data", f"{accession}*.fastq.gz"))
        if not gz_files:
            raise RuntimeError(f"[process_accession] No .gz after compression for {accession}")
        for gz in gz_files:
            if not is_valid_gzip(gz):
                raise RuntimeError(f"[process_accession] Invalid gzip file: {gz}")

        # 4) Classify
        classify_fastq_by_read_type(accession, os.path.join(workdir, "fastq_data"),
                                    debug_file, debug_lock)

        # 5) Cleanup invalid
        cleanup_invalid_fastqs(accession, os.path.join(workdir, "fastq_data"),
                               debug_file, debug_lock)

        # 6) Remove .sra
        remove_file_safely(os.path.join(workdir, "fastq_data", f"{accession}.sra"),
                           debug_file, debug_lock)

        # All success => append to completed
        append_line_with_lock(accession, completed_file, completed_lock)
        return f"[OK] {accession}"

    except Exception as e:
        # On any error => log debug, append to failed
        msg = f"[FAIL] {accession} => {str(e)}"
        log_debug_message(msg, debug_file, debug_lock)
        append_line_with_lock(accession, failed_file, failed_lock)
        return msg

##########################
# Main
##########################

def parse_args():
    parser = argparse.ArgumentParser(description="Download multiple SRA accessions in parallel.")
    parser.add_argument("--accessions-file", default="run_accessions.txt",
                        help="File of run accessions, one per line (default: run_accessions.txt)")
    parser.add_argument("--workdir", default=".",
                        help="Working directory for outputs/logs (default: current dir)")
    parser.add_argument("--debug-file", default="debug.log",
                        help="Path to a shared debug log file (default: debug.log)")
    parser.add_argument("--debug-lock", default="debug.lock",
                        help="Lock file for the debug log (default: debug.lock)")
    parser.add_argument("--completed-file", default="completed.txt",
                        help="File listing completed accessions (default: completed.txt)")
    parser.add_argument("--completed-lock", default="completed.lock",
                        help="Lock file for the completed-file (default: completed.lock)")
    parser.add_argument("--failed-file", default="failed.txt",
                        help="File listing failed accessions (default: failed.txt)")
    parser.add_argument("--failed-lock", default="failed.lock",
                        help="Lock file for the failed-file (default: failed.lock)")
    parser.add_argument("--combined-metadata", default="combined_metadata.tsv",
                        help="Path to the combined metadata TSV (default: combined_metadata.tsv)")
    parser.add_argument("--num-workers", type=int, default=4,
                        help="Number of parallel workers (default: 4)")
    parser.add_argument("--api-key", default=None,
                        help="NCBI API key (optional). If provided, exported as NCBI_API_KEY.")
    parser.add_argument("--email", default=None,
                        help="Email address for EDirect (optional). If provided, exported as EMAIL.")
    return parser.parse_args()

def cleanup_locks(args):
    """Remove lock files on exit."""
    lock_files = [args.completed_lock, args.debug_lock, args.failed_lock]
    for lock_file in lock_files:
        if os.path.exists(lock_file):
            try:
                os.remove(lock_file)
                # Since this runs on exit, we can't use log_debug_message reliably
                # (file descriptors might be closed), so print to stderr instead
                sys.stderr.write(f"[INFO] Removed lock file: {lock_file}\n")
            except Exception as e:
                sys.stderr.write(f"[WARN] Failed to remove lock file {lock_file}: {e}\n")

def main():
    args = parse_args()

    # Register cleanup function to run on exit
    atexit.register(cleanup_locks, args)

    # If the user provided an API key or email, export them for EDirect / sra-toolkit
    if args.api_key:
        os.environ["NCBI_API_KEY"] = args.api_key
    if args.email:
        os.environ["EMAIL"] = args.email

    # 1) Read all run accessions into a set (so duplicates are skipped)
    if not os.path.isfile(args.accessions_file):
        sys.stderr.write(f"[ERROR] --accessions-file not found: {args.accessions_file}\n")
        sys.exit(1)
    with open(args.accessions_file, "r") as f:
        all_accs = {line.strip() for line in f if line.strip()}

    # 2) Read 'completed' file into a set => skip these
    done_accs = set()
    if os.path.isfile(args.completed_file):
        with open(args.completed_file, "r") as cf:
            done_accs = {line.strip() for line in cf if line.strip()}

    # 3) Find accessions we still need
    to_download = all_accs - done_accs
    if not to_download:
        print("[INFO] All accessions are already completed. Nothing to do.")
        sys.exit(0)

    print(f"[INFO] Found {len(all_accs)} total, {len(done_accs)} completed, "
          f"so {len(to_download)} remaining.")

    # 4) Make sure directories exist
    os.makedirs(os.path.join(args.workdir, "fastq_data"), exist_ok=True)
    os.makedirs(os.path.join(args.workdir, "metadata"), exist_ok=True)
    os.makedirs(os.path.dirname(args.debug_file) or args.workdir, exist_ok=True)

    # 5) Concurrency: use a ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        future_map = {}
        # sorted() is optional, just to maintain consistent order
        for acc in sorted(to_download):
            fut = executor.submit(process_accession, acc, args)
            future_map[fut] = acc

        for fut in as_completed(future_map):
            acc = future_map[fut]
            try:
                res = fut.result()
                print(res)  # e.g. "[OK] SRR1234" or "[FAIL] SRR9999 => reason"
                log_debug_message(res, args.debug_file, args.debug_lock)
            except Exception as e:
                # Should be rare if process_accession handles exceptions
                msg = f"[FATAL] Worker error for {acc}: {e}"
                print(msg, file=sys.stderr)
                log_debug_message(msg, args.debug_file, args.debug_lock)

    print("[INFO] All done.")

if __name__ == "__main__":
    main()

