#!/usr/bin/env python3
import os
import re
import sys
import time
import glob
import gzip
import fcntl
import shutil
import psutil
import random
import atexit
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

##########################
# Locking & Logging Helpers
##########################

def flock_exclusive(fd):
    fcntl.flock(fd, fcntl.LOCK_EX)

def flock_release(fd):
    fcntl.flock(fd, fcntl.LOCK_UN)

def append_line_with_lock(line, out_file, lock_file):
    try:
        with open(lock_file, "a+") as lk:
            flock_exclusive(lk)
            with open(out_file, "a") as f:
                f.write(line + "\n")
            flock_release(lk)
    except Exception as e:
        sys.stderr.write(f"[ERROR] Could not append '{line}' to {out_file}: {e}\n")

def log_debug_message(msg, debug_file, debug_lock):
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
    return os.path.isfile(path) and os.path.getsize(path) >= min_size_bytes

def is_valid_gzip(path):
    try:
        with gzip.open(path, "rb") as f:
            f.read(1)
        return True
    except:
        return False

def remove_file_safely(path, debug_file, debug_lock):
    if os.path.exists(path):
        try:
            os.remove(path)
        except Exception as e:
            log_debug_message(f"[ERROR] remove_file_safely({path}): {e}", debug_file, debug_lock)

def get_read_lengths(fastq_gz):
    lengths = set()
    try:
        with gzip.open(fastq_gz, "rt") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    lengths.add(len(line.strip()))
                if i >= 400:
                    break
        return lengths
    except Exception as e:
        return None

##########################
# Command Runner with Retry
##########################

@retry(
    stop=stop_after_attempt(5),
    wait=wait_exponential(multiplier=1, min=2, max=10),
    retry=retry_if_exception_type(RuntimeError)
)
def run_command(cmd_list, err_msg, debug_file, debug_lock):
    try:
        result = subprocess.run(cmd_list, check=True, text=True, stderr=subprocess.PIPE)
        return result.stdout
    except subprocess.CalledProcessError as e:
        msg = f"{err_msg}:\n{e.stderr}"
        log_debug_message(msg, debug_file, debug_lock)
        raise RuntimeError(msg) from e

##########################
# FASTQ Steps
##########################

def compress_fastqs(accession, fastq_dir, debug_file, debug_lock):
    pattern = os.path.join(fastq_dir, f"{accession}*.fastq")
    fastq_files = glob.glob(pattern)
    compressor = shutil.which("pigz") or shutil.which("gzip")
    if not compressor:
        raise RuntimeError("No pigz or gzip found on system!")

    for fq in fastq_files:
        if not fq.endswith(".gz"):
            cmd = [compressor]
            if "pigz" in compressor:
                cmd.extend(["-p", "8"])
            cmd.append(fq)
            run_command(cmd, f"[compress_fastqs] {accession} compression failed on {fq}",
                        debug_file, debug_lock)

def cleanup_invalid_fastqs(accession, fastq_dir, debug_file, debug_lock):
    gz_files = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz"))
    for gz in gz_files:
        if not is_valid_fastq(gz) or not is_valid_gzip(gz):
            log_debug_message(f"[WARN] Removing invalid FASTQ: {gz}", debug_file, debug_lock)
            try:
                os.remove(gz)
            except Exception as e:
                log_debug_message(f"[ERROR] while removing {gz}: {e}", debug_file, debug_lock)

def classify_fastq_by_read_type(accession, fastq_dir, debug_file, debug_lock):
    vdb_dump = shutil.which("vdb-dump")
    if not vdb_dump:
        log_debug_message(f"[classify] vdb-dump not found; defaulting to biological for {accession}.",
                          debug_file, debug_lock)
        bio_dir = os.path.join(fastq_dir, "fastq_biologicaldata")
        os.makedirs(bio_dir, exist_ok=True)
        for gz in glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz")):
            dest_path = os.path.join(bio_dir, os.path.basename(gz))
            shutil.move(gz, dest_path)
            log_debug_message(f"[classify] Moved {gz} -> {dest_path} (default biological)", debug_file, debug_lock)
        return

    cmd_read_type = f"vdb-dump {accession} -C READ_TYPE 2>/dev/null | head -n 1"
    try:
        result = subprocess.run(cmd_read_type, shell=True, capture_output=True, text=True, check=True)
        line = result.stdout.strip()
        if "READ_TYPE:" not in line:
            log_debug_message(f"[classify] No READ_TYPE info for {accession}; defaulting to biological.", debug_file, debug_lock)
            bio_dir = os.path.join(fastq_dir, "fastq_biologicaldata")
            os.makedirs(bio_dir, exist_ok=True)
            for gz in glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz")):
                dest_path = os.path.join(bio_dir, os.path.basename(gz))
                shutil.move(gz, dest_path)
                log_debug_message(f"[classify] Moved {gz} -> {dest_path} (default biological)", debug_file, debug_lock)
            return

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
                gz_unsuffixed = f"{accession}.fastq.gz"
                src_unsuffixed = os.path.join(fastq_dir, gz_unsuffixed)
                if os.path.exists(src_unsuffixed):
                    src_path = src_unsuffixed
                    gz_filename = gz_unsuffixed

            if not src_path:
                raise FileNotFoundError(f"[classify] Missing FASTQ for read #{i+1} in {accession}")

            dest_dir = bio_dir if "BIOLOGICAL" in read_type else tech_dir
            dest_path = os.path.join(dest_dir, gz_filename)
            log_debug_message(f"[classify] Move {src_path} -> {dest_path}", debug_file, debug_lock)
            shutil.move(src_path, dest_path)

    except subprocess.CalledProcessError as e:
        log_debug_message(f"[classify] vdb-dump failed for {accession}:\n{e.stderr}; defaulting to biological.", debug_file, debug_lock)
        bio_dir = os.path.join(fastq_dir, "fastq_biologicaldata")
        os.makedirs(bio_dir, exist_ok=True)
        for gz in glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz")):
            dest_path = os.path.join(bio_dir, os.path.basename(gz))
            shutil.move(gz, dest_path)
            log_debug_message(f"[classify] Moved {gz} -> {dest_path} (default biological)", debug_file, debug_lock)
    except Exception as e:
        log_debug_message(f"[classify] Error for {accession}: {e}; defaulting to biological.", debug_file, debug_lock)
        bio_dir = os.path.join(fastq_dir, "fastq_biologicaldata")
        os.makedirs(bio_dir, exist_ok=True)
        for gz in glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz")):
            dest_path = os.path.join(bio_dir, os.path.basename(gz))
            shutil.move(gz, dest_path)
            log_debug_message(f"[classify] Moved {gz} -> {dest_path} (default biological)", debug_file, debug_lock)

##########################
# Metadata Fetch
##########################

def fetch_sra_metadata(accession, metadata_dir, combined_meta, debug_file, debug_lock):
    csv_path = os.path.join(metadata_dir, f"{accession}-run-info.csv")
    cmd = f'esearch -db sra -query "{accession}" | efetch -format runinfo > "{csv_path}"'
    ret = subprocess.call(cmd, shell=True)
    if ret != 0 or not os.path.isfile(csv_path) or os.path.getsize(csv_path) == 0:
        log_debug_message(f"[metadata] Failed to fetch runinfo for {accession} (ret={ret})",
                          debug_file, debug_lock)
        return

    tsv_path = csv_path.replace(".csv", ".tsv")
    try:
        with open(csv_path, "r") as inf, open(tsv_path, "w") as outf:
            for line in inf:
                outf.write(line.replace(",", "\t"))

        os.makedirs(os.path.dirname(combined_meta), exist_ok=True)
        with open(tsv_path, "r") as t_in:
            lines = t_in.readlines()
        with open(debug_lock, "a+") as lk:
            flock_exclusive(lk)
            if not os.path.isfile(combined_meta):
                with open(combined_meta, "w") as cm:
                    cm.write(lines[0])  # Header
            with open(combined_meta, "a") as cm:
                cm.writelines(lines[1:])  # Data rows
            flock_release(lk)

        remove_file_safely(csv_path, debug_file, debug_lock)
        remove_file_safely(tsv_path, debug_file, debug_lock)

    except Exception as e:
        log_debug_message(f"[metadata] Error merging runinfo for {accession}: {e}", debug_file, debug_lock)

##########################
# SRA Download Route with Custom tmp_dir
##########################

def sra_download_route(accession, workdir, debug_file, debug_lock, combined_meta, n_threads_per_worker, tmp_dir=None):
    fastq_dir = os.path.join(workdir, "fastq_data")
    metadata_dir = os.path.join(workdir, "metadata")
    tmp_fastq_dir = os.path.join(tmp_dir or fastq_dir, f"tmp/fastq_{accession}")
    os.makedirs(tmp_fastq_dir, exist_ok=True)

    mem_bytes = psutil.virtual_memory().available // n_threads_per_worker
    mem_gb = max(1, int(mem_bytes / (1024**3)))
    mem_str = f"{mem_gb}G"

    sra_file = os.path.join(tmp_fastq_dir, f"{accession}.sra")

    try:
        run_command(["prefetch", accession, "--max-size", "200G", "--temp-location", tmp_dir or fastq_dir, "--output-file", sra_file],
                    f"[sra_download_route] prefetch failed {accession}", debug_file, debug_lock)
    except RuntimeError:
        shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
        return False

    try:
        subprocess.run(["vdb-validate", sra_file], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        log_debug_message(f"[sra_download_route] vdb-validate failed {accession}:\n{e.stderr}",
                          debug_file, debug_lock)
        shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
        return False

    try:
        run_command([
            "fasterq-dump", sra_file, 
            "--temp", tmp_dir or tmp_fastq_dir,
            "--outdir", fastq_dir,
            "--threads", str(n_threads_per_worker), 
            "--mem", mem_str,
            "--split-files", "--include-technical"
        ], f"[sra_download_route] fasterq-dump failed {accession}", debug_file, debug_lock)
    except RuntimeError:
        shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
        return False

    fetch_sra_metadata(accession, metadata_dir, combined_meta, debug_file, debug_lock)

    fastqs = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq"))
    if not fastqs:
        log_debug_message(f"[sra_download_route] No FASTQs for {accession}", debug_file, debug_lock)
        shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
        return False
    for fq in fastqs:
        if not is_valid_fastq(fq, 1024):
            log_debug_message(f"[sra_download_route] Invalid or tiny FASTQ: {fq}", debug_file, debug_lock)
            shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
            return False

    compressor = shutil.which("pigz") or shutil.which("gzip")
    if not compressor:
        raise RuntimeError("No pigz or gzip found on system!")
    for fq in fastqs:
        cmd = [compressor, "-p", "8"] if "pigz" in compressor else [compressor]
        cmd.append(fq)
        run_command(cmd, f"[sra_download_route] Compression failed for {fq}", debug_file, debug_lock)

    shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
    return True

##########################
# GSA Download Route with Custom tmp_dir
##########################

def gsa_download_route(accession, workdir, debug_file, debug_lock):
    fastq_dir = os.path.join(workdir, "fastq_data")
    tmp_fastq_dir = os.path.join(fastq_dir, f"tmp/fastq_{accession}")
    os.makedirs(tmp_fastq_dir, exist_ok=True)

    iseq = shutil.which("iseq")
    if not iseq:
        log_debug_message(f"[gsa_download_route] 'iseq' not found in PATH for {accession}", debug_file, debug_lock)
        shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
        return False

    try:
        run_command([
            iseq, "-i", accession, "-o", tmp_fastq_dir, "-a"
        ], f"[gsa_download_route] iseq failed for {accession}", debug_file, debug_lock)
    except RuntimeError:
        shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
        return False

    fastqs = glob.glob(os.path.join(tmp_fastq_dir, f"{accession}*.fastq"))
    if not fastqs:
        log_debug_message(f"[gsa_download_route] No FASTQs found after iseq for {accession}", debug_file, debug_lock)
        shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
        return False
    for fq in fastqs:
        if not is_valid_fastq(fq):
            log_debug_message(f"[gsa_download_route] Invalid FASTQ found: {fq}", debug_file, debug_lock)
            shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
            return False

    compressor = shutil.which("pigz") or shutil.which("gzip")
    for fq in fastqs:
        cmd = [compressor, "-p", "8"] if "pigz" in compressor else [compressor]
        cmd.append(fq)
        run_command(cmd, f"[gsa_download_route] Compression failed for {fq}", debug_file, debug_lock)

    classify_fastq_by_read_type(accession, tmp_fastq_dir, debug_file, debug_lock)

    shutil.rmtree(tmp_fastq_dir, ignore_errors=True)
    return True

##########################
# Updated Worker Function
##########################

def process_accession(accession, args, n_threads_per_worker, tmp_dir=None, final_retry=False):
    debug_lock     = args.debug_lock
    debug_file     = args.debug_file
    completed_file = args.completed_file
    completed_lock = args.completed_lock
    failed_file    = args.failed_file
    failed_lock    = args.failed_lock
    workdir        = args.workdir
    combined_meta  = args.combined_metadata

    try:
        is_gsa = accession.startswith("C") or accession.startswith("GSA")
        ok = False
        if is_gsa:
            ok = gsa_download_route(accession, workdir, debug_file, debug_lock)
        else:
            ok = sra_download_route(accession, workdir, debug_file, debug_lock, combined_meta, n_threads_per_worker, tmp_dir)

        if not ok:
            raise RuntimeError(f"[process_accession] Download route failed for {accession}")

        fastq_dir = os.path.join(workdir, "fastq_data")
        cleanup_invalid_fastqs(accession, fastq_dir, debug_file, debug_lock)
        classify_fastq_by_read_type(accession, fastq_dir, debug_file, debug_lock)

        append_line_with_lock(accession, completed_file, completed_lock)
        return f"[OK] {accession}"

    except Exception as e:
        msg = f"[FAIL] {accession} => {str(e)}"
        log_debug_message(msg, debug_file, debug_lock)
        if not final_retry:
            append_line_with_lock(accession, failed_file, failed_lock)
        return msg

##########################
# Post-Processing Checks
##########################

def validate_accessions_file(accessions_file, debug_file, debug_lock):
    if not os.path.isfile(accessions_file):
        log_debug_message(f"[validate_accessions_file] File not found: {accessions_file}", debug_file, debug_lock)
        return set()
    
    accession_pattern = re.compile(r"^[ESDR]{1}RR\d+$")
    valid_accessions = set()
    with open(accessions_file, "r") as f:
        for i, line in enumerate(f, 1):
            acc = line.strip()
            if not acc:
                log_debug_message(f"[validate_accessions_file] Empty line at {i}", debug_file, debug_lock)
                continue
            if not accession_pattern.match(acc):
                log_debug_message(f"[validate_accessions_file] Invalid accession at line {i}: {acc}", debug_file, debug_lock)
                continue
            valid_accessions.add(acc)
    return valid_accessions

def check_and_sort_fastqs(args, n_threads_per_worker):
    fastq_dir = os.path.join(args.workdir, "fastq_data")
    bio_dir = os.path.join(fastq_dir, "fastq_biologicaldata")
    tech_dir = os.path.join(fastq_dir, "fastq_technicaldata")
    os.makedirs(bio_dir, exist_ok=True)
    os.makedirs(tech_dir, exist_ok=True)

    expected = validate_accessions_file(args.accessions_file, args.debug_file, args.debug_lock)
    if not expected:
        log_debug_message("[check_and_sort_fastqs] No valid accessions to check", args.debug_file, args.debug_lock)
        return

    bio_files = {re.match(r"(.+?)(?:_\d+)?\.fastq\.gz$", f.name).group(1) for f in Path(bio_dir).glob("*.fastq.gz") if re.match(r"(.+?)(?:_\d+)?\.fastq\.gz$", f.name)}
    missing = expected - bio_files

    if missing:
        log_debug_message(f"[check_and_sort_fastqs] Missing {len(missing)} accessions in biological data: {', '.join(sorted(missing))}", args.debug_file, args.debug_lock)

    for acc in missing:
        tmp_fastq_dir = os.path.join(args.tmp_dir or fastq_dir, f"tmp/fastq_{acc}")
        gz_files = glob.glob(os.path.join(fastq_dir, f"{acc}*.fastq.gz"))
        sra_file = os.path.join(fastq_dir, f"{acc}.sra")

        if gz_files:
            all_valid = True
            for gz in gz_files:
                if not (is_valid_fastq(gz) and is_valid_gzip(gz)):
                    all_valid = False
                    log_debug_message(f"[check_and_sort_fastqs] Incomplete/invalid file: {gz}", args.debug_file, args.debug_lock)
                    remove_file_safely(gz, args.debug_file, args.debug_lock)

            if all_valid:
                lengths = set()
                for gz in gz_files:
                    lens = get_read_lengths(gz)
                    if lens:
                        lengths.update(lens)
                    else:
                        all_valid = False
                        log_debug_message(f"[check_and_sort_fastqs] Could not read lengths from {gz}", args.debug_file, args.debug_lock)
                        break

                dest_dir = bio_dir if len(lengths) > 1 else tech_dir
                for gz in gz_files:
                    shutil.move(gz, os.path.join(dest_dir, os.path.basename(gz)))
                    log_debug_message(f"[check_and_sort_fastqs] Moved {gz} to {dest_dir}", args.debug_file, args.debug_lock)
                continue

        if os.path.exists(sra_file):
            try:
                subprocess.run(["vdb-validate", sra_file], check=True, capture_output=True, text=True)
                log_debug_message(f"[check_and_sort_fastqs] Validated {sra_file}", args.debug_file, args.debug_lock)

                for pattern in [f"{acc}*.fastq", f"{acc}*.fastq.gz"]:
                    for f in glob.glob(os.path.join(fastq_dir, pattern)):
                        remove_file_safely(f, args.debug_file, args.debug_lock)

                mem_bytes = psutil.virtual_memory().available
                mem_str = f"{int(mem_bytes / (1024**3))}G"
                run_command([
                    "fasterq-dump", sra_file, 
                    "--temp", args.tmp_dir or tmp_fastq_dir,
                    "--outdir", fastq_dir,
                    "--threads", str(n_threads_per_worker), 
                    "--mem", mem_str,
                    "--split-files", 
                    "--include-technical"
                ], f"[check_and_sort_fastqs] fasterq-dump failed for {acc}", args.debug_file, args.debug_lock)

                compress_fastqs(acc, fastq_dir, args.debug_file, args.debug_lock)

                gz_files = glob.glob(os.path.join(fastq_dir, f"{acc}*.fastq.gz"))
                if gz_files and all(is_valid_gzip(gz) for gz in gz_files):
                    lengths = set()
                    for gz in gz_files:
                        lens = get_read_lengths(gz)
                        if lens:
                            lengths.update(lens)
                    dest_dir = bio_dir if len(lengths) > 1 else tech_dir
                    for gz in gz_files:
                        shutil.move(gz, os.path.join(dest_dir, os.path.basename(gz)))
                        log_debug_message(f"[check_and_sort_fastqs] Moved {gz} to {dest_dir}", args.debug_file, args.debug_lock)
                    remove_file_safely(sra_file, args.debug_file, args.debug_lock)
                else:
                    append_line_with_lock(acc, args.failed_file, args.failed_lock)
            except Exception as e:
                log_debug_message(f"[check_and_sort_fastqs] Failed to reprocess {acc}: {e}", args.debug_file, args.debug_lock)
                append_line_with_lock(acc, args.failed_file, args.failed_lock)
        else:
            log_debug_message(f"[check_and_sort_fastqs] No .sra or valid .fastq.gz for {acc}", args.debug_file, args.debug_lock)
            append_line_with_lock(acc, args.failed_file, args.failed_lock)

def final_validation(args, n_threads_per_worker):
    expected = validate_accessions_file(args.accessions_file, args.debug_file, args.debug_lock)
    bio_dir = os.path.join(args.workdir, "fastq_data", "fastq_biologicaldata")
    fastq_dir = os.path.join(args.workdir, "fastq_data")

    bio_files = {
        re.match(r"(.+?)(?:_\d+)?\.fastq\.gz$", f.name).group(1)
        for f in Path(bio_dir).glob("*.fastq.gz")
        if re.match(r"(.+?)(?:_\d+)?\.fastq\.gz$", f.name)
    }
    missing = expected - bio_files

    retry_these = []
    for acc in missing:
        sra_file = os.path.join(fastq_dir, f"{acc}.sra")
        fastqs = glob.glob(os.path.join(fastq_dir, f"{acc}*.fastq.gz"))
        if not fastqs and not os.path.exists(sra_file):
            log_debug_message(f"[final_validation] {acc} completely missing. Will retry.", 
                              args.debug_file, args.debug_lock)
            retry_these.append(acc)

    if not missing:
        log_debug_message("[final_validation] All accessions accounted for in biological data", 
                          args.debug_file, args.debug_lock)
        return
    else:
        log_debug_message(f"[final_validation] {len(missing)} accessions are missing. "
                          f"Retrying {len(retry_these)} of them.", 
                          args.debug_file, args.debug_lock)

    if not retry_these:
        return

    if os.path.exists(args.failed_file):
        with open(args.failed_file, "r") as f:
            old_failures = {line.strip() for line in f if line.strip()}
        updated_failures = old_failures - set(retry_these)
        with open(args.failed_file, "w") as f:
            for acc_fail in sorted(updated_failures):
                f.write(acc_fail + "\n")

    with ProcessPoolExecutor(max_workers=min(os.cpu_count() or 4, len(retry_these))) as executor:
        future_map = {}
        for acc in retry_these:
            time.sleep(random.uniform(5, 30))
            fut = executor.submit(process_accession, acc, args, n_threads_per_worker, args.tmp_dir, True)
            future_map[fut] = acc

        for fut in as_completed(future_map):
            acc = future_map[fut]
            try:
                result = fut.result()
                if "[FAIL]" in result:
                    append_line_with_lock(acc, args.failed_file, args.failed_lock)
                log_debug_message(f"[final_validation] Retry result for {acc}: {result}", 
                                  args.debug_file, args.debug_lock)
            except Exception as e:
                msg = f"[final_validation] Retry crashed for {acc}: {e}"
                log_debug_message(msg, args.debug_file, args.debug_lock)
                append_line_with_lock(acc, args.failed_file, args.failed_lock)

def write_sorted_unique_completed(completed_file):
    if os.path.isfile(completed_file):
        with open(completed_file, "r") as f:
            entries = sorted(set(line.strip() for line in f if line.strip()))
        with open(completed_file, "w") as f:
            for line in entries:
                f.write(line + "\n")

##########################
# Main with Dynamic Workers
##########################

def parse_args():
    parser = argparse.ArgumentParser(description="Download multiple SRA accessions in parallel with dynamic workers.")
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
    parser.add_argument("--num-workers", type=int, default=None,
                        help="Number of parallel workers (default: dynamic, based on CPU and tasks, max 64)")
    parser.add_argument("--api-key", default=None,
                        help="NCBI API key (optional). If provided, exported as NCBI_API_KEY.")
    parser.add_argument("--email", default=None,
                        help="Email address for EDirect (optional). If provided, exported as EMAIL.")
    parser.add_argument("--tmp-dir", default=None,
                        help="Directory for temporary files (default: within workdir/fastq_data)")
    return parser.parse_args()

def cleanup_locks(args):
    lock_files = [args.completed_lock, args.debug_lock, args.failed_lock]
    for lock_file in lock_files:
        if os.path.exists(lock_file):
            try:
                os.remove(lock_file)
                sys.stderr.write(f"[INFO] Removed lock file: {lock_file}\n")
            except Exception as e:
                sys.stderr.write(f"[WARN] Failed to remove lock file {lock_file}: {e}\n")

def main():
    args = parse_args()
    atexit.register(cleanup_locks, args)

    if args.api_key:
        os.environ["NCBI_API_KEY"] = args.api_key
        log_debug_message("[main] Set NCBI_API_KEY environment variable.", args.debug_file, args.debug_lock)

    if args.email:
        os.environ["EMAIL"] = args.email
        log_debug_message(f"[main] Set EMAIL environment variable to {args.email}.", args.debug_file, args.debug_lock)

    log_debug_message("[main] Validating accessions file...", args.debug_file, args.debug_lock)
    all_accs = validate_accessions_file(args.accessions_file, args.debug_file, args.debug_lock)
    if not all_accs:
        msg = f"[ERROR] No valid accessions in {args.accessions_file}"
        sys.stderr.write(msg + "\n")
        log_debug_message(msg, args.debug_file, args.debug_lock)
        sys.exit(1)

    log_debug_message(f"[main] Found {len(all_accs)} valid accessions in {args.accessions_file}.", 
                      args.debug_file, args.debug_lock)

    done_accs = set()
    if os.path.isfile(args.completed_file):
        with open(args.completed_file, "r") as cf:
            done_accs = {line.strip() for line in cf if line.strip()}
    log_debug_message(f"[main] Found {len(done_accs)} previously completed accessions.", 
                      args.debug_file, args.debug_lock)

    to_download = all_accs - done_accs
    log_debug_message(f"[main] {len(to_download)} accessions remain to be processed.", 
                      args.debug_file, args.debug_lock)

    if not to_download:
        print("[INFO] All accessions are already completed. Nothing to do.")
        log_debug_message("[main] All accessions are already completed. Exiting early.", 
                          args.debug_file, args.debug_lock)
    else:
        print(f"[INFO] Found {len(all_accs)} total, {len(done_accs)} completed, "
              f"so {len(to_download)} remaining.")
        
        total_cores = os.cpu_count() or 4
        min_threads_per_worker = 4
        max_workers_cap = 64
        if args.num_workers is None:
            max_workers_by_cores = total_cores // min_threads_per_worker
            num_workers = min(max_workers_by_cores, len(to_download), max_workers_cap)
        else:
            num_workers = min(args.num_workers, max_workers_cap)
        n_threads_per_worker = max(min_threads_per_worker, total_cores // num_workers)
        log_debug_message(f"[main] Dynamically allocated {num_workers} workers with {n_threads_per_worker} threads each "
                          f"(total cores: {total_cores}, tasks: {len(to_download)}).",
                          args.debug_file, args.debug_lock)

        log_debug_message(f"[main] Launching parallel download for {len(to_download)} accessions "
                          f"using {num_workers} workers.", args.debug_file, args.debug_lock)

        os.makedirs(os.path.join(args.workdir, "fastq_data"), exist_ok=True)
        os.makedirs(os.path.join(args.workdir, "metadata"), exist_ok=True)
        os.makedirs(os.path.dirname(args.debug_file) or args.workdir, exist_ok=True)
        if args.tmp_dir:
            os.makedirs(args.tmp_dir, exist_ok=True)

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            future_map = {}
            for acc in sorted(to_download):
                time.sleep(random.uniform(0.1, 1.0))
                fut = executor.submit(process_accession, acc, args, n_threads_per_worker, args.tmp_dir)
                future_map[fut] = acc

            for fut in as_completed(future_map):
                acc = future_map[fut]
                try:
                    res = fut.result()
                    print(res)
                    log_debug_message(res, args.debug_file, args.debug_lock)
                except Exception as e:
                    msg = f"[FATAL] Worker error for {acc}: {e}"
                    print(msg, file=sys.stderr)
                    log_debug_message(msg, args.debug_file, args.debug_lock)
                    append_line_with_lock(acc, args.failed_file, args.failed_lock)

    log_debug_message("[main] Post-processing: Checking and sorting FASTQs...", args.debug_file, args.debug_lock)
    check_and_sort_fastqs(args, n_threads_per_worker)

    log_debug_message("[main] Post-processing: Final validation of downloaded files...", args.debug_file, args.debug_lock)
    final_validation(args, n_threads_per_worker)

    log_debug_message("[main] Writing sorted list of unique completed accessions...", args.debug_file, args.debug_lock)
    write_sorted_unique_completed(args.completed_file)

    log_debug_message("[main] Pipeline finished successfully.", args.debug_file, args.debug_lock)
    print("[INFO] All done.")

if __name__ == "__main__":
    main()
