#!/usr/bin/env python3
#SBATCH --export=ALL
#SBATCH --job-name=fastq_{ACCESSION}
#SBATCH --output={WORKDIR}/logs/{ACCESSION}.out
#SBATCH --error={WORKDIR}/logs/{ACCESSION}.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

import os
import sys
import subprocess
import time
import fcntl
import glob
import atexit
import argparse
import shutil
import xml.etree.ElementTree as ET

##########################
# Helper Functions
##########################

def check_required_tools(tools):
    missing = []
    for tool in tools:
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        sys.stderr.write(f"ERROR: Missing required tools in PATH: {', '.join(missing)}\n")
        sys.exit(1)

def parse_cli_args():
    parser = argparse.ArgumentParser(
        description="Download a single ENA/SRA accession via Slurm, with fallback and metadata merging."
    )
    parser.add_argument('--accession', help='Accession to download (e.g., ERR123456)')
    parser.add_argument('--workdir', help='Working directory')
    parser.add_argument('--dry-run', action='store_true', help='Only print actions without running them')
    parser.add_argument('--metadata', help='Combined Metadata file (e.g., combined_metadata.csv)') 
    parser.add_argument('--debug', action='store_true', help='Enable extra debug output')
    return parser.parse_args()

def read_env_var(name, required=True):
    """Reads an environment variable or exits if missing."""
    val = os.environ.get(name)
    if required and (val is None or val.strip() == ""):
        sys.stderr.write(f"ERROR: Missing required env var: {name}\n")
        sys.exit(1)
    return val

def flock_exclusive(fd):
    """Acquire an exclusive lock on a file descriptor."""
    fcntl.flock(fd, fcntl.LOCK_EX)

def flock_release(fd):
    """Release the lock on a file descriptor."""
    fcntl.flock(fd, fcntl.LOCK_UN)

def run_command(cmd_list, err_msg):
    """
    Runs a shell command with subprocess.run(check=True).
    Raises RuntimeError if it fails.
    """
    try:
        subprocess.run(cmd_list, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"{err_msg}: {e}")

def log_debug_message(debug_lock_path, message):
    """
    Appends a message to debug.log under a file lock.
    Reads the WORKDIR from the environment.
    """
    workdir = read_env_var("WORKDIR")
    debug_log_path = os.path.join(workdir, "logs", "debug.log")
    try:
        # Open the lock file (could be the same as the log or separate)
        with open(debug_lock_path, "a+") as lk:
            # (Optional: acquire lock with flock_exclusive(lk) if needed)
            lk.write(f"{time.ctime()} {message}\n")
            # (Optional: release lock with flock_release(lk))
    except Exception as e:
        sys.stderr.write(f"Could not log debug message '{message}' due to: {e}\n")

def remove_file_safely(path, accession, debug_lock_path):
    """Remove a file and log an error if it fails."""
    if os.path.exists(path):
        try:
            os.remove(path)
        except Exception as e:
            log_debug_message(debug_lock_path, f"❌ ERROR: {accession} Failed to delete {path}: {e}")
            raise

def append_to_checkpoint(accession, checkpoint_file, checkpoint_lock):
    """Atomically appends the accession to the checkpoint file."""
    try:
        with open(checkpoint_lock, "r+") as cplock:
            flock_exclusive(cplock)
            with open(checkpoint_file, "a") as cpf:
                cpf.write(accession + "\n")
            flock_release(cplock)
    except Exception as e:
        sys.stderr.write(f"Failed to write to checkpoint for {accession}: {e}\n")

def release_token():
    """
    Mimics the bash 'release_token' function:
    - Locks token file
    - Increments the token count
    - Logs the release using the current ACCESSION from env.
    """
    accession = read_env_var("ACCESSION")
    token_file = read_env_var("TOKEN_FILE")
    tokenlock_file = read_env_var("TOKEN_LOCK_FILE")
    workdir = read_env_var("WORKDIR")

    try:
        with open(tokenlock_file, "r+") as lf:
            flock_exclusive(lf)
            with open(token_file, "r") as tf:
                current_tokens = int(tf.read().strip())
            new_tokens = current_tokens + 1
            with open(token_file, "w") as tf:
                tf.write(str(new_tokens) + "\n")
            # Log the token release
            with open(os.path.join(workdir, "token_audit.log"), "a") as ta:
                ta.write(f"[{time.ctime()}] RELEASED TOKEN FOR {accession} (NOW {new_tokens})\n")
    except Exception as e:
        with open(os.path.join(workdir, "lock_errors.log"), "a") as errf:
            errf.write(f"[{time.ctime()}] FAILED LOCK FOR {accession}: {str(e)}\n")

def compress_fastqs(fastq_dir, accession, debug_lock_path):
    """Find all uncompressed *.fastq files in fastq_dir matching accession and gzip them."""
    pattern = os.path.join(fastq_dir, f"{accession}*.fastq")
    fastq_files = glob.glob(pattern)
    for fq in fastq_files:
        if not fq.endswith(".gz"):
            try:
                run_command(["gzip", fq], f"ERROR: {accession} gzip failed on {fq}")
            except RuntimeError as e:
                log_debug_message(debug_lock_path, str(e))

def guess_ena_ftp_paths(accession):
    """
    Return a list of possible ftp paths for a given ENA accession.
    """
    base_ftp = "ftp.sra.ebi.ac.uk/vol1/fastq"
    prefix = accession[:6]  # e.g. ERR119
    suffix = accession[6:]   # e.g. 3450

    paths = [
        f"{base_ftp}/{prefix}/{suffix}/{accession}",
        f"{base_ftp}/{prefix}/00{suffix[-1]}/{accession}",
        f"{base_ftp}/{prefix}/{suffix[-3:]}/{accession}",
        f"{base_ftp}/{prefix}/{accession}",
    ]
    unique_paths = list(dict.fromkeys(paths))
    return unique_paths

def flatten_xml_to_dict(element, path=""):
    """
    Recursively flatten an XML tree to {path: value} pairs.
    """
    flattened = {}
    tag_path = f"{path}/{element.tag}" if path else element.tag

    if len(element) == 0:
        text = element.text.strip() if element.text else "NA"
        flattened[tag_path] = text
    else:
        for child in element:
            flattened.update(flatten_xml_to_dict(child, tag_path))
    return flattened

def try_enaDataGet(accession, fastq_dir, metadata_dir, debug_lock_path):
    """
    Attempt to download using enaDataGet. If successful:
    - Convert XML to flattened TSV
    - Append metadata to project-level TSV file
    - Move FASTQ files into fastq_dir
    - Remove the accession subdirectory
    Return True if successful, False otherwise.
    """
    cmd = ["enaDataGet", "-f", "fastq", "-d", fastq_dir, "-m", "True", accession]
    accession_dir = os.path.join(fastq_dir, accession)
    xml_path = os.path.join(accession_dir, f"{accession}.xml")

    try:
        run_command(cmd, f"❌ ERROR: {accession} enaDataGet failed")

        if not os.path.exists(xml_path):
            log_debug_message(debug_lock_path, f"⚠️ WARNING: XML metadata file missing for {accession}")
            return False

        # Parse XML and flatten to dictionary
        tree = ET.parse(xml_path)
        root = tree.getroot()
        flat_data = flatten_xml_to_dict(root)

        # Extract study/project accession
        study_accession = root.findtext(".//STUDY/IDENTIFIERS/PRIMARY_ID", default=None)
        if not study_accession:
            study_accession = root.findtext(".//STUDY_REF/IDENTIFIERS/PRIMARY_ID", default="unknown_project")
        if not study_accession:
            study_accession = "unknown_project"

        project_tsv = os.path.join(metadata_dir, f"{study_accession}-project-metadata.tsv")
        os.makedirs(metadata_dir, exist_ok=True)

        # Write header if needed and append TSV row
        write_header = not os.path.exists(project_tsv)
        with open(project_tsv, "a") as out:
            if write_header:
                out.write("\t".join(flat_data.keys()) + "\n")
            out.write("\t".join(flat_data.values()) + "\n")

        print(f"📄 Appended metadata to project file: {project_tsv}")

        # Move FASTQ files to fastq_dir
        fastqs = glob.glob(os.path.join(accession_dir, "*.fastq.gz"))
        for fq in fastqs:
            dest = os.path.join(fastq_dir, os.path.basename(fq))
            shutil.move(fq, dest)
            print(f"📦 Moved FASTQ: {fq} → {dest}")

        # Delete the accession subfolder
        shutil.rmtree(accession_dir)
        print(f"🗑️ Deleted {accession_dir}")

        return True

    except Exception as e:
        log_debug_message(debug_lock_path, f"❌ ERROR in try_enaDataGet for {accession}: {e}")
        return False

def try_direct_wget(accession, fastq_dir, debug_lock_path, max_retries=2):
    """
    Attempt direct wget from guessed ENA ftp paths.
    Return True if any file is successfully downloaded, else False.
    """
    success = False
    possible_dirs = guess_ena_ftp_paths(accession)
    candidates = [
        f"{accession}.fastq.gz",
        f"{accession}_1.fastq.gz",
        f"{accession}_2.fastq.gz",
    ]
    for dir_url in possible_dirs:
        for f in candidates:
            url = f"https://{dir_url}/{f}"
            for attempt in range(1, max_retries+1):
                cmd = ["wget", "-O", os.path.join(fastq_dir, f), url]
                try:
                    run_command(cmd, f"Wget attempt {attempt} for {accession} failed: {url}")
                    success = True
                    break  # exit retry loop if successful
                except RuntimeError as e:
                    log_debug_message(debug_lock_path, str(e))
                    time.sleep(2 * attempt)
            # Continue to try other candidate files even if one succeeds
    downloaded = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz"))
    return len(downloaded) > 0

def try_ftp_fallback(accession, fastq_dir, debug_lock_path):
    """
    Fallback method using ftp URLs with 'curl' or 'wget'.
    Return True if any file is downloaded, else False.
    """
    success = False
    possible_dirs = guess_ena_ftp_paths(accession)
    candidates = [
        f"{accession}.fastq.gz",
        f"{accession}_1.fastq.gz",
        f"{accession}_2.fastq.gz",
    ]
    for dir_path in possible_dirs:
        for f in candidates:
            url = f"ftp://{dir_path}/{f}"
            out_path = os.path.join(fastq_dir, f)
            cmd = ["curl", "-o", out_path, url]
            try:
                run_command(cmd, f"curl ftp fallback for {accession} failed: {url}")
                success = True
            except RuntimeError as e:
                log_debug_message(debug_lock_path, str(e))
                cmd2 = ["wget", "-O", out_path, url]
                try:
                    run_command(cmd2, f"wget ftp fallback for {accession} failed: {url}")
                    success = True
                except RuntimeError as e2:
                    log_debug_message(debug_lock_path, str(e2))
    downloaded = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz"))
    return len(downloaded) > 0

def sra_route(accession, fastq_dir, metadata_dir, combined_meta, debug_lock_path, checkpoint_lock):
    """
    Attempt SRA route: prefetch + fasterq-dump + SRA metadata.
    Return True if successful, else False.
    """
    sra_file = os.path.join(fastq_dir, f"{accession}.sra")

    print(f"🔹 [Fallback] Prefetching SRA file for {accession}")
    cmd_prefetch = ["prefetch", accession, "--output-file", sra_file]
    try:
        run_command(cmd_prefetch, f"❌ ERROR: {accession} prefetch failed")
    except RuntimeError as e:
        log_debug_message(debug_lock_path, str(e))
        return False
    time.sleep(2)

    print(f"🔹 [Fallback] Converting SRA to FASTQ for {accession}")
    cmd_fasterq = [
        "fasterq-dump",
        sra_file,
        "--outdir", fastq_dir,
        "--threads", "4",
        "--mem", "8G",
        "--split-3"
    ]
    try:
        run_command(cmd_fasterq, f"❌ ERROR: {accession} fasterq-dump failed")
    except RuntimeError as e:
        log_debug_message(debug_lock_path, str(e))
        return False
    time.sleep(2)

    # Fetch SRA metadata via esearch/efetch
    max_retries  = 3
    success_meta = False
    csv_path = os.path.join(metadata_dir, f"{accession}-run-info.csv")

    for i in range(1, max_retries + 1):
        cmd_esearch = f'esearch -db sra -query "{accession}" | efetch -format runinfo > {csv_path}'
        ret = subprocess.call(cmd_esearch, shell=True)
        if ret == 0 and os.path.isfile(csv_path) and os.path.getsize(csv_path) > 0:
            success_meta = True
            tsv_file = os.path.join(metadata_dir, f"{accession}-run-info.tsv")
            try:
                with open(csv_path, "r") as inf, open(tsv_file, "w") as outf:
                    for line in inf:
                        outf.write(line.replace(",", "\t"))
                with open(checkpoint_lock, "r+") as lock_fd:
                    flock_exclusive(lock_fd)
                    if not os.path.isfile(combined_meta):
                        with open(tsv_file, "r") as tsv_in:
                            header_line = tsv_in.readline()
                        with open(combined_meta, "w") as cm:
                            cm.write(header_line)
                    with open(tsv_file, "r") as tsv_in, open(combined_meta, "a") as cm:
                        _ = tsv_in.readline()  # skip header
                        for row in tsv_in:
                            cm.write(row)
                    flock_release(lock_fd)
                remove_file_safely(tsv_file, accession, debug_lock_path)
                remove_file_safely(csv_path, accession, debug_lock_path)
            except Exception as e:
                log_debug_message(debug_lock_path, f"⚠️ WARNING: {accession} SRA metadata handling error: {e}")
            break
        else:
            log_debug_message(debug_lock_path, f"⚠️ WARNING: {accession} SRA metadata attempt {i} failed")
            time.sleep(i)

    if not success_meta:
        log_debug_message(debug_lock_path, f"❌ ERROR: {accession} All SRA metadata attempts failed")

    all_fastqs = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq"))
    if not all_fastqs:
        log_debug_message(debug_lock_path, f"❌ ERROR: {accession} No .fastq from fallback SRA route.")
        return False

    return True

##########################
# Main entry point
##########################
if __name__ == "__main__":

    check_required_tools([
        "enaDataGet", "prefetch", "fasterq-dump",
        "esearch", "efetch", "wget", "curl", "gzip"
    ])
    
    args = parse_cli_args()

    # Allow CLI args to override environment variables (useful for testing)
    if args.accession:
        os.environ["ACCESSION"] = args.accession
    if args.workdir:
        os.environ["WORKDIR"] = args.workdir
    if args.metadata:
        os.environ["COMBINED_METADATA"] = args.metadata

    # Required variables from environment
    accession       = read_env_var("ACCESSION")
    workdir         = read_env_var("WORKDIR")
    checkpoint_file = read_env_var("CHECKPOINT_FILE")
    checkpoint_lock = read_env_var("CHECKPOINT_LOCK_FILE")
    combined_meta   = read_env_var("COMBINED_METADATA")
    debug_lock      = read_env_var("DEBUG_LOCK")
    token_file      = read_env_var("TOKEN_FILE")
    tokenlock_file  = read_env_var("TOKEN_LOCK_FILE")

    fastq_dir    = os.path.join(workdir, "fastq_data")
    metadata_dir = os.path.join(workdir, "metadata")
    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(metadata_dir, exist_ok=True)

    @atexit.register
    def on_exit():
        release_token()

    print(f"[DEBUG] Accession={accession}", flush=True)

    # Determine provider by accession prefix
    prefix = accession[:3].upper()
    if prefix in ("SRR", "SRX", "SRS", "SRP"):
        provider = "sra"
    elif prefix in ("ERR", "ERX", "ERS", "ERP", "DRR", "DRX", "DRS", "DRP"):
        provider = "ena"
    else:
        log_debug_message(debug_lock, f"❌ ERROR: {accession} Unrecognized prefix => no known route.")
        sys.exit(1)

    if provider == "sra":
        print(f"[INFO] Using SRA route for {accession}")
        ok = sra_route(accession, fastq_dir, metadata_dir, combined_meta, debug_lock, checkpoint_lock)
        if not ok:
            sys.exit(1)
    else:
        print(f"[INFO] Using ENA route(s) for {accession}")
        # Try enaDataGet first (pass metadata_dir as well)
        ena_ok = try_enaDataGet(accession, fastq_dir, metadata_dir, debug_lock)
        if not ena_ok:
            print(f"[WARN] enaDataGet failed => trying direct wget for {accession}")
            ena_ok = try_direct_wget(accession, fastq_dir, debug_lock, max_retries=2)
            if not ena_ok:
                print(f"[WARN] direct wget failed => trying ftp fallback for {accession}")
                ena_ok = try_ftp_fallback(accession, fastq_dir, debug_lock)
                if not ena_ok:
                    print(f"[WARN] All ENA methods failed => fallback to SRA route for {accession}")
                    sra_ok = sra_route(accession, fastq_dir, metadata_dir, combined_meta, debug_lock, checkpoint_lock)
                    if not sra_ok:
                        log_debug_message(debug_lock, f"❌ ERROR: {accession} ENA + SRA fallback all failed.")
                        sys.exit(1)
                else:
                    print(f"[INFO] ftp fallback succeeded for {accession}")
            else:
                print(f"[INFO] direct wget succeeded for {accession}")
        else:
            print(f"[INFO] enaDataGet succeeded for {accession}")

    # Compress any leftover .fastq files
    compress_fastqs(fastq_dir, accession, debug_lock)

    # Check for .fastq.gz files
    gz_pattern = os.path.join(fastq_dir, f"{accession}*.fastq.gz")
    found_gz = glob.glob(gz_pattern)
    if not found_gz:
        log_debug_message(debug_lock, f"❌ ERROR: {accession} No FASTQ files found after process.")
        sys.exit(1)

    # Append accession to checkpoint file
    append_to_checkpoint(accession, checkpoint_file, checkpoint_lock)

    # Cleanup any leftover .sra file
    sra_file_path = os.path.join(fastq_dir, f"{accession}.sra")
    remove_file_safely(sra_file_path, accession, debug_lock)

    print(f"✅ Successfully processed {accession}")
