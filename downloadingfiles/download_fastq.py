#!/usr/bin/env python3
import os
import sys
import subprocess
import time
import fcntl
import glob
import atexit
import argparse
import shutil
import psutil
import xml.etree.ElementTree as ET
import gzip
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

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
    val = os.environ.get(name)
    if required and (val is None or val.strip() == ""):
        sys.stderr.write(f"ERROR: Missing required env var: {name}\n")
        sys.exit(1)
    return val

def flock_exclusive(fd):
    fcntl.flock(fd, fcntl.LOCK_EX)

def flock_release(fd):
    fcntl.flock(fd, fcntl.LOCK_UN)

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10), retry=retry_if_exception_type(subprocess.CalledProcessError))
def run_command(cmd_list, err_msg, debug_lock_path, accession):
    try:
        subprocess.run(cmd_list, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        error_msg = f"{err_msg}: {e.stderr}"
        log_debug_message(debug_lock_path, error_msg)
        raise RuntimeError(error_msg) from e

def log_debug_message(debug_lock_path, message):
    workdir = read_env_var("WORKDIR")
    debug_log_path = os.path.join(workdir, "logs", "debug.log")
    try:
        with open(debug_lock_path, "a") as lk:
            flock_exclusive(lk)
            with open(debug_log_path, "a") as logf:
                logf.write(f"{time.ctime()} {message}\n")
            flock_release(lk)
    except Exception as e:
        sys.stderr.write(f"Could not log debug message '{message}' due to: {e}\n")

def remove_file_safely(path, accession, debug_lock_path):
    if os.path.exists(path):
        try:
            os.remove(path)
        except Exception as e:
            log_debug_message(debug_lock_path, f"ERROR: {accession} Failed to delete {path}: {e}")
            raise

def append_to_checkpoint(accession, checkpoint_file, checkpoint_lock):
    try:
        with open(checkpoint_lock, "a+") as cplock:
            flock_exclusive(cplock)
            with open(checkpoint_file, "a") as cpf:
                cpf.write(accession + "\n")
            flock_release(cplock)
    except Exception as e:
        sys.stderr.write(f"Failed to write to checkpoint for {accession}: {e}\n")
        raise

def release_token():
    accession = read_env_var("ACCESSION")
    token_file = read_env_var("TOKEN_FILE")
    tokenlock_file = read_env_var("TOKEN_LOCK_FILE")
    workdir = read_env_var("WORKDIR")

    try:
        with open(tokenlock_file, "a+") as lf:
            flock_exclusive(lf)
            try:
                with open(token_file, "r") as tf:
                    contents = tf.read().strip()
                    current_tokens = int(contents) if contents.isdigit() else 0
            except Exception as read_err:
                current_tokens = 0
                log_debug_message(tokenlock_file, f"WARNING: Failed to read token file: {read_err}")
            new_tokens = current_tokens + 1
            with open(token_file, "w") as tf:
                tf.write(f"{new_tokens}\n")
            with open(os.path.join(workdir, "logs", "token_audit.log"), "a") as ta:
                ta.write(f"[{time.ctime()}] RELEASED TOKEN FOR {accession} (NOW {new_tokens})\n")
            flock_release(lf)
    except Exception as e:
        with open(os.path.join(workdir, "logs", "lock_errors.log"), "a") as errf:
            errf.write(f"[{time.ctime()}] FAILED LOCK FOR {accession}: {str(e)}\n")

def compress_fastqs(fastq_dir, accession, debug_lock_path):
    pattern = os.path.join(fastq_dir, f"{accession}*.fastq")
    fastq_files = glob.glob(pattern)
    for fq in fastq_files:
        if not fq.endswith(".gz"):
            try:
                threads = max(1, os.cpu_count() - 1)
                run_command(["pigz", "-p", str(threads), fq], f"ERROR: {accession} pigz failed on {fq}", debug_lock_path, accession)
            except RuntimeError as e:
                log_debug_message(debug_lock_path, str(e))

def flatten_xml_to_dict(element, path=""):
    flattened = {}
    tag_path = f"{path}/{element.tag}" if path else element.tag
    if len(element) == 0:
        text = element.text.strip() if element.text else "NA"
        flattened[tag_path] = text
    else:
        for child in element:
            flattened.update(flatten_xml_to_dict(child, tag_path))
    return flattened

def is_valid_fastq(path, min_size_bytes=1024):
    return os.path.isfile(path) and os.path.getsize(path) >= min_size_bytes

def is_valid_gzip(path):
    try:
        with gzip.open(path, 'rb') as f:
            f.read(1)
        return True
    except:
        return False

def cleanup_invalid_fastqs(fastq_dir, accession):
    removed = False
    for f in glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz")):
        if not is_valid_fastq(f) or not is_valid_gzip(f):
            print(f"[WARN] Removing invalid or empty file: {f}")
            os.remove(f)
            removed = True
    return removed

def sra_route(accession, fastq_dir, metadata_dir, combined_meta, debug_lock_path, checkpoint_lock, threads, mem):
    max_retries = 3
    sra_file = os.path.join(fastq_dir, f"{accession}.sra")
    
    for i in range(1, max_retries + 1):
        print(f"Prefetching SRA file for {accession}: attempt {i}")
        cmd_prefetch = ["prefetch", str(accession), "--max-size", "200G", "--output-file", str(sra_file)]
        try:
            run_command(cmd_prefetch, f"ERROR: {accession} prefetch failed", debug_lock_path, accession)
            break  # success: exit the loop
        except RuntimeError as e:
            log_debug_message(debug_lock_path, f"Attempt {i} failed: {e}")
            if i == max_retries:
                return False  # only return False after final failed attempt
            time.sleep(2)


    print(f"Running vdb-validate for {accession}")
    cmd_validate = ["vdb-validate", sra_file]
    try:
        result = subprocess.run(cmd_validate, capture_output=True, text=True, check=True)
        print(" vdb-validate passed")
    except subprocess.CalledProcessError as e:
        log_debug_message(debug_lock_path, f"ERROR: {accession} vdb-validate failed:\n{e.stderr}")
        return False
    except Exception as e:
        log_debug_message(debug_lock_path, f"ERROR: {accession} vdb-validate exception: {e}")
        return False

    print(f"Converting SRA to FASTQ for {accession}")
    cmd_fasterq = [
        "fasterq-dump", sra_file, "--outdir", str(fastq_dir),
        "--threads", str(threads), "--mem", str(mem), "--split-files", "--include-technical"
    ]
    try:
        run_command(cmd_fasterq, f"ERROR: {accession} fasterq-dump failed", debug_lock_path, accession)
    except RuntimeError as e:
        log_debug_message(debug_lock_path, str(e))
        return False
    time.sleep(2)

    success_meta = False
    csv_path = os.path.join(metadata_dir, f"{accession}-run-info.csv")
    for i in range(1, max_retries + 1):
        try:
            cmd_esearch = f'esearch -db sra -query "{accession}" | efetch -format runinfo > {csv_path}'
            ret = subprocess.call(cmd_esearch, shell=True)
            if ret == 0 and os.path.isfile(csv_path) and os.path.getsize(csv_path) > 0:
                success_meta = True
                tsv_file = os.path.join(metadata_dir, f"{accession}-run-info.tsv")
                try:
                    with open(csv_path, "r") as inf, open(tsv_file, "w") as outf:
                        for line in inf:
                            outf.write(line.replace(",", "\t"))
                    with open(checkpoint_lock, "a+") as lock_fd:
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
                    log_debug_message(debug_lock_path, f"WARNING: {accession} SRA metadata handling error: {e}")
                break
            else:
                raise RuntimeError(f"Metadata fetch failed with return code {ret}")
        except Exception as e:
            log_debug_message(debug_lock_path, f"WARNING: {accession} SRA metadata attempt {i} failed: {e}")
            if i < max_retries:
                time.sleep(2 ** i)  # Exponential backoff
            else:
                log_debug_message(debug_lock_path, f"ERROR: {accession} All SRA metadata attempts failed")
                success_meta = False

    if not success_meta:
        return False

    all_fastqs = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq"))
    if not all_fastqs:
        log_debug_message(debug_lock_path, f"ERROR: {accession} No .fastq from SRA route.")
        return False

    for f in all_fastqs:
        if not is_valid_fastq(f, min_size_bytes=1024):
            log_debug_message(debug_lock_path, f"WARNING: {accession} Invalid FASTQ file: {f}")
            return False

    return True

def classify_fastq_by_read_type(accession, fastq_dir, debug_lock_path):
    try:
        print(f"[Classify] Getting READ_TYPE for {accession}")
        cmd_read_type = f"vdb-dump {accession} -C READ_TYPE 2>/dev/null | head -n 1"
        result = subprocess.run(cmd_read_type, shell=True, capture_output=True, text=True, check=True)
        line = result.stdout.strip()

        if "READ_TYPE:" in line:
            types_part = line.split(":", 1)[1].strip()
            read_types = [t.strip() for t in types_part.split(",")]

            bio_dir = os.path.join(fastq_dir, "fastq_biologicaldata")
            tech_dir = os.path.join(fastq_dir, "fastq_technicaldata")
            os.makedirs(bio_dir, exist_ok=True)
            os.makedirs(tech_dir, exist_ok=True)

            for i, read_type in enumerate(read_types):
                # First try suffixed filename (e.g., accession_1.fastq.gz)
                gz_suffixed = f"{accession}_{i+1}.fastq.gz"
                src_suffixed = os.path.join(fastq_dir, gz_suffixed)
            
                # Initialize source path as None
                src_path = None
                gz_filename = None
            
                if os.path.exists(src_suffixed):
                    src_path = src_suffixed
                    gz_filename = gz_suffixed
                elif i == 0:
                    # Only try unsuffixed for first read type
                    gz_unsuffixed = f"{accession}.fastq.gz"
                    src_unsuffixed = os.path.join(fastq_dir, gz_unsuffixed)
                    if os.path.exists(src_unsuffixed):
                        src_path = src_unsuffixed
                        gz_filename = gz_unsuffixed
            
                # If neither file was found, raise an error
                if src_path is None:
                    raise FileNotFoundError(
                        f"Missing FASTQ file for accession {accession}, read_type {read_type}. "
                        f"Tried: {src_suffixed}" + (f" and {src_unsuffixed}" if i == 0 else "")
                    )
            
                # Determine destination directory (biological or technical)
                dest_dir = bio_dir if "BIOLOGICAL" in read_type else tech_dir
                dest_path = os.path.join(dest_dir, gz_filename)
                print(f"  Moving {src_path} → {dest_path}")
                shutil.move(src_path, dest_path)

                
        else:
            log_debug_message(debug_lock_path, f"WARNING: {accession} Could not parse READ_TYPE info")
    except subprocess.CalledProcessError as e:
        log_debug_message(debug_lock_path, f"WARNING: {accession} vdb-dump failed: {e.stderr}")
    except Exception as e:
        log_debug_message(debug_lock_path, f"WARNING: {accession} FASTQ classification failed: {e}")

##########################
# Main entry point
##########################
if __name__ == "__main__":
    check_required_tools([
        "prefetch", "fasterq-dump", "esearch", "efetch", "wget", "curl", "gzip", "pigz", "vdb-dump"
    ])

    args = parse_cli_args()

    if args.accession:
        os.environ["ACCESSION"] = args.accession
    if args.workdir:
        os.environ["WORKDIR"] = args.workdir
    if args.metadata:
        os.environ["COMBINED_METADATA"] = args.metadata

    accession       = read_env_var("ACCESSION")
    workdir         = read_env_var("WORKDIR")
    checkpoint_file = read_env_var("CHECKPOINT_FILE")
    checkpoint_lock = read_env_var("CHECKPOINT_LOCK_FILE")
    combined_meta   = read_env_var("COMBINED_METADATA")
    debug_lock      = read_env_var("DEBUG_LOCK")
    token_file      = read_env_var("TOKEN_FILE")
    tokenlock_file  = read_env_var("TOKEN_LOCK_FILE")
    threads = os.cpu_count()
    mem_bytes = psutil.virtual_memory().available
    mem_gb = int(mem_bytes / (1024**3))
    mem = f"{mem_gb}G"

    fastq_dir    = os.path.join(workdir, "fastq_data")
    metadata_dir = os.path.join(workdir, "metadata")
    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(metadata_dir, exist_ok=True)

    # Ensure logs directory exists
    os.makedirs(os.path.join(workdir, "logs"), exist_ok=True)


    print(f"[DEBUG] Accession={accession}", flush=True)

    try:
        try:
            ok = sra_route(accession, fastq_dir, metadata_dir, combined_meta, debug_lock, checkpoint_lock, threads, mem)
            if not ok:
                msg = f"[sra_route] SRA route failed for {accession}"
                log_debug_message(debug_lock, msg)
                raise RuntimeError(msg)
        except Exception as e:
            log_debug_message(debug_lock, f"[sra_route] Exception: {e}")
            raise
    
        try:
            compress_fastqs(fastq_dir, accession, debug_lock)
        except Exception as e:
            log_debug_message(debug_lock, f"[compress_fastqs] Exception: {e}")
            raise
    
        try:
            gz_pattern = os.path.join(fastq_dir, f"{accession}*.fastq.gz")
            found_gz = glob.glob(gz_pattern)
            if not found_gz:
                msg = f"[glob] No FASTQ files found after compressing {accession}"
                log_debug_message(debug_lock, msg)
                raise RuntimeError(msg)
        except Exception as e:
            log_debug_message(debug_lock, f"[glob] Exception: {e}")
            raise
    
        try:
            for gz_file in found_gz:
                if not is_valid_gzip(gz_file):
                    msg = f"[is_valid_gzip] Invalid gzip file: {gz_file}"
                    log_debug_message(debug_lock, msg)
                    raise RuntimeError(msg)
        except Exception as e:
            log_debug_message(debug_lock, f"[is_valid_gzip] Exception: {e}")
            raise
    
        try:
            classify_fastq_by_read_type(accession, fastq_dir, debug_lock)
        except Exception as e:
            log_debug_message(debug_lock, f"[classify_fastq_by_read_type] Exception: {e}")
            raise
    
        try:
            cleanup_invalid_fastqs(fastq_dir, accession)
        except Exception as e:
            log_debug_message(debug_lock, f"[cleanup_invalid_fastqs] Exception: {e}")
            raise
    
        try:
            sra_file_path = os.path.join(fastq_dir, f"{accession}.sra")
            remove_file_safely(sra_file_path, accession, debug_lock)
        except Exception as e:
            log_debug_message(debug_lock, f"[remove_file_safely] Exception: {e}")
            raise
    
        try:
            append_to_checkpoint(accession, checkpoint_file, checkpoint_lock)
        except Exception as e:
            log_debug_message(debug_lock, f"[append_to_checkpoint] Exception: {e}")
            raise
    
        print(f"Success: {accession}")
    
    except Exception as e:
        log_debug_message(debug_lock, f"[FAILURE] Accession {accession} processing failed: {e}")
        sys.exit(1)
