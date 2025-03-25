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

##########################
# Helper Functions
##########################

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
    In practice, you might prefer a separate lock file. This is a demo.
    """
    debug_log_path = os.path.join(workdir, "logs", "debug.log")
    try:
        # We open the lock file (which might be the same file or a separate .lock)
        with open(debug_lock_path, "a+") as lk:
            # If you want to truly lock it, do flock_exclusive(lk).
            # flock_exclusive(lk)
            lk.write(f"{time.ctime()} {message}\n")
            # flock_release(lk)
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
    - Unlocks
    """
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
    """Find all uncompressed *.fastq and gzip them."""
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
    This is a heuristic. Adjust for your naming scheme if needed.
    e.g. ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/...
    """
    base_ftp = "ftp.sra.ebi.ac.uk/vol1/fastq"
    # We might guess subdirectories based on partial expansions of the accession
    # e.g. ERR1193450 -> "ERR119", "3450" => There's variation in how ENA splits directories.
    # We'll produce some plausible guesses.
    # For example, subdir approach:
    #   - ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/003/ERR1193450
    #   - ftp.sra.ebi.ac.uk/vol1/fastq/ERR119/450/ERR1193450
    # We'll try a few. This is a *demo*; real logic can be more elaborate.

    # Common pattern: The accession splitted after the first 6 or so digits, or last 3
    # e.g. ERR119/3450 => "ERR1193450" or "ERR119" + "3450"
    # Another approach is: (ERR119)(3450) => /ERR119/3450/ERR1193450. We'll try both.
    # We'll also try no sub-subdir, in case the run is stored in a simpler path.

    prefix = accession[:6]  # e.g. ERR119
    suffix = accession[6:]   # e.g. 3450

    # Potential subfolders to guess
    paths = [
        f"{base_ftp}/{prefix}/{suffix}/{accession}",
        f"{base_ftp}/{prefix}/00{suffix[-1]}/{accession}",   # 3450 -> last digit = 0 => /003
        f"{base_ftp}/{prefix}/{suffix[-3:]}/{accession}",    # last 3 digits
        f"{base_ftp}/{prefix}/{accession}",                  # no sub-subdir
    ]
    # We'll remove duplicates
    unique_paths = list(dict.fromkeys(paths))
    return unique_paths

def try_enaDataGet(accession, fastq_dir, debug_lock):
    """
    Attempt to download using enaDataGet.
    Return True if successful, False otherwise.
    """
    cmd = ["enaDataGet", "-f", "fastq", "-d", fastq_dir, accession]
    try:
        run_command(cmd, f"❌ ERROR: {accession} enaDataGet failed")
        return True
    except RuntimeError as e:
        log_debug_message(debug_lock, str(e))
        return False

def try_direct_wget(accession, fastq_dir, debug_lock, max_retries=2):
    """
    Attempt direct wget from guessed ENA ftp paths.
    We'll guess a few possible subdirectory patterns, each with .fastq.gz filenames.
    Return True if any file is successfully downloaded, else False.
    """
    success = False
    possible_dirs = guess_ena_ftp_paths(accession)
    # We assume the actual files might be like: {accession}_1.fastq.gz, {accession}_2.fastq.gz, or {accession}.fastq.gz
    candidates = [
        f"{accession}.fastq.gz",
        f"{accession}_1.fastq.gz",
        f"{accession}_2.fastq.gz",
    ]
    for dir_url in possible_dirs:
        for f in candidates:
            url = f"https://{dir_url}/{f}"
            # We'll do a few retries
            for attempt in range(1, max_retries+1):
                cmd = ["wget", "-O", os.path.join(fastq_dir, f), url]
                try:
                    run_command(cmd, f"Wget attempt {attempt} for {accession} failed: {url}")
                    success = True
                    break  # break out of the attempt loop
                except RuntimeError as e:
                    log_debug_message(debug_lock, str(e))
                    time.sleep(2 * attempt)
            if success:
                # If we got one file, let's continue with the next candidate
                # But maybe we want to see if there's also a _2?
                # We'll not break from the entire loop, so we continue to get both _1 and _2 if they exist.
                pass
        # After finishing candidates in one dir
    # We'll consider success if we at least got one .fastq.gz
    downloaded = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz"))
    return len(downloaded) > 0

def try_ftp_fallback(accession, fastq_dir, debug_lock):
    """
    Another fallback with 'curl' or 'wget' from an ftp: path
    Usually the URL is ftp://ftp.sra.ebi.ac.uk/vol1/fastq/<subdir>/...
    Return True if any file is downloaded, else False.
    """
    # We'll just re-use the guess_ena_ftp_paths, but with ftp://
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
                log_debug_message(debug_lock, str(e))
                # last resort: try wget ftp
                cmd2 = ["wget", "-O", out_path, url]
                try:
                    run_command(cmd2, f"wget ftp fallback for {accession} failed: {url}")
                    success = True
                except RuntimeError as e2:
                    log_debug_message(debug_lock, str(e2))
            # If success, keep going to try to get _2 if it exists
    downloaded = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq.gz"))
    return len(downloaded) > 0

def fetch_ena_metadata_and_merge(accession, fastq_dir, metadata_dir, combined_meta, checkpoint_lock, debug_lock):
    """
    Since we're not using fastq-dl's built-in metadata, we might rely on
    'enaDataGet' generating some metadata or do custom steps. Here is a placeholder
    that you can adapt to fetch ENA metadata. For now, let's just do a no-op or
    a hypothetical 'enaBrowserTools' run info extraction.
    """
    # If enaDataGet succeeded, it usually places a directory with run data and maybe
    # a readme. The official 'enaDataGet' doesn't automatically produce a single TSV
    # with metadata. You can parse the 'run.xml' or other info if you want to
    # unify with your combined_meta. This is just a placeholder.

    # lock combined meta if needed
    pass

def sra_route(accession, fastq_dir, metadata_dir, combined_meta, debug_lock, checkpoint_lock):
    """
    Attempt SRA route: prefetch + fasterq-dump + SRA metadata.
    Return True if successful, else False.
    """
    sra_files = os.path.join(fastq_dir, f"{accession}.sra")

    # Prefetch
    print(f"🔹 [Fallback] Prefetching SRA file for {accession}")
    cmd_prefetch = ["prefetch", accession, "--output-file", sra_files]
    try:
        run_command(cmd_prefetch, f"❌ ERROR: {accession} prefetch failed")
    except RuntimeError as e:
        log_debug_message(debug_lock, str(e))
        return False
    time.sleep(2)

    # fasterq-dump
    print(f"🔹 [Fallback] Converting SRA to FASTQ for {accession}")
    cmd_fasterq = [
        "fasterq-dump",
        sra_files,
        "--outdir", fastq_dir,
        "--threads", "4",
        "--mem", "8G",
        "--split-3"
    ]
    try:
        run_command(cmd_fasterq, f"❌ ERROR: {accession} fasterq-dump failed")
    except RuntimeError as e:
        log_debug_message(debug_lock, str(e))
        return False
    time.sleep(2)

    # fetch SRA metadata
    max_retries  = 3
    success_meta = False
    csv_path = os.path.join(metadata_dir, f"{accession}-run-info.csv")

    for i in range(1, max_retries + 1):
        cmd_esearch = (
            f'esearch -db sra -query "{accession}" | efetch -format runinfo > {csv_path}'
        )
        ret = subprocess.call(cmd_esearch, shell=True)
        if ret == 0 and os.path.isfile(csv_path) and os.path.getsize(csv_path) > 0:
            success_meta = True
            # Convert CSV -> TSV
            tsv_file = os.path.join(metadata_dir, f"{accession}-run-info.tsv")
            try:
                with open(csv_path, "r") as inf, open(tsv_file, "w") as outf:
                    for line in inf:
                        outf.write(line.replace(",", "\t"))
                # Merge into combined
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
                # Cleanup
                remove_file_safely(tsv_file, accession, debug_lock)
                remove_file_safely(csv_path, accession, debug_lock)
            except Exception as e:
                log_debug_message(debug_lock, f"⚠️ WARNING: {accession} SRA metadata handling error: {e}")
            break
        else:
            log_debug_message(debug_lock, f"⚠️ WARNING: {accession} SRA metadata attempt {i} failed")
            time.sleep(i)

    # If we got here, success_meta might be True or False
    if not success_meta:
        log_debug_message(debug_lock, f"❌ ERROR: {accession} All SRA metadata attempts failed")

    # Check if we got any .fastq files
    all_fastqs = glob.glob(os.path.join(fastq_dir, f"{accession}*.fastq"))
    if not all_fastqs:
        log_debug_message(debug_lock, f"❌ ERROR: {accession} No .fastq from fallback SRA route.")
        return False

    return True

##########################
# Main entry point
##########################
if __name__ == "__main__":
    # 1) Parse env vars
    accession       = read_env_var("ACCESSION")
    workdir         = read_env_var("WORKDIR")
    checkpoint_file = read_env_var("CHECKPOINT_FILE")
    checkpoint_lock = read_env_var("CHECKPOINT_LOCK_FILE")
    combined_meta   = read_env_var("COMBINED_METADATA")
    debug_lock      = read_env_var("DEBUG_LOCK")
    token_file      = read_env_var("TOKEN_FILE")
    tokenlock_file  = read_env_var("TOKEN_LOCK_FILE")

    fastq_dir       = os.path.join(workdir, "fastq_data")
    metadata_dir    = os.path.join(workdir, "metadata")
    os.makedirs(fastq_dir, exist_ok=True)
    os.makedirs(metadata_dir, exist_ok=True)

    @atexit.register
    def on_exit():
        release_token()

    global accession
    print(f"[DEBUG] Accession={accession}", flush=True)

    # 2) Determine if it's SRA or ENA by prefix
    prefix = accession[:3].upper()
    if prefix in ("SRR","SRX","SRS","SRP"):
        provider = "sra"
    elif prefix in ("ERR","ERX","ERS","ERP","DRR","DRX","DRS","DRP"):
        provider = "ena"
    else:
        log_debug_message(debug_lock, f"❌ ERROR: {accession} Unrecognized prefix => no known route.")
        sys.exit(1)

    # 3) If SRA => do the SRA route
    if provider == "sra":
        print(f"[INFO] Using SRA route for {accession}")
        ok = sra_route(accession, fastq_dir, metadata_dir, combined_meta, debug_lock, checkpoint_lock)
        if not ok:
            sys.exit(1)

    # 4) If ENA => attempt multiple methods:
    else:
        print(f"[INFO] Using ENA route(s) for {accession}")
        # 4a) Try enaDataGet
        ena_ok = try_enaDataGet(accession, fastq_dir, debug_lock)
        if not ena_ok:
            # 4b) If that fails => direct wget
            print(f"[WARN] enaDataGet failed => trying direct wget for {accession}")
            ena_ok = try_direct_wget(accession, fastq_dir, debug_lock, max_retries=2)
            if not ena_ok:
                # 4c) If that fails => ftp fallback with curl/wget
                print(f"[WARN] direct wget failed => trying ftp fallback for {accession}")
                ena_ok = try_ftp_fallback(accession, fastq_dir, debug_lock)
                if not ena_ok:
                    # 4d) If that fails => fallback to SRA route
                    print(f"[WARN] All ENA methods failed => fallback to SRA route for {accession}")
                    sra_ok = sra_route(accession, fastq_dir, metadata_dir, combined_meta, debug_lock, checkpoint_lock)
                    if not sra_ok:
                        log_debug_message(debug_lock, f"❌ ERROR: {accession} ENA + SRA fallback all failed.")
                        sys.exit(1)
                    else:
                        # If SRA fallback worked, we continue below
                        pass
                else:
                    print(f"[INFO] ftp fallback succeeded for {accession}")
            else:
                print(f"[INFO] direct wget succeeded for {accession}")
        else:
            print(f"[INFO] enaDataGet succeeded for {accession}")

        # If the final ENA route we used was enaDataGet, directWget, or ftp fallback
        # and it succeeded, we can do metadata merges if you want:
        # fetch_ena_metadata_and_merge(accession, fastq_dir, metadata_dir, combined_meta, checkpoint_lock, debug_lock)

    # 5) Compress any leftover .fastq
    compress_fastqs(fastq_dir, accession, debug_lock)

    # 6) Check we have .fastq.gz
    gz_pattern = os.path.join(fastq_dir, f"{accession}*.fastq.gz")
    found_gz = glob.glob(gz_pattern)
    if not found_gz:
        log_debug_message(debug_lock, f"❌ ERROR: {accession} No FASTQ files found after process.")
        sys.exit(1)

    # 7) Append to checkpoint
    append_to_checkpoint(accession, checkpoint_file, checkpoint_lock)

    # 8) Cleanup any leftover .sra if it existed
    sra_file_path = os.path.join(fastq_dir, f"{accession}.sra")
    remove_file_safely(sra_file_path, accession, debug_lock)

    print(f"✅ Successfully processed {accession}")
