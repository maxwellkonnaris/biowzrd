#!/usr/bin/env python3

"""
counts_16s.py
=============

A SLURM‑aware pipeline for validating, repairing and processing 16S‑rRNA
FASTQ files prior to DADA2 denoising.  The script can either

1.  Perform **one‑shot** processing on a full accession list, or
2.  Act as a **job‑submitter**, spawning a separate SLURM job per
    BioProject (via the ``--submit-only`` flag).

Key features
------------
* **Checksum‑driven integrity checks** with BLAKE3 (``b3sum``).
* Parallel FASTQ validation / repair using ``concurrent.futures``.
* Automatic CSV checkpoint updates guarded with file locks.
* Resource‑aware SLURM submission that scales memory / CPU per project.
* Environment autodetection for Micromamba‑managed R/DADA2 envs.

Typical usage
-------------
.. code-block:: bash

    python counts_16s.py \
        -i metadata.csv \
        -d fastq_data/fastq_biologicaldata \
        --submit-only          # submit one SLURM job per BioProject

Or run the whole pipeline locally:

.. code-block:: bash

    python counts_16s.py -i metadata.csv -w 12
"""

import argparse
import csv
from datetime import datetime
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import time
import fcntl
import tempfile
import shlex
from typing import List, Dict, Optional, Tuple

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

#try:
#    import fireducks.pandas as pd
#    logger.info("Using fireducks.pandas for faster multithreaded I/O")
#except ImportError:
import pandas as pd
logger.info("Falling back to pandas for single threaded I/O")
    
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# SLURM-like configuration 
SLURM_CONFIG = {
    "job-name": "counts_16s",
    "output": "slurm_16s-%j.out",
    "error": "slurm_16s-%j.err",
    "time": "48:00:00",
    "nodes": 1,
    "ntasks": 1,
    "cpus-per-task": 16,
    "mem": "64G",
    "account": "open"
}

def submit_slurm_job(dependency=None):
    """
    Submit the pipeline as a SLURM job if not already running under SLURM.
    """
    if "SLURM_JOB_ID" in os.environ:
        logger.info("Already running under SLURM, skipping submission")
        return

    # Build command without --submit-slurm to avoid recursion
    cmd = [sys.executable] + [arg for arg in sys.argv if arg != "--submit-slurm"]

    # Write temporary SLURM script
    with tempfile.NamedTemporaryFile("w", suffix=".slurm", delete=False) as sb:
        sb.write("#!/bin/bash\n")
        for key, val in SLURM_CONFIG.items():
            sb.write(f"#SBATCH --{key}={val}\n")
        if dependency:
            sb.write(f"#SBATCH --dependency=afterok:{dependency}\n")
        sb.write("\n")
        sb.write("source ~/.bashrc\n")
        active_env = os.environ.get("CONDA_DEFAULT_ENV") or os.environ.get("MAMBA_DEFAULT_ENV") or "dada2"  # Fallback to 'dada2'
        if not active_env:
            logger.error("Could not determine active micromamba environment")
            sys.exit(1)
        sb.write(f"micromamba activate {active_env}\n\n")
        sb.write(" ".join(cmd) + "\n")
        script_path = sb.name

    # Submit the job
    try:
        res = subprocess.run(
            ["sbatch", script_path],
            capture_output=True, text=True, check=True
        )
        logger.info(f"Submitted SLURM job: {res.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        logger.error(f"SLURM submission failed: {e.stderr.strip()}")
        logger.error(f"Command output: {e.stdout.strip()}")
        sys.exit(1)
    finally:
        os.unlink(script_path)

# Default variables and directories
DEFAULT_DIR = "fastq_data/fastq_biologicaldata"
OUTPUT_BASE = Path("dada2_16s")
LOCK_DIR = Path("locks")
FAILED_FILE = Path("failed_16s.log")
RDP_DATABASE = "rdp_19_toGenus_trainset.fa.gz"

# these must already be in df, or you’ll raise an error
HARD_REQUIRED = ["Bioproject", "RunAccession", "SequencingType"]

# these will be added if missing, with “Validated”→"0", fastq cols→""
OPTIONAL_COLUMNS = ["Fastq1", "Fastq2", "Validated"]
ZERO_DEFAULTS = {"Validated"}

def ensure_required_columns(df: pd.DataFrame,
                            required: list[str],
                            zero_cols: set[str]) -> pd.DataFrame:
    """
    Guarantee that `df` has every column in `required`.
    For any missing col in zero_cols, fill with "0"; otherwise fill with "".
    """
    for col in required:
        if col not in df.columns:
            df[col] = "0" if col in zero_cols else ""
    return df

def prepare_16s_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    1) Error if any HARD_REQUIRED columns are missing.
    2) Auto-add any OPTIONAL_COLUMNS that aren’t present,
       defaulting “Validated” to "0", fastq paths to "".
    """
    missing = set(HARD_REQUIRED) - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {missing!r}")
    return ensure_required_columns(df, OPTIONAL_COLUMNS, ZERO_DEFAULTS)

def setup_directories():
    """Create required directories."""
    for d in [LOCK_DIR, OUTPUT_BASE]:
        d.mkdir(parents=True, exist_ok=True)

def append_with_flock(line: str, file_path: Path):
    """
    Append a line to a file using an exclusive file lock
    so multiple processes can’t write simultaneously.
    """
    with file_path.open("a") as f:
        # Acquire an exclusive lock
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(line + "\n")
        # Release the lock
        fcntl.flock(f, fcntl.LOCK_UN)

def parse_arguments():
    """Parse command-line arguments."""
    cpu_count = os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 8)
    if isinstance(cpu_count, str):
        cpu_count = int(cpu_count)
    if cpu_count is None:
        logger.warning("CPU count unavailable; falling back to 8 CPUs")
        cpu_count = 8
    default_workers = cpu_count * 3 // 4  # Default to 3/4ths of CPUs

    parser = argparse.ArgumentParser(description="16S-rRNA FASTQ processing pipeline")
    parser.add_argument("-i", "--input-file", required=True, help="Input CSV file")
    parser.add_argument("-d", "--fastq-dir", default=DEFAULT_DIR, help="FASTQ directory")
    parser.add_argument("-w", "--num-workers", type=int, default=default_workers,
                       help=f"Number of workers (default: 3/4ths of CPUs = {default_workers})")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("--no-repair-fastqs", action="store_false", dest="repair_fastqs",
                       help="Disable FASTQ repair")
    parser.add_argument("--submit-only", action="store_true",
                        help="Skip all processing and just submit per-Bioproject jobs")
    parser.add_argument("--submit-slurm", action="store_true",
                        help="Submit this script as a SLURM job and exit")
    parser.add_argument("--dependency", type=str, default=None,
                        help="SLURM job ID to wait for (afterok)")
    return parser.parse_args()

def setup_environment():
    try:
        result = subprocess.run(["micromamba", "env", "list"], capture_output=True, text=True, check=True)
        envs = result.stdout.splitlines()
        for line in envs:
            if "*" in line:
                active_env = line.split()[0]
                break
        else:
            logger.error("No active micromamba environment detected")
            sys.exit(1)
        logger.info(f"Using activated DADA2 environment: {active_env}")
        return active_env
    except subprocess.CalledProcessError:
        logger.error("micromamba not in PATH")
        sys.exit(1)

def run_command(cmd: str, description: str, env_name: str = None) -> None:
    """
    Run a shell command (e.g. Rscript) directly, assuming the environment is already active.
    """
    logger.info(description)
    try:
        subprocess.run(cmd, shell=True, check=True, text=True)
    except subprocess.CalledProcessError:
        logger.error(f"Command failed: {description}")
        raise

def validate_input_file(file_path: Path) -> str:
    """
    Validate & normalize CSV header; return delimiter.

    1) Detects comma vs tab delimiter.
    2) Reads header, normalizes each column name:
       - all → lowercase, then capitalize first character
       - columns starting with "run" → "Run" + capitalize(rest)
       - columns starting with "sequencing" → "Sequencing" + capitalize(rest)
    3) If header changed, rewrite the file with the normalized header.
    4) Returns delimiter (',' or '\\t').

    Args
    ----
    file_path : pathlib.Path
        Path to the metadata CSV/TSV to validate.

    Returns
    -------
    str
        Delimiter string: either ',' or '\\t'.

    Raises
    ------
    SystemExit
        If no comma or tab can be found in the header.
    """
    text = file_path.read_text().splitlines()
    if not text:
        logger.error("Input file is empty")
        sys.exit(1)

    orig_header, *rest = text

    # 1) detect delimiter
    if ',' in orig_header:
        delim = ','
    elif '\t' in orig_header:
        delim = '\t'
    else:
        logger.error("Could not detect delimiter; header has neither comma nor tab")
        sys.exit(1)

    # 2) split into column names & normalize
    cols = [c.strip() for c in orig_header.split(delim)]

    def normalize(col: str) -> str:
        """lower → capitalize first char; special-cases run*/sequencing* prefixes."""
        cl = col.lower().replace(' ', '')
        if cl.startswith('run'):
            # e.g. "runaccession" → "Run" + "accession".capitalize()
            return 'Run' + cl[3:].capitalize()
        if cl.startswith('sequencing'):
            # e.g. "sequencingtype" → "Sequencing" + "type".capitalize()
            return 'Sequencing' + cl[10:].capitalize()
        # default: just Title-case
        return cl.capitalize()

    new_cols = [normalize(c) for c in cols]
    new_header = delim.join(new_cols)

    # 3) rewrite file if header changed
    if new_header != orig_header:
        logger.info(f"Normalizing header: {orig_header!r} → {new_header!r}")
        with file_path.open('w', newline='') as f:
            f.write(new_header + '\n')
            f.write('\n'.join(rest))

    # 4) return the delimiter so downstream code knows how to read the file
    return delim

def check_rdp_database():
    """Verify RDP database exists."""
    if not Path(RDP_DATABASE).is_file():
        logger.error(f"Missing RDP database: {RDP_DATABASE}")
        sys.exit(1)
    logger.info(f"RDP DB ok: {RDP_DATABASE}")

def build_validated_set(input_file: Path, delim: str) -> list[str]:
    """
    Read the input CSV and return all RunAccession values
    for which the 'Validated' field == '1'.
    """
    logger.info(f"Building validated set from {input_file}")
    try:
        df = pd.read_csv(
            input_file,
            sep=delim,
            usecols=["RunAccession", "Validated"],
            dtype={"RunAccession": str, "Validated": str}
        )
        validated = df.loc[df["Validated"] == "1", "RunAccession"].tolist()
        logger.info(f"Found {len(validated)} validated accessions")
        return validated
    except ValueError as e:
        logger.warning("No 'Validated' column yet — skipping.")
        return []

def process_fastq_global(fq_path_str: str, input_accs: set) -> Optional[Tuple[str, str]]:
    fq = Path(fq_path_str)
    acc = re.sub(r"(_[12])?\.fastq(\.gz)?$", "", fq.name)
    return (acc, fq_path_str) if acc in input_accs else None
    

def compute_checksum(f: Path) -> str:
    """Return the b3sum output for this file."""
    res = subprocess.run(
        ["b3sum", str(f)],
        capture_output=True, text=True, check=True
    )
    return res.stdout.strip()


def generate_fastq_checksums(fastq_dir: Path, num_workers: int):
    """Generate BLAKE3 checksums for FASTQ files, with a progress bar."""
    fastq_files = sorted(fastq_dir.glob("*.fastq.gz"))
    checksum_file = fastq_dir / "checksums.b3"
    checksums = []
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        for checksum in tqdm(
            executor.map(compute_checksum, fastq_files),
            total=len(fastq_files),
            desc="Checksumming FASTQs"
        ):
            checksums.append(checksum)
    with checksum_file.open("w") as f:
        f.write("\n".join(checksums) + "\n")
    logger.info(f"Checksums saved to {checksum_file}")

def validate_fastq(fq: Path) -> bool:
    """
    Validate a FASTQ file.
    Logs failures at INFO; successes at DEBUG.
    """
    logger.debug(f"[VAL] Checking {fq}")
    if not fq.is_file():
        logger.error(f"[VAL FAIL] {fq} does not exist")
        return False
    try:
        with gzip.open(fq, "rt") as f:
            lines = [f.readline().strip() for _ in range(4)]
            if len(lines) < 4 or not lines[0].startswith("@") or not lines[2].startswith("+"):
                logger.error(f"[VAL FAIL] {fq} → Invalid FASTQ structure")
                return False
    except (gzip.BadGzipFile, EOFError):
        logger.error(f"[VAL FAIL] {fq} is not valid gzip")
        return False
    logger.debug(f"[VAL PASS] {fq}")
    return True

def is_biological_fastq(fq: Path) -> bool:
    """Check if a FASTQ file contains biological sequences."""
    with gzip.open(fq, "rt") as f:
        lengths = []
        for i, line in enumerate(f):
            if i % 4 == 1:
                lengths.append(len(line.strip()))
            if len(lengths) >= 20:
                break
    uniq_lengths = len(set(lengths))
    max_length = max(lengths) if lengths else 0
    logger.debug(f"is_biological_fastq: {fq}, uniq={uniq_lengths}, maxlen={max_length}")
    return uniq_lengths > 1 or max_length > 29

def repair_fastq_if_needed(
    fastq_dir: Path,
    num_workers: int,
    validated_accessions: list[str]
) -> Path:
    """Repair corrupted FASTQ files and return new FASTQ dir."""
    repaired_dir = fastq_dir.parent / f"{fastq_dir.name}_repaired"
    checksum_file = fastq_dir / "checksums.b3"
    failed_log = fastq_dir / "failed_checksums.txt"

    repaired_dir.mkdir(exist_ok=True)
    failed_log.write_text("")

    # 1) Ensure we have an initial checksum list
    if not checksum_file.exists():
        generate_fastq_checksums(fastq_dir, num_workers)

    # 2) Build list of (hash, filename) to check, skipping already-validated
    skip = {f"{acc}.fastq.gz" for acc in validated_accessions}
    to_check: list[tuple[str, str]] = []
    with checksum_file.open() as f:
        for line in f:
            hash_val, fname = line.strip().split(maxsplit=1)
            if Path(fname).name not in skip:
                to_check.append((hash_val, fname))

    # 3) Verify checksums in parallel, collecting any failures
    def verify(h: str, fn: str) -> str | None:
        cmd = f"echo '{h}  {fn}' | b3sum -c - --quiet"
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            return None
        except subprocess.CalledProcessError:
            return fn

    failed_files: list[str] = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(verify, h, fn) for h, fn in to_check]
        for future in tqdm(as_completed(futures),
                           total=len(futures),
                           desc="Verifying FASTQs"):
            bad = future.result()
            if bad:
                failed_files.append(bad)

    # 4) Write out the failures all at once
    failed_log.write_text("\n".join(failed_files) + "\n")
    logger.info(f"{len(failed_files)} files failed checksum")

    # 5) Repair corrupted files if any
    if failed_files:
        logger.info("Repairing corrupted FASTQ files …")
        def do_repair(fn: str):
            inp = Path(fn)
            out = repaired_dir / inp.name
            with gzip.open(inp, "rt") as inf, gzip.open(out, "wt") as outf:
                for line in inf:
                    outf.write(line.rstrip("\r") + "\n")

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            list(tqdm(executor.map(do_repair, failed_files),
                      total=len(failed_files),
                      desc="Repairing FASTQs"))

    # 6) Symlink any untouched files
    for fq in fastq_dir.glob("*.fastq.gz"):
        dest = repaired_dir / fq.name
        if not dest.exists():
            dest.symlink_to(fq.resolve())

    # 7) Regenerate checksums on the repaired dir
    generate_fastq_checksums(repaired_dir, num_workers)
    logger.info(f"Repair complete → FASTQ_DIR updated to {repaired_dir}")

    return repaired_dir


def validate_and_check(fq: Path, acc: str) -> Tuple[str, Optional[Path], Optional[str]]:
    """
    Validate a FASTQ file and determine biological validity.
    Returns:
      - acc (str): the accession string
      - fq (Path or None): the Path if valid biological FASTQ, else None
      - reason (str or None): None if valid, else "TECHNICAL" or "CORRUPT"
    """
    if validate_fastq(fq):
        if is_biological_fastq(fq):
            return acc, fq, None
        else:
            return acc, None, "TECHNICAL"
    else:
        return acc, None, "CORRUPT"

def update_input_with_fastq_paths(
    input_file: Path,
    fastq_dir: Path,
    delim: str,
    num_workers: int
):
    logger.info("Updating FASTQ paths")

    # 1) Load input CSV
    df = pd.read_csv(input_file, sep=delim, dtype=str)
    df = prepare_16s_df(df)
    
    input_accs = set(df.loc[df["SequencingType"] == "16S", "RunAccession"].str.strip())
    logger.info(f"{len(input_accs)} target accessions from CSV")

    # 2) Scan FASTQ dir
    fastq_files = list(fastq_dir.glob("*.fastq.gz"))
    fq_paths = [str(p) for p in fastq_files]
    fastq_map: Dict[str, List[Path]] = {}

    # 2a) Map accession→file in parallel, with progress
    map_futures = []
    with ProcessPoolExecutor(max_workers=num_workers) as exe:
        for fq in fq_paths:
            map_futures.append(exe.submit(process_fastq_global, fq, input_accs))

        for future in tqdm(as_completed(map_futures),
                           total=len(map_futures),
                           desc="Scanning FASTQ files"):
            res = future.result()
            if res:
                acc, fq_path = res
                fastq_map.setdefault(acc, []).append(Path(fq_path))
    logger.info(f"Found FASTQs for {len(fastq_map)} of {len(input_accs)} accessions")

    # 3) Validate each FASTQ in parallel, with progress
    validated_map: Dict[str, List[Path]] = {}
    val_futures = []
    with ProcessPoolExecutor(max_workers=num_workers) as exe:
        for acc, fqs in fastq_map.items():
            for fq in fqs:
                val_futures.append(exe.submit(validate_and_check, fq, acc))

        for future in tqdm(as_completed(val_futures),
                           total=len(val_futures),
                           desc="Validating FASTQs"):
            acc, fq, reason = future.result()
            if reason:
                append_with_flock(f"{acc}:{reason}", FAILED_FILE)
            elif fq:
                validated_map.setdefault(acc, []).append(fq)

    # 4) Build a tiny DataFrame of Fastq1/Fastq2/Validated
    validated_data = []
    for acc, vfs in validated_map.items():
        if   not vfs:
            validated_data.append({"RunAccession": acc, "Fastq1":"MISSING","Fastq2":"MISSING","Validated":"0"})
        elif len(vfs)==1:
            validated_data.append({"RunAccession": acc, "Fastq1":str(vfs[0]),"Fastq2":"","Validated":"1"})
        elif len(vfs) == 2:
            fq_sorted = sorted(vfs, key=lambda p: ("_2" in p.name, p.name))
            validated_data.append({"RunAccession": acc,"Fastq1": str(fq_sorted[0]),"Fastq2": str(fq_sorted[1]),
                "Validated": "1"
            })
        else:
            fq_sorted = sorted(vfs, key=lambda p: ("_2" in p.name, p.name))
            validated_data.append({"RunAccession": acc,"Fastq1": str(fq_sorted[0]),"Fastq2": str(fq_sorted[1]),
                "Validated": "1"
            })
    validated_df = pd.DataFrame(validated_data)

    # 5) Merge back onto the original DataFrame
    #    Drop any old Fastq1/2/Validated columns if they exist,
    #    then left‑merge so that everything lines up.
    df = df.drop(columns=["Fastq1","Fastq2","Validated"], errors="ignore")
    if not validated_df.empty:
        df = df.merge(validated_df, on="RunAccession", how="left")
    else:
        df["Fastq1"] = df["Fastq2"] = df["Validated"] = "0"

    # Fill in any missing values
    df["Fastq1"]    = df["Fastq1"   ].fillna("MISSING")
    df["Fastq2"]    = df["Fastq2"   ].fillna("MISSING")
    df["Validated"] = df["Validated"].fillna("0")

    # And finally overwrite the CSV
    logger.info(f"Writing updated CSV → {input_file}")
    df[REQUIRED_COLUMNS_16S].to_csv(input_file, sep=delim, index=False)
    logger.info("FASTQ path update done")


def update_checkpoint(
    input_file: Path,
    accession: str,
    field: str,
    value: str,
    delim: str
):
    """
    Safely update a single checkpoint cell using a file lock
    so multiple processes can’t clobber each other.
    """
    lock_path = LOCK_DIR / (input_file.name + ".lock")
    # Ensure lock file exists
    lock_path.touch(exist_ok=True)

    # Acquire an exclusive OS‑level lock
    with lock_path.open("r+") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            df = pd.read_csv(input_file, sep=delim, dtype=str)

            if field not in df.columns:
                logger.warning(f"[CHECKPOINT] Field '{field}' not in columns — skipping")
                return

            match = df["RunAccession"] == accession
            if not match.any():
                logger.warning(f"[CHECKPOINT] Accession '{accession}' not found in CSV — skipping")
                return

            df.loc[match, field] = value

            # Backup and write
            bak = input_file.with_suffix(".bak")
            shutil.copy2(input_file, bak)
            df.to_csv(input_file, sep=delim, index=False)
            logger.debug(f"[CHECKPOINT] Updated {field}={value} for {accession}")
        except Exception as e:
            logger.error(f"[CHECKPOINT ERROR] Failed to update {field} for {accession}: {e}")
        finally:
            # Release the lock
            fcntl.flock(lk, fcntl.LOCK_UN)

def process_bioprojects(
    bioprojects: List[str],
    inputfile: str,
    threads: int
):
    """
    Submit one SLURM job per bioproject, passing:
      Rscript run_dada2.R <biop> <inputfile> <threads> <outdir>

    Each job will have:
      --job-name=<biop>
      --output=<biop>-%j.out
      --error=<biop>-%j.err

    Automatically estimates FASTQ size to suggest resource adjustments.
    """

    LOGS_DIR = OUTPUT_BASE / "logs"
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    # Read input CSV to get paths
    df = pd.read_csv(inputfile)

    for biop in bioprojects:
        outdir = OUTPUT_BASE / biop
        outdir.mkdir(parents=True, exist_ok=True)

        # Subset paths for the current bioproject
        subdf = df[df["Bioproject"] == biop]
        large_file_flag = False
        missing_paths = 0
        total_size = 0
        
        for path in subdf["Fastq1"]:
            try:
                size = Path(path).stat().st_size
                total_size += size
                if size > 10 * 1e9:  # Single file >10 GB
                    large_file_flag = True
            except FileNotFoundError:
                missing_paths += 1
        
        total_size_gb = total_size / 1e9
        num_files = len(subdf)
        
        # Decide resource allocation
        if large_file_flag:
            mem = "256G"
            cpus = 16
            time = "48:00:00"
            account = "one"
            print(f"Bioproject {biop} has at least one FASTQ > 10 GB. Increasing resources: {cpus} CPUs, {mem}, {time}.")
        elif total_size_gb > 150:
            mem = "256G"
            cpus = "32"
            time = "48:00:00"
            account = "one"
            print(f"Bioproject {biop} has {num_files} FASTQ files totaling {total_size_gb:.1f} GB. Increasing resources: {cpus} CPUs, {mem}, {time}.")
        else:
            mem = "128G"
            cpus = 8
            time = "48:00:00"
            account = "open"
        
        # write a temporary SLURM script
        with tempfile.NamedTemporaryFile("w", suffix=".slurm", delete=False) as sb:
            sb.write("#!/bin/bash\n")
            sb.write(f"#SBATCH --job-name={biop}\n")
            sb.write(f"#SBATCH --output={LOGS_DIR}/{biop}-%j.out\n")
            sb.write(f"#SBATCH --error={LOGS_DIR}/{biop}-%j.err\n")
            sb.write(f"#SBATCH --cpus-per-task={cpus}\n")
            sb.write(f"#SBATCH --mem={mem}\n")
            sb.write(f"#SBATCH --time={time}\n")
            sb.write(f"#SBATCH --account={account}\n")
            for key, val in SLURM_CONFIG.items():
                if key in ("job-name", "output", "error", "cpus-per-task", "mem", "time", "account"):
                    continue
                sb.write(f"#SBATCH --{key}={val}\n")
            sb.write("\n")
            sb.write("source ~/.bashrc\n")
            sb.write(
                f"Rscript {Path.cwd()}/run_dada2.R "
                f"{shlex.quote(biop)} {inputfile} {threads} {outdir}\n"
            )

        subprocess.run(["sbatch", sb.name], check=True)

    
def main():
    args = parse_arguments()
    from pathlib import Path

    # Determine delimiter early (needed for submit-only)
    delim = validate_input_file(Path(args.input_file))
    os.environ["DELIM"] = delim

    # ── SUBMIT-ONLY MODE ───────────────────────────────────────
    if args.submit_only:
        # 1) read your CSV (must include a 'Bioproject' column)
        df = pd.read_csv(args.input_file, sep=delim, dtype=str)
        # 2) extract unique BioProjects
        bioprojects = df["Bioproject"].unique().tolist()
        # 3) ensure output/log directories exist
        setup_directories()
        # 4) submit one Slurm job per BioProject
        process_bioprojects(bioprojects, args.input_file, args.num_workers)
        return
    # ────────────────────────────────────────────────────────────

    # Existing --submit-slurm behavior
    if args.submit_slurm:
        submit_slurm_job(dependency=args.dependency)
        return

    # Normal pipeline execution
    if args.debug:
        logger.setLevel(logging.DEBUG)

    # Cap workers to 3/4 of available CPUs
    cpu_count = os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 8)
    if isinstance(cpu_count, str):
        cpu_count = int(cpu_count)
    args.num_workers = max(1, min(args.num_workers, cpu_count * 3 // 4))
    logger.info(f"Using {args.num_workers} workers for non-DADA2 steps")

    # Setup
    setup_directories()
    env_name = setup_environment()
    os.environ["INPUT_FILE"] = args.input_file
    os.environ["FAILED_FILE"] = str(FAILED_FILE)
    os.environ["OUTPUT_BASE"] = str(OUTPUT_BASE)

    check_rdp_database()
    fastq_dir = Path(args.fastq_dir).resolve()

    # Step 1: Generate checksums if needed
    if not (fastq_dir / "checksums.b3").exists():
        generate_fastq_checksums(fastq_dir, args.num_workers)

    # Step 2: Build validated set
    validated_accessions = build_validated_set(Path(args.input_file), delim)

    # Step 3: Repair FASTQs if needed
    if args.repair_fastqs:
        fastq_dir = repair_fastq_if_needed(fastq_dir, args.num_workers, validated_accessions)

    # Step 4: Update CSV with FASTQ paths + validation
    update_input_with_fastq_paths(Path(args.input_file), fastq_dir, delim, args.num_workers)

    # Step 5: Submit per-BioProject jobs
    full_df = pd.read_csv(args.input_file, sep=delim, dtype=str)
    full_df = full_df[full_df['Validated'] == '1']
    bioprojects = full_df["Bioproject"].unique().tolist()
    inputfile = OUTPUT_BASE / "16s_validated.csv"
    full_df.to_csv(inputfile, sep=delim, index=False)
    process_bioprojects(bioprojects, str(inputfile), args.num_workers)

if __name__ == "__main__":
    main()

