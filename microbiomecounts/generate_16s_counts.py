#!/usr/bin/env python3

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
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# SLURM-like configuration (adapt as needed for your environment)
SLURM_CONFIG = {
    "job-name": "counts_16s",
    "output": "slurm_16s-%j.out",
    "error": "slurm_16s-%j.err",
    "time": "48:00:00",
    "nodes": 1,
    "ntasks": 1,
    "cpus-per-task": 32,
    "mem": "256G",
    "account": "one",
    "mail-user": "mak6930@psu.edu",
}

def submit_slurm_job():
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

# at the very top, after imports
REQUIRED_COLUMNS_16S = [
    "Bioproject", "RunAccession", "SequencingType",
    "Fastq1", "Fastq2", "Validated", "Dada2", "Completed"
]
ZERO_DEFAULTS_16S = {"Validated", "Dada2", "Completed"}

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
    parser.add_argument("--submit-slurm", action="store_true",
                        help="Submit this script as a SLURM job and exit")
    parser.add_argument("--process-sample", action="store_true", help="Process a single sample")
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
    """Validate CSV header and determine delimiter."""
    with file_path.open("r") as f:
        first_line = f.readline().strip()
    if first_line.startswith("Bioproject,RunAccession,SequencingType"):
        return ","
    if first_line.startswith("Bioproject\tRunAccession\tSequencingType"):
        return "\t"
    logger.error("Malformed header in input file")
    sys.exit(1)

def check_rdp_database():
    """Verify RDP database exists."""
    if not Path(RDP_DATABASE).is_file():
        logger.error(f"Missing RDP database: {RDP_DATABASE}")
        sys.exit(1)
    logger.info(f"RDP DB ok: {RDP_DATABASE}")

def build_validated_set(input_file: Path, delim: str) -> List[str]:
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


def generate_fastq_checksums(fastq_dir: Path, num_workers: int):
    """Generate BLAKE3 checksums for FASTQ files, with a progress bar."""
    fastq_files = sorted(fastq_dir.glob("*.fastq.gz"))
    checksum_file = fastq_dir / "checksums.b3"
    checksums = []
    def compute_checksum(f: Path) -> str:
        result = subprocess.run(["b3sum", str(f)], capture_output=True, text=True, check=True)
        return result.stdout.strip()
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
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


def _check_asv(path: Path) -> bool:
    """
    Return True if R can read `path` to a non‑empty table.
    """
    rscript = f"""
    args <- commandArgs(trailingOnly=TRUE)
    x <- readRDS(args[1])
    if ((is.matrix(x) || is.data.frame(x)) && ncol(x) > 0) {{
      quit(status=0)
    }} else {{
      quit(status=1)
    }}
    """
    try:
        subprocess.run(
            ["Rscript", "-e", rscript, str(path)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
        )
        return True
    except subprocess.CalledProcessError:
        return False

def regenerate_input_with_validation(
    input_file: Path,
    output_dir: Path,
    delim: str = ","
) -> pd.DataFrame:
    df = pd.read_csv(input_file, sep=delim, dtype=str)
    df["Dada2"]     = "0"
    df["Completed"] = "0"

    # 1) Mark merges done immediately
    for idx, row in df.iterrows():
        biop = row["Bioproject"]
        merged_file = output_dir/biop/f"{biop}_merged.rds"
        if merged_file.exists():
            df.at[idx,"Dada2"] = df.at[idx,"Completed"] = "1"

    # 2) Build list of per‑sample ASV files to verify
    tasks = []
    for idx, row in df.iterrows():
        if df.at[idx,"Dada2"] == "0":  # skip those already merged
            biop   = row["Bioproject"]
            acc    = row["RunAccession"]
            asv_fn = output_dir/biop/f"asv_{acc}.rds"
            if asv_fn.exists() and asv_fn.stat().st_size >= 100:
                tasks.append((idx, asv_fn))

    # 3) Verify in parallel with a tqdm progress bar
    with ProcessPoolExecutor(max_workers=os.cpu_count()//2 or 4) as exe:
        futures = {exe.submit(_check_asv, path): idx for idx, path in tasks}
        for future in tqdm(as_completed(futures),
                           total=len(futures),
                           desc="Verifying ASV files"):
            idx = futures[future]
            if future.result():
                df.at[idx, "Dada2"] = "1"

    # 4) Any of those with Dada2==1 but not merged are “in progress” → leave Completed=0

    # backup & write out…
    timestamp   = datetime.now().strftime("%Y%m%d-%H%M%S")
    bak_path    = input_file.with_suffix(f".bak.{timestamp}.csv")
    df.to_csv(bak_path, sep=delim, index=False)
    logger.info(f"[BACKUP] {bak_path}")
    df.to_csv(input_file, sep=delim, index=False)
    logger.info(f"[REGEN] Updated {input_file}")
    return df

def update_input_with_fastq_paths(
    input_file: Path,
    fastq_dir: Path,
    delim: str,
    num_workers: int
):
    logger.info("Updating FASTQ paths")

    # 1) Load input CSV
    df = pd.read_csv(input_file, sep=delim)
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

    df = ensure_required_columns(df, REQUIRED_COLUMNS_16S,ZERO_DEFAULTS_16S)

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


def process_sample(biop, acc, st, fq1, fq2, val, d2, comp, env_name):
    if st != "16S":
        return
    if val == "0":
        return
    if not (fq1 or fq2):
        append_with_flock(f"{acc}:NO_FASTQ", FAILED_FILE)
        return
    if comp == "1":
        return

    outdir = OUTPUT_BASE / biop
    outdir.mkdir(parents=True, exist_ok=True)
    seqtab = outdir / f"asv_{acc}.rds"

    # Only do DADA2 *and* mark Completed if d2 != "1"
    if d2 != "1":
        cmd = (
            f"Rscript {Path.cwd()}/run_dada2_partial.R {fq1} {fq2} {seqtab}"
            if fq2 else
            f"Rscript {Path.cwd()}/run_dada2_partial.R {fq1} {seqtab}"
        )
        env = os.environ.copy()
        env.update(OMP_NUM_THREADS="4", MKL_NUM_THREADS="4")
        try:
            subprocess.run(cmd, shell=True, check=True, text=True, env=env)
            if not seqtab.exists():
                append_with_flock(f"{acc}:DADA2_FAIL", FAILED_FILE)
                return
            # mark Dada2 done
            update_checkpoint(Path(os.environ["INPUT_FILE"]), acc, "Dada2", "1", os.environ["DELIM"])
            # mark Completed now that DADA2 really ran
            update_checkpoint(Path(os.environ["INPUT_FILE"]), acc, "Completed", "1", os.environ["DELIM"])
            logger.info(f"Finished DADA2 for {acc}")
        except subprocess.CalledProcessError:
            append_with_flock(f"{acc}:DADA2_FAIL", FAILED_FILE)
            return
    else:
        # we already had d2 == "1", so skip entirely
        logger.debug(f"Skipping {acc}: already has Dada2=1")


def merge_profiles(biop: str, env_name: str):
    """Merge DADA2 profiles for a bioproject."""
    odir = OUTPUT_BASE/biop
    with Path(os.environ["INPUT_FILE"]).open() as f:
        reader = csv.reader(f, delimiter=os.environ["DELIM"])
        next(reader)
        for row in reader:
            if row[0]==biop and row[2]=="16S" and row[7]!="1":
                return
    tabs = list(odir.glob("asv_*.rds"))
    if tabs:
        cmd = f"Rscript {Path.cwd()}/merge_dada2.R {biop} {' '.join(map(str,tabs))}"
        run_command(cmd, f"Merging {biop}", env_name)

def final_validation_and_merge(input_file: Path, delim: str, env_name: str):
    """Run final validation and merge profiles."""
    bioprojects = set()
    with input_file.open() as f:
        reader = csv.reader(f, delimiter=delim)
        next(reader)
        for row in reader:
            if row[2]=="16S":
                bioprojects.add(row[0])
    for bp in sorted(bioprojects):
        merge_profiles(bp, env_name)
    logger.info("Final validation done")

def process_samples(input_file: Path, delim: str, num_workers: int, env_name: str):
    """Process samples in parallel."""
    with input_file.open() as f:
        reader = csv.reader(f, delimiter=delim)
        next(reader)
        samples = [row for row in reader if row[2] == "16S"]
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_sample, *sample, env_name) for sample in samples]
        for f in as_completed(futures):
            f.result()  # will raise if any failed
    logger.info("All samples finished")

def main():
    args = parse_arguments()
    from pathlib import Path
    if args.submit_slurm:
        submit_slurm_job()
        return
    if args.debug:
        logger.setLevel(logging.DEBUG)
        
    # Cap num_workers at 3/4ths of available CPUs for non-DADA2 steps
    cpu_count = os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 8)
    if isinstance(cpu_count, str):
        cpu_count = int(cpu_count)
    if cpu_count is None:
        logger.warning("CPU count unavailable; falling back to 8 CPUs")
        cpu_count = 8
    args.num_workers = max(1, min(args.num_workers, cpu_count * 3 // 4))  # Use 3/4ths of CPUs
    logger.info(f"Using {args.num_workers} workers for non-DADA2 parallel steps ({args.num_workers} CPUs)")

    setup_directories()
    env_name = setup_environment()
    delim = validate_input_file(Path(args.input_file))
    os.environ["DELIM"] = delim
    os.environ["INPUT_FILE"] = args.input_file
    os.environ["FAILED_FILE"] = str(FAILED_FILE)
    os.environ["OUTPUT_BASE"] = str(OUTPUT_BASE)

    check_rdp_database()
    fastq_dir = Path(args.fastq_dir).resolve()

    # Step 1: Generate checksums if needed
    start = time.time()
    if not (fastq_dir / "checksums.b3").exists():
        generate_fastq_checksums(fastq_dir, args.num_workers)
    logger.info(f"Checksum generation took {time.time() - start:.2f} seconds")

    # Step 2: Build validated set (may be empty initially)
    start = time.time()
    validated_accessions = build_validated_set(Path(args.input_file), delim)
    logger.info(f"Validated set building took {time.time() - start:.2f} seconds")

    # Step 3: Repair FASTQs if needed
    start = time.time()
    if args.repair_fastqs:
        fastq_dir = repair_fastq_if_needed(fastq_dir, args.num_workers, validated_accessions)
    logger.info(f"FASTQ repair took {time.time() - start:.2f} seconds")

    # Step 4: Update input CSV with FASTQ paths + validation
    start = time.time()
    update_input_with_fastq_paths(Path(args.input_file), fastq_dir, delim, args.num_workers)
    logger.info(f"FASTQ path update took {time.time() - start:.2f} seconds")

    # NEW: Step 4.5 – Restore/checkpoint DADA2 results if output dir exists and is non-empty
    if OUTPUT_BASE.exists() and any(OUTPUT_BASE.iterdir()):
        logger.info(f"DADA2 output detected {time.time() - start:.2f} seconds")
        from pathlib import Path
        df = regenerate_input_with_validation(Path(args.input_file), OUTPUT_BASE, delim)
        df.to_csv(Path(args.input_file), sep=delim, index=False)
        logger.info(f"Checkpoint fields restored from existing DADA2 output {time.time() - start:.2f} seconds")

    # Re‑read the CSV you just overwrote
    full_df = pd.read_csv(args.input_file, sep=delim, dtype=str)
    
    # Filter to only those still needing DADA2
    to_run = full_df.loc[
        (full_df["SequencingType"] == "16S") & (full_df["Dada2"] == "0"),
        full_df.columns
    ]
    
    if to_run.empty:
        logger.info("All samples already have DADA2 results — nothing to do.")
    else:
        # Write a temporary CSV of just the remaining rows
        subset_path = Path(args.input_file).with_name("pending_DADA2.csv")
        to_run.to_csv(subset_path, sep=delim, index=False)
        logger.info(f"{len(to_run)} samples pending DADA2 → {subset_path.name}")
    
    # Step 5: Run DADA2 processing
    start = time.time()
    cpu_count = os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 8)
    if isinstance(cpu_count, str):
        cpu_count = int(cpu_count)
    if cpu_count is None:
        logger.warning("CPU count unavailable; falling back to 8 CPUs")
        cpu_count = 8
    args.num_workers = max(1, min(args.num_workers, cpu_count // 4, 8))
    logger.info(f"Using {args.num_workers} workers for DADA2 processing")
    process_samples(subset_path, delim, args.num_workers, env_name)
    logger.info(f"DADA2 processing took {time.time() - start:.2f} seconds")
    
    # Step 6: Final validation and merging
    start = time.time()
    final_validation_and_merge(Path(args.input_file), delim, env_name)
    logger.info(f"Final validation took {time.time() - start:.2f} seconds")

    logger.info("Pipeline complete.")


if __name__ == "__main__":
    if "--process-sample" in sys.argv:
        logger.error("process-sample mode not fully supported in Python version")
        sys.exit(1)
    main()
