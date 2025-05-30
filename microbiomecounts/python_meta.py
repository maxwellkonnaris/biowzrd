#!/usr/bin/env python3

import argparse, csv, gzip, logging, os, re, shutil, subprocess, sys, time, fcntl, tempfile
from datetime import datetime
from pathlib import Path
from typing import List, Set, Dict, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import shlex
from tqdm import tqdm

# set up logging
logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger(__name__)

SLURM_CONFIG = {
    "job-name": "counts_meta",
    "output": "slurm_meta-%j.out",
    "error": "slurm_meta-%j.err",
    "time": "48:00:00",
    "nodes": 1,
    "ntasks": 1,
    "cpus-per-task": 16,
    "mem": "256G",
    "account": "open",
    "mail-user": "mak6930@psu.edu",
}

# Default variables and directories
DEFAULT_DIR = "fastq_data/fastq_biologicaldata"
OUTPUT_BASE = Path("metagenomics")
LOCK_DIR = Path("locks")
FAILED_FILE = Path("failed_meta.log")
METAPHLAN_DB = "/storage/work/mak6930/applicationstorage/micromamba/envs/mpa/lib/python3.7/site-packages/metaphlan/metaphlan_databases"

# these must already be in df, or you’ll raise an error
HARD_REQUIRED = ["Bioproject", "RunAccession", "SequencingType"]

# these will be added if missing, with “Validated”→"0", fastq cols→""
OPTIONAL_COLUMNS = ["Fastq1", "Fastq2", "Validated", "MetaPhlAn", "Motus",  "Completed"]
ZERO_DEFAULTS = {"Validated", "MetaPhlAn", "Motus",  "Completed"}

REQUIRED_COLUMNS_META = HARD_REQUIRED + OPTIONAL_COLUMNS

def ensure_required_columns(df: pd.DataFrame,
                            required: List[str],
                            zero_cols: Set[str]) -> pd.DataFrame:
    """
    Guarantee that `df` has every column in `required`.
    For any missing col in zero_cols, fill with "0"; otherwise fill with "".
    """
    for col in required:
        if col not in df.columns:
            df[col] = "0" if col in zero_cols else ""
    return df

def prepare_meta_df(df: pd.DataFrame) -> pd.DataFrame:
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
    so multiple processes can't write simultaneously.
    """
    with file_path.open("a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(line + "\n")
        fcntl.flock(f, fcntl.LOCK_UN)
        
def get_unique_bioprojects(csv_path: Path, delim: str) -> List[str]:
    df = pd.read_csv(csv_path, sep=delim, usecols=["Bioproject"], dtype=str)
    return sorted(df["Bioproject"].dropna().unique())

def submit_slurm(cmd: List[str], job_name: str, dependency: Optional[str] = None) -> None:
    """
    Submit the given command as a SLURM job named `job_name`.
    If `dependency` is provided, the job will wait for it to finish successfully.
    """
    # Build a temporary sbatch script
    with tempfile.NamedTemporaryFile("w", suffix=".slurm", delete=False) as sb:
        sb.write("#!/bin/bash\n")
        # SBATCH directives
        for key, val in SLURM_CONFIG.items():
            sb.write(f"#SBATCH --{key}={val}\n")
        sb.write(f"#SBATCH --job-name={job_name}\n")
        sb.write(f"#SBATCH --output=slurm_{job_name}-%j.out\n")
        sb.write(f"#SBATCH --error=slurm_{job_name}-%j.err\n")
        if dependency:
            sb.write(f"#SBATCH --dependency=afterok:{dependency}\n")
        sb.write("\nsource ~/.bashrc\n")
        # Activate the current conda/micromamba environment
        active_env = os.getenv("CONDA_DEFAULT_ENV") or os.getenv("MAMBA_DEFAULT_ENV") or "mpa"
        sb.write(f"micromamba activate {active_env}\n\n")
        # Write the actual command
        sb.write(" ".join(shlex.quote(x) for x in cmd) + "\n")
        script_path = sb.name

    # Submit and clean up
    try:
        result = subprocess.run([
            "sbatch", script_path
        ], capture_output=True, text=True, check=True)
        logger.info(f"[SLURM] Submitted {job_name}: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        logger.error(f"SLURM submission failed for {job_name}: {e.stderr.strip()}")
        sys.exit(1)
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            logger.warning(f"Could not remove temporary sbatch script: {script_path}")


def submit_self_as_slurm(dependency: Optional[str] = None) -> None:
    """
    If not already running under SLURM, submit this script to Slurm.
    Removes the --submit-slurm flag to avoid recursion.
    """
    # Skip if already under SLURM
    if "SLURM_JOB_ID" in os.environ:
        logger.info("Already running under SLURM; skipping submission.")
        return

    # Build the command to re-launch without --submit-slurm
    cmd = [sys.executable] + [arg for arg in sys.argv if arg != "--submit-slurm"]
    # Derive a job name from the script filename
    job_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
    submit_slurm(cmd, job_name, dependency)
    sys.exit(0)

def parse_arguments():
    """Parse command-line arguments."""
    cpu_count = os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 8)
    if isinstance(cpu_count, str):
        cpu_count = int(cpu_count)
    if cpu_count is None:
        logger.warning("CPU count unavailable; falling back to 8 CPUs")
        cpu_count = 8
    default_workers = cpu_count * 3 // 4  # Default to 3/4ths of CPUs

    parser = argparse.ArgumentParser(description="Metagenomics FASTQ processing pipeline")
    parser.add_argument("-i", "--input-file", required=True, help="Input CSV file")
    parser.add_argument("-d", "--fastq-dir", default=DEFAULT_DIR, help="FASTQ directory")
    parser.add_argument("-w", "--num-workers", type=int, default=default_workers,
                       help=f"Number of workers (default: 3/4ths of CPUs = {default_workers})")
    parser.add_argument("--mode",
                        choices=["both", "metaphlan", "motus"],
                        default="metaphlan",
                        help="Run only MetaPhlAn, only mOTUs, or both")
    parser.add_argument("--metaphlan-threads", type=int, default=4,
                        help="Threads for MetaPhlAn bowtie2 (default: 4)")
    parser.add_argument("--motus-threads",    type=int, default=4,
                        help="Threads for mOTUs profiling (default: 4)")
    parser.add_argument("--skip-preprocessing", action="store_true",
                        help="Skip checksum generation, FASTQ repair, path "
                             "updates, and checkpoint regeneration. "
                             "Assumes the input CSV already contains valid "
                             "Fastq1/Fastq2/Validated and MetaPhlAn/mOTUs "
                             "checkpoint columns.")
    parser.add_argument("--bioproject",
                        help="Run only this BioProject ID. "
                             "If omitted the script becomes a *launcher* that "
                             "fires one SLURM job per BioProject.")
    parser.add_argument("--no-repair-fastqs", action="store_false", dest="repair_fastqs",
                       help="Disable FASTQ repair")
    parser.add_argument("--submit-slurm", action="store_true",
                        help="Submit this script as a SLURM job and exit")
    parser.add_argument("--dependency", type=str, default=None,
                        help="SLURM job ID to wait for (afterok)")
    parser.add_argument("--debug", action="store_false",
                        help="add logging verbose")
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
        logger.info(f"Using activated MPA environment: {active_env}")
        return active_env
    except subprocess.CalledProcessError:
        logger.error("micromamba not in PATH")
        sys.exit(1)

def run_command(cmd: str, description: str, env_name: str = None) -> None:
    """
    Run a shell command (e.g. Rscript) directly, assuming the environment is already active.
    Captures stdout and stderr, and logs them on failure.
    """
    logger.info(description)
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        return result.stdout  # if callers expect stdout
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {description}")
        logger.error(f"[stdout]\n{e.stdout.strip()}")
        logger.error(f"[stderr]\n{e.stderr.strip()}")
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

def check_metaphlan_database() -> Path:
    """Detect or install MetaPhlAn database with fallback."""
    try:
        result = subprocess.run(
            ["python", "-c", "import metaphlan, pathlib; print(pathlib.Path(metaphlan.__file__).parent / 'metaphlan_databases')"],
            capture_output=True, text=True, check=True
        )
        db_path = Path(result.stdout.strip())
    except subprocess.CalledProcessError:
        logger.debug(f"MetaPhlAn DB auto-detection failed; using fallback: {METAPHLAN_DB}")
        db_path = Path(METAPHLAN_DB)

    if not db_path.is_dir():
        logger.error(f"MetaPhlAn DB not found at {db_path}")
        sys.exit(1)

    if not list(db_path.glob("*.bt2l")):
        logger.info("MetaPhlAn DB empty – installing...")
        run_command(f"metaphlan --install --bowtie2db \"{db_path}\"", "Installing MetaPhlAn DB")
    else:
        logger.info(f"MetaPhlAn DB found at {db_path}")

    return db_path

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
    validated_accessions: List[str]
) -> Path:
    """Repair corrupted FASTQ files and return new FASTQ dir."""
    repaired_dir = fastq_dir.parent / f"{fastq_dir.name}_repaired"
    checksum_file = fastq_dir / "checksums.b3"
    failed_log    = fastq_dir / "failed_checksums.txt"

    repaired_dir.mkdir(exist_ok=True)
    failed_log.write_text("")

    # 1) Ensure we have an initial checksum list
    if not checksum_file.exists():
        generate_fastq_checksums(fastq_dir, num_workers)

    # 2) Build list of (hash, filename) to check, skipping already-validated
    skip = {f"{acc}.fastq.gz" for acc in validated_accessions}
    to_check: List[Tuple[str, str]] = []
    with checksum_file.open() as f:
        for line in f:
            hash_val, fname = line.strip().split(maxsplit=1)
            if Path(fname).name not in skip:
                to_check.append((hash_val, fname))

    # 3) Verify checksums in parallel, collecting any failures
    def verify(h: str, fn: str) -> Optional[str]:
        cmd = f"echo '{h}  {fn}' | b3sum -c - --quiet"
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            return None
        except subprocess.CalledProcessError:
            return fn

    failed_files: List[str] = []
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

def _check_file_valid(path: Path, min_size: int = 100) -> bool:
    """Return True if file exists and is at least `min_size` bytes."""
    return path.exists() and path.stat().st_size >= min_size

def regenerate_input_with_validation(
    input_file: Path,
    output_dir: Path,
    delim: str = ","
) -> pd.DataFrame:
    """
    Read a CSV mapping and update MetaPhlAn/mOTUs status columns.
    """
    df = pd.read_csv(input_file, sep=delim, dtype=str)
    df["MetaPhlAn"] = "0"
    df["Motus"] = "0"
    df["Completed"] = "0"

    # Check merged outputs
    for idx, row in df.iterrows():
        biop = row["Bioproject"]
        mp_counts_merge = output_dir / biop / f"{biop}_MetaPhlAn_counts_merged.tsv"
        mp_abund_merge = output_dir / biop / f"{biop}_MetaPhlAn_abundances_merged.tsv"
        mt_merge = output_dir / biop / f"{biop}_motus_merged.tsv"
        done_mp = _check_file_valid(mp_counts_merge) and _check_file_valid(mp_abund_merge)
        done_mt = _check_file_valid(mt_merge)
        df.at[idx, "MetaPhlAn"] = "1" if done_mp else "0"
        df.at[idx, "Motus"] = "1" if done_mt else "0"
        if done_mp and done_mt:
            df.at[idx, "Completed"] = "1"

    # Check per-accession outputs
    tasks = []
    for idx, row in df.iterrows():
        biop, acc = row["Bioproject"], row["RunAccession"]
        if df.at[idx, "MetaPhlAn"] == "0":
            tasks.append((idx, output_dir / biop / f"{acc}_MetaPhlAn_counts.txt", "MetaPhlAn"))
            tasks.append((idx, output_dir / biop / f"{acc}_profiled.tsv", "MetaPhlAn"))
        if df.at[idx, "Motus"] == "0":
            tasks.append((idx, output_dir / biop / f"{acc}_motus.out", "Motus"))

    with ProcessPoolExecutor(max_workers=os.cpu_count() // 2 or 4) as exe:
        futures = {exe.submit(_check_file_valid, path): (idx, method) for idx, path, method in tasks}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Verifying per-accession outputs"):
            idx, method = futures[future]
            if future.result():
                df.at[idx, method] = "1"

    df["Completed"] = df.apply(
        lambda r: "1" if r["MetaPhlAn"] == "1" and r["Motus"] == "1" else "0",
        axis=1
    )

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    bak_path = input_file.with_suffix(f".bak.{timestamp}.csv")
    df.to_csv(bak_path, sep=delim, index=False)
    logger.info(f"[BACKUP] {bak_path}")
    df.to_csv(input_file, sep=delim, index=False)
    logger.info(f"[REGEN] Updated {input_file} with MetaPhlAn/mOTUs checkpoints")
    return df

def update_input_with_fastq_paths(
    input_file: Path,
    fastq_dir: Path,
    delim: str,
    num_workers: int
):
    logger.info("Updating FASTQ paths")

    # 1) Load input CSV
    df = pd.read_csv(input_file, sep=delim, dtype=str)
    df = prepare_meta_df(df)
    input_accs = set(df.loc[df["SequencingType"] == "meta", "RunAccession"].str.strip())
    logger.info(f"{len(input_accs)} target accessions from CSV")

    # 2) Scan FASTQ dir
    # sanity-check the FASTQ directory up front
    if not fastq_dir.is_dir():
        logger.error(f"FASTQ directory not found: {fastq_dir}")
        sys.exit(1)
    fastq_files = list(fastq_dir.glob("*.fastq.gz"))
    if not fastq_files:
        logger.error(f"No .fastq.gz files found in: {fastq_dir}")
        sys.exit(1)
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

    df = ensure_required_columns(df, REQUIRED_COLUMNS_META,ZERO_DEFAULTS)

    # And finally overwrite the CSV
    logger.info(f"Writing updated CSV → {input_file}")
    df[REQUIRED_COLUMNS_META].to_csv(input_file, sep=delim, index=False)
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
    so multiple processes can't clobber each other.
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

def convert_metaphlan_to_counts(
    profile: Path,
    output: Path
) -> bool:
    """
    Convert a MetaPhlAn profile (with commented header lines) into a simple
    clade→relative_abundance→read count TSV.

    - Skips all lines starting with '#'
    - Expects columns: clade_name, NCBI_tax_id, relative_abundance, [additional_species]
    - Uses relative_abundance (3rd field) and total mapped reads from the log.

    Returns True on success (and writes `output`), False on any error.
    """
    # 1) parse total reads from the log
    mapped = None
    with profile.open() as f:
        for line in f:
            if "reads processed" in line.lower():
                try:
                    # e.g. "#5420638 reads processed"
                    mapped = float(line.strip().lstrip("#").split()[0])
                except ValueError:
                    logger.error(f"Could not parse reads from {profile}: {line.strip()}")
                    return False
                break
    if mapped is None:
        logger.error(f"Total reads not found in {profile}")
        return False

    # 2) open profile, skip comments, parse data lines
    try:
        with profile.open() as inf, output.open("w") as outf:
            outf.write("#clade\tabundance\treads\n")
            for line in inf:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 3:
                    logger.warning(f"Skipping malformed line in {profile}: {line.strip()}")
                    continue
                clade     = fields[0]
                abundance = float(fields[2])  # third column
                reads     = mapped * abundance / 100.0
                outf.write(f"{clade}\t{abundance:.6f}\t{reads:.2f}\n")
        logger.debug(f"Converted MetaPhlAn profile → counts at {output}")
        return True

    except Exception as e:
        logger.error(f"Error converting {profile} to counts: {e}")
        return False

def merge_profiles(biop: str, tool: str):
    """Merge MetaPhlAn or mOTUs profiles for a bioproject."""
    odir = OUTPUT_BASE / biop

    # locate MetaPhlAn merge script
    try:
        import metaphlan
        merge_script = Path(metaphlan.__file__).parent / "utils" / "merge_metaphlan_tables.py"
        if not merge_script.exists():
            logger.error(f"MetaPhlAn merge script not found: {merge_script}")
            return
    except ImportError:
        logger.error("Cannot import 'metaphlan'; skipping merge step.")
        return

    if tool == "metaphlan":
        BATCH_SIZE = 95

        # --- PROPORTIONS ---
        prop_files = sorted(odir.glob("*_profiled.tsv"))
        if prop_files:
            prop_out = odir / f"{biop}_proportions.tsv"
            batches = []
            for i in range(0, len(prop_files), BATCH_SIZE):
                batch = prop_files[i : i + BATCH_SIZE]
                num = i // BATCH_SIZE + 1
                batch_out = odir / f"{biop}_proportions_batch{num}.tsv"

                if len(batch) == 1:
                    shutil.copy2(batch[0], batch_out)
                else:
                    files_str = " ".join(f'"{p}"' for p in batch)
                    cmd = f'python "{merge_script}" {files_str} > "{batch_out}"'
                    run_command(cmd, f"Merging proportions batch {num} for {biop}")

                batches.append(batch_out)

            # final merge or move
            if len(batches) == 1:
                shutil.move(str(batches[0]), str(prop_out))
            else:
                files_str = " ".join(f'"{p}"' for p in batches)
                final_cmd = f'python "{merge_script}" {files_str} > "{prop_out}"'
                run_command(final_cmd, f"Final merge of proportions for {biop}")

            # — show first 10 lines of the merged proportions —
            logger.info(f"First 10 lines of {prop_out.name}:")
            with open(prop_out, 'r') as fh:
                for idx, line in enumerate(fh):
                    if idx >= 10:
                        break
                    logger.info(f"  {line.rstrip()}")

        else:
            logger.debug(f"No MetaPhlAn proportions files found for {biop}")

        # --- COUNTS (pandas concat) ---
        count_files = sorted(odir.glob("*_MetaPhlAn_counts.txt"))
        if count_files:
            count_out = odir / f"{biop}_MetaPhlAn_merged_counts.tsv"
            series_list = []
            for f in count_files:
                sample = f.stem.replace("_MetaPhlAn_counts", "")
                df = pd.read_csv(
                    f,
                    sep="\t",
                    comment="#",
                    header=None,
                    names=["clade", "abundance", "reads"],
                    dtype={"clade": str, "abundance": float, "reads": float},
                )
                df = df.set_index("clade")
                reads = df["reads"].rename(sample)
                series_list.append(reads)

            merged = pd.concat(series_list, axis=1).fillna(0).astype(int)
            merged.to_csv(count_out, sep="\t")

            # — show first 10 lines of the merged counts —
            logger.info(f"First 10 lines of {count_out.name}:")
            with open(count_out, 'r') as fh:
                for idx, line in enumerate(fh):
                    if idx >= 10:
                        break
                    logger.info(f"  {line.rstrip()}")

        else:
            logger.debug(f"No MetaPhlAn counts files found for {biop}")

    else:  # motus
        files = list(odir.glob("*_motus.out"))
        if not files:
            logger.debug(f"No mOTUs files found for {biop}")
            return
        out = odir / f"{biop}_motus_merged.tsv"
        file_list = ",".join(map(str, files))
        cmd = f"motus merge -i \"{file_list}\" -o \"{out}\""
        run_command(cmd, f"Merging {biop} {tool}")

        # — show first 10 lines of the mOTUs merged file —
        logger.info(f"First 10 lines of {out.name}:")
        with open(out, 'r') as fh:
            for idx, line in enumerate(fh):
                if idx >= 10:
                    break
                logger.info(f"  {line.rstrip()}")


def final_validation_and_merge(input_file: Path, delim: str):
    """Run final validation and merge MetaPhlAn and mOTUs profiles."""
    bioprojects = set()
    with input_file.open() as f:
        reader = csv.reader(f, delimiter=delim)
        next(reader) 
        for row in reader:
            if row[2] == "meta":
                bioprojects.add(row[0])
    for bp in sorted(bioprojects):
        merge_profiles(bp, "metaphlan")
        merge_profiles(bp, "motus")
    logger.info("Final merge of MetaPhlAn counts, abundances, and mOTUs done")

def process_sample(
    biop: str,
    acc: str,
    st: str,
    fq1: str,
    fq2: str,
    val: str,
    mp_flag: str,
    mt_flag: str,
    comp: str,
    mode: str,
    metaphlan_threads,
    motus_threads
):
    from time import perf_counter

    # skip non-meta or already‐done samples
    if st != "meta" or val == "0" or comp == "1":
        return

    outdir = OUTPUT_BASE / biop
    outdir.mkdir(exist_ok=True, parents=True)
    
    input_file = Path(os.environ["INPUT_FILE"])
    delim = os.environ["DELIM"]
    extra_mp_args = "--read_min_len 30"

    # -------------------------------------------------------------------
    # MetaPhlAn
    # -------------------------------------------------------------------
    if mode in ("both", "metaphlan") and mp_flag != "1":
        start = perf_counter()
        mp_out        = outdir / f"{acc}_profiled.tsv"
        mp_counts     = outdir / f"{acc}_MetaPhlAn_counts.txt"
        mp_bowtie2out = outdir / f"{acc}_bowtie2out.txt"
    
        # decide whether to reuse an existing bowtie2out
        use_bowtie2out = mp_bowtie2out.exists() and mp_bowtie2out.stat().st_size > 1000
    
        if use_bowtie2out:
            cmd = (
                f"metaphlan {mp_bowtie2out} "
                f"--input_type bowtie2out "
                f"--nproc {metaphlan_threads} "
                f"--bowtie2db {os.environ['METAPHLAN_DB']} "
                f"--unclassified_estimation "
                f"{extra_mp_args} " 
                f"-o {mp_out}"
            )
        else:
            cmd = (
                f"metaphlan {fq1} "
                f"--input_type fastq "
                f"--nproc {metaphlan_threads} "
                f"--bowtie2db {os.environ['METAPHLAN_DB']} "
                f"--bowtie2out {mp_bowtie2out} "
                f"--unclassified_estimation "
                f"{extra_mp_args} " 
                f"-o {mp_out}"
            )
    
        logger.info(f"{acc}: MetaPhlAn → start")
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
            # sanity‐check output
            if not mp_out.exists() or mp_out.stat().st_size < 100:
                raise RuntimeError("MetaPhlAn output missing or too small")
            # convert to counts
            if not convert_metaphlan_to_counts(mp_out, mp_counts):
                raise RuntimeError("MetaPhlAn counts conversion failed")
    
            logger.info(f"{acc}: MetaPhlAn → success ({perf_counter()-start:.1f}s)")
            update_checkpoint(input_file, acc, "MetaPhlAn", "1", delim)
    
        except subprocess.CalledProcessError as e:
            logger.info(f"{acc}: MetaPhlAn → FAIL (exit {e.returncode})")
            logger.error(f"[MetaPhlAn stderr]\n{e.stderr.strip()}")
            append_with_flock(f"{acc}:MP_FAIL", FAILED_FILE)
            return
    
        except Exception as e:
            logger.info(f"{acc}: MetaPhlAn → FAIL ({e})")
            logger.exception("Unexpected error in MetaPhlAn block")
            append_with_flock(f"{acc}:MP_FAIL", FAILED_FILE)
            return


    # -------------------------------------------------------------------
    # mOTUs
    # -------------------------------------------------------------------
    if mode in ("both", "motus") and mt_flag != "1":
        start = perf_counter()
        mt_out = outdir / f"{acc}_motus.out"
        cmd = (
            f"motus profile -s '{fq1}' "
            f"-c "
            f"-t {motus_threads} "
            f"-o '{mt_out}' "
            f"-A "
        )
    
        logger.info(f"{acc}: mOTUs → start")
        try:
            subprocess.run(cmd, shell=True, check=True, text=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if not mt_out.exists() or mt_out.stat().st_size < 100:
                raise RuntimeError("mOTUs output missing/too small")
            logger.info(f"{acc}: mOTUs → success ({perf_counter()-start:.1f}s)")
            update_checkpoint(input_file, acc, "Motus", "1", delim)
        except Exception as e:
            logger.info(f"{acc}: mOTUs → FAIL ({e})")
            append_with_flock(f"{acc}:MT_FAIL", FAILED_FILE)
            return


    # -------------------------------------------------------------------
    # Mark overall completion
    # -------------------------------------------------------------------
    update_checkpoint(input_file, acc, "Completed", "1", delim)
    logger.info(f"{acc}: Completed")
    
def process_samples(input_file: Path, delim: str, num_workers: int, mode: str,
    metaphlan_threads: str,
    motus_threads: str):
    """Process samples in parallel."""
    with input_file.open() as f:
        reader = csv.reader(f, delimiter=delim)
        header = next(reader)  # Skip header
        expected_columns = len(REQUIRED_COLUMNS_META)  # Should be 9
        if len(header) != expected_columns:
            logger.error(f"CSV header has {len(header)} columns, expected {expected_columns}: {header}")
            sys.exit(1)
        samples = [row for row in reader if row[2] == "meta"]
        for sample in samples:
            if len(sample) != expected_columns:
                logger.error(f"Row has {len(sample)} columns, expected {expected_columns}: {sample}")
                sys.exit(1)
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_sample, *sample, mode, metaphlan_threads, motus_threads) for sample in samples]
        for f in as_completed(futures):
            f.result()  # Will raise if any failed
    logger.info("All samples finished")

def main():
    args = parse_arguments()

    # 1) If flagged for Slurm submission, self-submit and exit
    if args.submit_slurm:
        submit_self_as_slurm(dependency=args.dependency)

    # 2) Launcher mode: no --bioproject -> fire one job per BioProject
    if args.bioproject is None:
        delim = validate_input_file(Path(args.input_file))
        bps = get_unique_bioprojects(Path(args.input_file), delim)
        logger.info(f"Launching {len(bps)} BioProjects as separate SLURM jobs")

        base_cmd = [sys.executable, sys.argv[0]] + [
            a for a in sys.argv[1:]
            if a != "--submit-slurm" and not a.startswith("--bioproject")
        ]
        for bp in bps:
            cmd = base_cmd + ["--bioproject", bp, "--submit-slurm"]
            submit_slurm(cmd, job_name=f"{bp}_meta", dependency=args.dependency)
        return

    # 3) Worker mode: continue locally under Slurm (or without Slurm)
    if args.debug:
        logger.setLevel(logging.DEBUG)

    setup_directories()
    env_name = setup_environment()
    delim = validate_input_file(Path(args.input_file))
    os.environ.update({
        "DELIM": delim,
        "INPUT_FILE": args.input_file,
        "FAILED_FILE": str(FAILED_FILE),
        "OUTPUT_BASE": str(OUTPUT_BASE),
    })
    db = check_metaphlan_database()
    os.environ["METAPHLAN_DB"] = str(db)
    fastq_dir = Path(args.fastq_dir).resolve()

    # 4) Pre-processing (checksums, repair, CSV updates)
    if not args.skip_preprocessing:
        # 4a) checksums
        if not (fastq_dir / "checksums.b3").exists():
            generate_fastq_checksums(fastq_dir, args.num_workers)
        # 4b) build validated set
        validated = build_validated_set(Path(args.input_file), delim)
        # 4c) repair FASTQs
        if args.repair_fastqs:
            fastq_dir = repair_fastq_if_needed(fastq_dir, args.num_workers, validated)
    else:
        logger.info("⚡ Skipping preprocessing (--skip-preprocessing)")
        
    # 4d) update CSV with Fastq1/Fastq2/Validated
    update_input_with_fastq_paths(Path(args.input_file), fastq_dir, delim, args.num_workers)
    
    # 4e) regenerate MetaPhlAn/mOTUs checkpoints if outputs already exist
    if OUTPUT_BASE.exists() and any(OUTPUT_BASE.iterdir()):
        regenerate_input_with_validation(Path(args.input_file), OUTPUT_BASE, delim)

    # 5) Load CSV, filter to this BioProject (if any), select pending samples
    full_df = pd.read_csv(args.input_file, sep=delim, dtype=str)
    full_df = prepare_meta_df(full_df) 
    if args.bioproject:
        full_df = full_df[full_df["Bioproject"] == args.bioproject]

    to_run = full_df.loc[
        (full_df["SequencingType"] == "meta") &
        ((full_df["MetaPhlAn"] == "0") | (full_df["Motus"] == "0"))
    ]

    if to_run.empty:
        logger.info("All samples already processed—nothing to do.")
        return

    subset_name = f"pending_{args.bioproject or 'all'}.csv"
    subset_path = Path(args.input_file).with_name(subset_name)
    to_run.to_csv(subset_path, sep=delim, index=False)
    logger.info(f"{len(to_run)} samples pending → {subset_name}")

    # 6) Profile in parallel
    cpu_env    = os.environ.get("SLURM_CPUS_PER_TASK")
    cpu_count  = int(cpu_env) if cpu_env and cpu_env.isdigit() else (os.cpu_count() or 4)
    args.num_workers = max(1, min(args.num_workers, cpu_count // 4, 4))
    logger.info(f"Using {args.num_workers} workers for profiling")
    
    try:
        process_samples(
            subset_path,
            delim,
            args.num_workers,
            args.mode,
            args.metaphlan_threads,
            args.motus_threads
        )
    except Exception as e:
        logger.error("Error during sample profiling!", exc_info=True)
        # Optionally: sys.exit(1) if you want to stop here
        return


    # 7) Final merge: for a single BioProject, merge just that one;
    #    otherwise fall back to the old “merge all” helper
    if args.bioproject:
        merge_profiles(args.bioproject, "metaphlan")
        merge_profiles(args.bioproject, "motus")
    else:
        final_validation_and_merge(Path(args.input_file), delim)

    logger.info("Pipeline complete.")


if __name__ == "__main__":
    if "--process-sample" in sys.argv:
        logger.error("process-sample mode not fully supported in Python version")
        sys.exit(1)
    main()
