#!/usr/bin/env python3
# utils_counts.py

import argparse
import csv
import fcntl
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

try:
    import fireducks.pandas as pd
    logger = logging.getLogger(__name__)
    logger.info("Using fireducks.pandas for faster multithreaded I/O")
except ImportError:
    import pandas as pd
    logger = logging.getLogger(__name__)
    logger.info("Falling back to pandas for single-threaded I/O")

from tqdm import tqdm

def setup_logging(level: int) -> None:
    """Configure logging with the specified level."""
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )

def compute_worker_counts(cpu_count: Optional[int] = None, processing_step: bool = False, seq_type: str = "16S") -> int:
    """Compute number of workers based on CPU count and step."""
    if cpu_count is None:
        cpu_count = os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 8)
        if isinstance(cpu_count, str):
            cpu_count = int(cpu_count)
        if cpu_count is None:
            logger.warning("CPU count unavailable; falling back to 8 CPUs")
            cpu_count = 8
    if processing_step:
        if seq_type == "16S":
            return max(1, min(cpu_count // 4, 8))  # DADA2: 4 CPUs per worker, max 8 workers
        else:  # Meta
            return max(1, min(cpu_count // 2, 16))  # MetaPhlAn/mOTUs: adjust as needed
    return max(1, cpu_count * 3 // 4)  # Non-processing: 3/4ths of CPUs

def submit_slurm_job(slurm_config: Dict[str, str], conda_env: str, args: argparse.Namespace) -> None:
    """Submit a SLURM job, re-invoking the script with --submit-slurm."""
    cmd = [sys.executable, sys.argv[0], f"-i={args.input_file}"]
    if args.fastq_dir != slurm_config.get("fastq_dir", "fastq_data/fastq_biologicaldata"):
        cmd.extend([f"-d={args.fastq_dir}"])
    if args.debug:
        cmd.append("--debug")
    if not args.repair_fastqs:
        cmd.append("--no-repair-fastqs")
    if hasattr(args, "env") and args.env != conda_env:
        cmd.extend([f"--env={args.env}"])
    if not args.submit_slurm:
        cmd.append("--submit-slurm")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".slurm", delete=False) as sb:
        sb.write("#!/bin/bash\n")
        for key, val in slurm_config.items():
            sb.write(f"#SBATCH --{key}={val}\n")
        sb.write("\n")
        sb.write(f"micromamba activate {conda_env}\n")
        sb.write("\n")
        sb.write(" ".join(cmd) + "\n")
        script_path = sb.name

    try:
        res = subprocess.run(
            ["sbatch", script_path], capture_output=True, text=True, check=True
        )
        logger.info(f"Submitted SLURM job: {res.stdout.strip()}")
        Path(script_path).unlink()
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed to submit SLURM job: {e.stderr}")
        raise

def ensure_required_columns(df: pd.DataFrame, required: List[str], zero_cols: Set[str]) -> pd.DataFrame:
    """Guarantee that `df` has every column in `required`."""
    for col in required:
        if col not in df.columns:
            df[col] = "0" if col in zero_cols else ""
    return df

def setup_directories(dirs: List[Path]) -> None:
    """Create required directories."""
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)

def append_with_flock(line: str, file_path: Path) -> None:
    """Append a line to a file with an exclusive lock."""
    with file_path.open("a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(line + "\n")
        fcntl.flock(f, fcntl.LOCK_UN)

def validate_input_file(file_path: Path) -> str:
    """Validate CSV header and determine delimiter."""
    with file_path.open("r") as f:
        first_line = f.readline().strip()
    for delim in [",", "\t", ";"]:
        if len(first_line.split(delim)) > 1:
            df = pd.read_csv(file_path, sep=delim, nrows=1)
            if all(col in df.columns for col in ["Bioproject", "RunAccession", "SequencingType"]):
                return delim
    logger.error("Malformed header in input file")
    sys.exit(1)

def setup_environment(conda_env: str) -> str:
    """Check micromamba and active environment."""
    try:
        result = subprocess.run(
            ["micromamba", "env", "list"], capture_output=True, text=True, check=True
        )
        envs = result.stdout.splitlines()
        active_env = next((line.split()[0] for line in envs if "*" in line), None)
        if active_env != conda_env:
            logger.error(f"Expected micromamba environment {conda_env}, but {active_env} is active")
            sys.exit(1)
        logger.info(f"Using activated micromamba environment: {active_env}")
        return active_env
    except subprocess.CalledProcessError:
        logger.error("micromamba not in PATH")
        sys.exit(1)

def check_database(database: Path) -> None:
    """Verify database exists."""
    if not database.is_file():
        logger.error(f"Database {database} not found")
        sys.exit(1)
    logger.info(f"Database {database} verified")

def process_fastq_global(fq_path_str: str, input_accs: Set[str]) -> Optional[Tuple[str, str]]:
    """Map FASTQ files to accessions."""
    fq = Path(fq_path_str)
    acc = re.sub(r"(_[12])?\.fastq(\.gz)?$", "", fq.name)
    return (acc, fq_path_str) if acc in input_accs else None

def generate_fastq_checksums(fastq_dir: Path, num_workers: int) -> None:
    """Generate BLAKE3 checksums for FASTQ files."""
    fastq_files = sorted(fastq_dir.glob("*.fastq.gz"))
    checksum_file = fastq_dir / "checksums.b3"
    checksums = []

    def compute_checksum(f: Path) -> str:
        try:
            result = subprocess.run(
                ["b3sum", str(f)], capture_output=True, text=True, check=True
            )
            return result.stdout.strip()
        except subprocess.CalledProcessError:
            logger.error(f"Failed to compute checksum for {f}")
            return ""

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        checksums = list(tqdm(
            executor.map(compute_checksum, fastq_files),
            total=len(fastq_files),
            desc="Checksumming FASTQs"
        ))
    with checksum_file.open("w") as f:
        f.write("\n".join([c for c in checksums if c]) + "\n")
    logger.info(f"Checksums saved to {checksum_file}")

def validate_fastq(fq: Path) -> bool:
    """Validate FASTQ file structure."""
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

def is_biological_fastq(fq: Path, seq_type: str = "16S") -> bool:
    """Check if FASTQ contains biological sequences."""
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
    if seq_type == "16S":
        return uniq_lengths > 1 or max_length > 29
    else:  # Meta
        return max_length > 50

def repair_fastq_if_needed(fastq_dir: Path, num_workers: int, validated_accessions: List[str]) -> Path:
    """Repair corrupted FASTQ files and return new FASTQ dir."""
    repaired_dir = fastq_dir.parent / f"{fastq_dir.name}_repaired"
    checksum_file = fastq_dir / "checksums.b3"
    failed_log = fastq_dir / "failed_checksums.txt"

    repaired_dir.mkdir(exist_ok=True)
    failed_log.write_text("")

    if not checksum_file.exists():
        generate_fastq_checksums(fastq_dir, num_workers)

    skip = {f"{acc}.fastq.gz" for acc in validated_accessions}
    to_check: List[Tuple[str, str]] = []
    with checksum_file.open() as f:
        for line in f:
            if line.strip():
                hash_val, fname = line.strip().split(maxsplit=1)
                if Path(fname).name not in skip:
                    to_check.append((hash_val, fname))

    def verify_checksum(hash_val: str, fname: str) -> Optional[str]:
        cmd = f"echo '{hash_val}  {fname}' | b3sum -c - --quiet"
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            return None
        except subprocess.CalledProcessError:
            return fname

    failed_files: List[str] = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(verify_checksum, h, fn) for h, fn in to_check]
        for future in tqdm(as_completed(futures), total=len(futures), desc="Verifying FASTQs"):
            bad = future.result()
            if bad:
                failed_files.append(bad)

    failed_log.write_text("\n".join(failed_files) + "\n")
    logger.info(f"{len(failed_files)} files failed checksum")

    if failed_files:
        logger.info("Repairing corrupted FASTQ files")
        def do_repair(fn: str):
            inp = fastq_dir / fn
            out = repaired_dir / inp.name
            try:
                with gzip.open(inp, "rt") as inf, gzip.open(out, "wt") as outf:
                    for line in inf:
                        outf.write(line.rstrip("\r") + "\n")
            except Exception as e:
                logger.error(f"Failed to repair {inp}: {e}")

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            list(tqdm(executor.map(do_repair, failed_files), total=len(failed_files), desc="Repairing FASTQs"))

    for fq in fastq_dir.glob("*.fastq.gz"):
        dest = repaired_dir / fq.name
        if not dest.exists():
            dest.symlink_to(fq.resolve())

    generate_fastq_checksums(repaired_dir, num_workers)
    logger.info(f"Repair complete → FASTQ_DIR updated to {repaired_dir}")
    return repaired_dir

def validate_and_check(fq: Path, acc: str, seq_type: str = "16S") -> Tuple[str, Optional[Path], Optional[str]]:
    """Validate FASTQ and check if biological."""
    if validate_fastq(fq):
        if is_biological_fastq(fq, seq_type):
            return acc, fq, None
        else:
            return acc, None, "TECHNICAL"
    else:
        return acc, None, "CORRUPT"

def update_input_with_fastq_paths(
    input_file: Path,
    fastq_dir: Path,
    delim: str,
    num_workers: int,
    seq_type: str,
    output_cols: List[str]
) -> None:
    """Update input CSV with FASTQ paths and validation status."""
    logger.info("Updating FASTQ paths")
    df = pd.read_csv(input_file, sep=delim, dtype=str)
    input_accs = set(df.loc[df["SequencingType"] == seq_type, "RunAccession"].str.strip())
    logger.info(f"{len(input_accs)} target accessions from CSV")

    fastq_map: Dict[str, List[Path]] = {}
    fastq_files = list(fastq_dir.glob("*.fastq.gz"))
    fq_paths = [str(p) for p in fastq_files]

    with ProcessPoolExecutor(max_workers=num_workers) as exe:
        futures = [exe.submit(process_fastq_global, fq, input_accs) for fq in fq_paths]
        for future in as_completed(futures):
            result = future.result()
            if result:
                acc, fq_path = result
                fq = Path(fq_path)
                fastq_map.setdefault(acc, []).append(fq)
    logger.info(f"Found FASTQs for {len(fastq_map)} of {len(input_accs)} accessions")

    validated_map: Dict[str, List[Path]] = {}
    with ProcessPoolExecutor(max_workers=num_workers) as exe:
        futures = [exe.submit(validate_and_check, fq, acc, seq_type) for acc, fqs in fastq_map.items() for fq in fqs]
        for future in as_completed(futures):
            acc, fq, reason = future.result()
            if reason:
                append_with_flock(f"{acc}:{reason}", fastq_dir / f"failed_{seq_type.lower()}.log")
            elif fq:
                validated_map.setdefault(acc, []).append(fq)

    validated_data = []
    for acc, vfs in validated_map.items():
        if not vfs:
            validated_data.append({"RunAccession": acc, output_cols[0]: "MISSING", output_cols[1]: "MISSING", output_cols[2]: "0"})
        elif len(vfs) == 1:
            validated_data.append({"RunAccession": acc, output_cols[0]: str(vfs[0]), output_cols[1]: "", output_cols[2]: "1"})
        else:
            validated_data.append({"RunAccession": acc, output_cols[0]: str(vfs[0]), output_cols[1]: str(vfs[1]), output_cols[2]: "1"})
    validated_df = pd.DataFrame(validated_data)

    df = df.drop(columns=output_cols, errors="ignore")
    if not validated_df.empty:
        df = df.merge(validated_df, on="RunAccession", how="left")
    else:
        for col in output_cols:
            df[col] = "MISSING" if col.startswith("Fastq") else "0"

    df[output_cols[0]] = df[output_cols[0]].fillna("MISSING")
    df[output_cols[1]] = df[output_cols[1]].fillna("MISSING")
    df[output_cols[2]] = df[output_cols[2]].fillna("0")

    df.to_csv(input_file, sep=delim, index=False)
    logger.info(f"FASTQ path update done → {input_file}")

def regenerate_flags(
    input_file: Path,
    output_base: Path,
    delim: str,
    specs: Dict[str, Dict],
    seq_type: str
) -> pd.DataFrame:
    """Regenerate flags based on output files."""
    df = pd.read_csv(input_file, sep=delim, dtype=str).fillna("")
    for spec in specs.values():
        for col in spec["columns"]:
            df[col] = "0"

    for idx, row in df.iterrows():
        if row["SequencingType"] != seq_type:
            continue
        biop, acc = row["Bioproject"], row["RunAccession"]
        bpdir = output_base / biop
        for tool, s in specs.items():
            merged = bpdir / s["merged"].format(biop=biop, acc=acc)
            pat, minsz = s["partial"]
            partial = bpdir / pat.format(biop=biop, acc=acc)
            if merged.exists():
                for col, val in zip(s["columns"], s["complete"]):
                    df.at[idx, col] = val
            elif partial.exists() and partial.stat().st_size >= minsz:
                cols = s.get("partial_cols", s["columns"])
                vals = s.get("partial_vals", s["complete"][:len(cols)])
                for col, val in zip(cols, vals):
                    df.at[idx, col] = val

    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    bak = input_file.with_suffix(f".bak.{stamp}.csv")
    df.to_csv(bak, sep=delim, index=False)
    logger.info(f"BACKUP → {bak}")
    df.to_csv(input_file, sep=delim, index=False)
    logger.info(f"REGEN → {input_file}")
    return df

def build_validated_set(input_file: Path, delim: str) -> List[str]:
    """Return RunAccession values where Validated == '1'."""
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
        logger.warning("No 'Validated' column yet — skipping")
        return []

def update_checkpoint(
    input_file: Path,
    lock_dir: Path,
    accession: str,
    field: str,
    value: str,
    delim: str
) -> None:
    """Update a single field in the CSV with a file lock."""
    lock_path = lock_dir / (input_file.name + ".lock")
    lock_path.touch(exist_ok=True)
    with lock_path.open("r+") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            df = pd.read_csv(input_file, sep=delim, dtype=str).fillna("")
            if field in df.columns and (df["RunAccession"] == accession).any():
                df.loc[df["RunAccession"] == accession, field] = value
                stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
                bak = input_file.with_suffix(f".bak.{stamp}.csv")
                df.to_csv(bak, sep=delim, index=False)
                df.to_csv(input_file, sep=delim, index=False)
                logger.debug(f"[CHECKPOINT] Updated {field}={value} for {accession}")
            else:
                logger.warning(f"[CHECKPOINT] Field '{field}' or accession '{accession}' not found")
        except Exception as e:
            logger.error(f"[CHECKPOINT ERROR] Failed to update {field} for {accession}: {e}")
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)

def run_command(cmd: str, description: str, env: Optional[Dict] = None, shell: bool = True) -> None:
    """Run a shell command and log output."""
    logger.info(description)
    try:
        result = subprocess.run(
            cmd, shell=shell, check=True, text=True, capture_output=True, env=env
        )
        logger.debug(f"✔ {cmd}: {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        logger.error(f"✘ {cmd}: {e.stderr.strip()}")
        raise
