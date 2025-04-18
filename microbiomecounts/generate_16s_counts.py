#!/usr/bin/env python3

import argparse
import csv
import datetime
import gzip
import hashlib
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from threading import Lock
from typing import Dict, List, Optional, Tuple

# SLURM-like configuration (adapt as needed for your environment)
SLURM_CONFIG = {
    "job_name": "counts_16s",
    "output": "slurm_16s-%j.out",
    "error": "slurm_16s-%j.err",
    "time": "48:00:00",
    "nodes": 1,
    "ntasks": 1,
    "cpus_per_task": 16,
    "mem": "256G",
    "account": "one",
    "mail_user": "mak6930@psu.edu",
}

# Default variables and directories
DEFAULT_DIR = "fastq_data/fastq_biologicaldata"
OUTPUT_BASE = Path("dada2_16s")
LOCK_DIR = Path("locks")
FAILED_FILE = Path("failed_16s.log")
FAILED_LOCK = LOCK_DIR / "failed.lock"
INPUT_LOCK = LOCK_DIR / "input.lock"
DEFAULT_WORKERS = 4
RDP_DATABASE = "rdp_19_toGenus_trainset.fa.gz"

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# File‑locking for thread‑safe operations
failed_lock = Lock()
input_lock = Lock()

def setup_directories():
    """Create required directories and lock files."""
    for d in [LOCK_DIR, OUTPUT_BASE]:
        d.mkdir(parents=True, exist_ok=True)
    for f in [FAILED_LOCK, INPUT_LOCK]:
        f.touch()

def append_with_lock(line: str, file_path: Path, lock: Lock):
    """Append a line to a file with locking."""
    with lock:
        with file_path.open("a") as f:
            f.write(line + "\n")

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="16S-rRNA FASTQ processing pipeline")
    parser.add_argument("-i", "--input-file", required=True, help="Input CSV file")
    parser.add_argument("-d", "--fastq-dir", default=DEFAULT_DIR, help="FASTQ directory")
    parser.add_argument("-w", "--num-workers", type=int, default=DEFAULT_WORKERS, help="Number of workers")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("--no-repair-fastqs", action="store_false", dest="repair_fastqs", help="Disable FASTQ repair")
    parser.add_argument("--process-sample", action="store_true", help="Process a single sample")
    return parser.parse_args()

def setup_environment():
    """Check micromamba and active environment."""
    try:
        result = subprocess.run(["micromamba", "env", "list"], capture_output=True, text=True, check=True)
        envs = result.stdout.splitlines()
        active_env = next((line.split()[0] for line in envs if "*" in line), None)
        if not active_env:
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
    validated: List[str] = []
    with input_file.open("r") as f:
        first = f.readline().rstrip("\n\r")
    fields = first.split(delim)
    try:
        val_idx = fields.index("Validated")
        acc_idx = fields.index("RunAccession")
    except ValueError:
        return validated
    with input_file.open("r") as f:
        next(f)
        for line in f:
            row = line.rstrip("\n\r").split(delim)
            if len(row) > val_idx and row[val_idx].strip() == "1":
                validated.append(row[acc_idx].strip())
    return validated


def generate_fastq_checksums(fastq_dir: Path, num_workers: int):
    """Generate BLAKE3 checksums for FASTQ files."""
    logger.info("Generating initial checksums for FASTQ files")
    checksum_file = fastq_dir / "checksums.b3"
    fastq_files = sorted(fastq_dir.glob("*.fastq.gz"))
    def compute_checksum(f: Path) -> str:
        result = subprocess.run(["b3sum", str(f)], capture_output=True, text=True, check=True)
        return result.stdout.strip()
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        checksums = list(executor.map(compute_checksum, fastq_files))
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
        subprocess.run(["gzip", "-t", str(fq)], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        logger.error(f"[VAL FAIL] {fq} is not valid gzip")
        return False
    with gzip.open(fq, "rt") as f:
        lines = [f.readline().strip() for _ in range(3)]
    if not lines[0].startswith("@") or (len(lines) > 2 and not lines[2].startswith("+")):
        logger.error(f"[VAL FAIL] {fq} → Invalid FASTQ structure")
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
    return uniq_lengths > 1 or max_length > 30

def repair_fastq_if_needed(fastq_dir: Path, num_workers: int, validated_accessions: List[str]) -> Path:
    """Repair corrupted FASTQ files and return new FASTQ dir."""
    repaired_dir = fastq_dir.parent / f"{fastq_dir.name}_repaired"
    checksum_file = fastq_dir / "checksums.b3"
    failed_log = fastq_dir / "failed_checksums.txt"
    repaired_dir.mkdir(exist_ok=True)
    failed_log.write_text("")

    if not checksum_file.exists():
        generate_fastq_checksums(fastq_dir, num_workers)

    # Build list of files to check
    skip = {f"{acc}.fastq.gz" for acc in validated_accessions}
    to_check = []
    with checksum_file.open() as f:
        for line in f:
            hash_val, fname = line.strip().split(maxsplit=1)
            if Path(fname).name not in skip:
                to_check.append((hash_val, fname))

    def verify(h, fn):
        cmd = f"echo '{h}  {fn}' | b3sum -c - --quiet"
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            return None
        except subprocess.CalledProcessError:
            return fn

    # Collect failures without walrus
    failed_files = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(verify, h, f) for h, f in to_check]
        for future in as_completed(futures):
            res = future.result()
            if res:
                failed_files.append(res)

    # Write out failures
    failed_log.write_text("\n".join(failed_files) + "\n")
    logger.info(f"{len(failed_files)} files failed checksum")

    # Repair corrupted FASTQs
    if failed_files:
        logger.info("Repairing corrupted FASTQ files …")
        def do_repair(fn):
            inp = Path(fn)
            out = repaired_dir / inp.name
            with gzip.open(inp, "rt") as inf, gzip.open(out, "wt") as outf:
                for line in inf:
                    outf.write(line.rstrip("\r") + "\n")
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            executor.map(do_repair, failed_files)

    # Symlink the rest
    for fq in fastq_dir.glob("*.fastq.gz"):
        target = repaired_dir / fq.name
        if not target.exists():
            target.symlink_to(fq.resolve())

    # Regenerate checksums on repaired_dir
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



def update_input_with_fastq_paths(input_file: Path, fastq_dir: Path, delim: str, num_workers: int):
    """
    Update the input CSV with Fastq1, Fastq2, and Validated columns
    based on validation of FASTQ files in fastq_dir.

    Writes output back to 'input_file', logs summary stats and previews first few lines.
    """
    logger.info("Updating FASTQ paths")

    # 1) Build map: accession → list of FASTQ paths
    fastq_map: Dict[str, List[Path]] = {}
    for fq in fastq_dir.glob("*.fastq.gz"):
        acc = re.sub(r"(_[1-4])?\.fastq(\.gz)?$", "", fq.name)
        fastq_map.setdefault(acc, []).append(fq)

    total_acc = len(fastq_map)
    logger.info(f"Found {total_acc} unique accessions")

    # 2) Validate & classify in parallel
    validated_map: Dict[str, List[Path]] = {}
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(validate_and_check, fq, acc)
            for acc, fqs in fastq_map.items()
            for fq in fqs
        ]
        for future in as_completed(futures):
            acc, fq, reason = future.result()
            if reason:
                append_with_lock(f"{acc}:{reason}", FAILED_FILE, failed_lock)
            elif fq:
                validated_map.setdefault(acc, []).append(fq)

    validated_count = sum(1 for acc, vfs in validated_map.items() if vfs)
    missing_count = sum(1 for acc in fastq_map if not validated_map.get(acc))
    logger.info(f"Validated accessions: {validated_count}")
    logger.info(f"Accessions with no valid FASTQ: {missing_count}")

    # 3) Write directly to a tmp file next to the input file
    tmp_output = input_file.with_suffix(".tmp")
    logger.info(f"Writing updated input CSV to {tmp_output}")

    with tmp_output.open("w") as out:
        writer = csv.writer(out, delimiter=delim, lineterminator="\n")
        writer.writerow([
            "Bioproject", "RunAccession", "SequencingType",
            "Fastq1", "Fastq2", "Validated", "Dada2", "Completed"
        ])
        with input_file.open("r") as inp:
            reader = csv.reader(inp, delimiter=delim)
            next(reader)  # skip header
            for row in reader:
                biop, acc, st = row[:3]
                val  = row[5] if len(row) > 5 else "0"
                d2   = row[6] if len(row) > 6 else "0"
                comp = row[7] if len(row) > 7 else "0"

                if st != "16S":
                    f1 = row[3] if len(row) > 3 else ""
                    f2 = row[4] if len(row) > 4 else ""
                    writer.writerow([biop, acc, st, f1, f2, val, d2, comp])
                    continue

                vfs = validated_map.get(acc, [])
                if not vfs:
                    append_with_lock(f"{acc}:NO_VALID", FAILED_FILE, failed_lock)
                    new1, new2, new_val = "MISSING", "MISSING", "0"
                elif len(vfs) == 1:
                    new1, new2, new_val = str(vfs[0]), "", "1"
                else:
                    new1, new2, new_val = str(vfs[0]), str(vfs[1]), "1"

                writer.writerow([
                    biop, acc, st,
                    new1, new2, new_val,
                    d2, comp
                ])

    # 4) Atomically replace the original CSV
    shutil.move(str(tmp_output), str(input_file))
    logger.info(f"Replaced original input CSV: {input_file}")

    # 5) Preview first few lines of the updated CSV
    logger.info("Preview of updated CSV:")
    with input_file.open() as f:
        for i, line in enumerate(f):
            if i < 6:
                logger.info("  " + line.strip())
            else:
                break


def update_checkpoint(input_file: Path, accession: str, field: str, value: str, delim: str):
    """Update a checkpoint field in the input CSV."""
    idx_map = {"Validated":5,"Dada2":6,"Completed":7}
    idx = idx_map.get(field)
    if idx is None:
        return
    tmp = input_file.with_suffix(".tmp")
    with input_file.open() as inp, tmp.open("w") as out:
        reader = csv.reader(inp, delimiter=delim)
        writer = csv.writer(out, delimiter=delim, lineterminator="\n")
        header = next(reader)
        writer.writerow(header)
        for row in reader:
            if row and row[1]==accession:
                row[idx]=value
            writer.writerow(row)
    shutil.move(str(tmp), str(input_file))

def process_sample(biop: str, acc: str, st: str, fq1: str, fq2: str, val: str, d2: str, comp: str, env_name: str):
    """Process a single sample."""
    if st!="16S" or val=="0":
        return
    if not (fq1 or fq2):
        append_with_lock(f"{acc}:NO_FASTQ", FAILED_FILE, failed_lock)
        return
    if comp=="1":
        return

    outdir = OUTPUT_BASE/biop
    outdir.mkdir(parents=True, exist_ok=True)
    seqtab = outdir/f"asv_{acc}.rds"

    if d2!="1":
        cmd = f"Rscript {Path.cwd()}/run_dada2_partial.R {fq1} {fq2} {seqtab}" if fq2 else f"Rscript {Path.cwd()}/run_dada2_partial.R {fq1} {seqtab}"
        try:
            run_command(cmd, f"DADA2 on {acc}", env_name)
            if not seqtab.exists():
                append_with_lock(f"{acc}:DADA2_FAIL", FAILED_FILE, failed_lock)
                return
            update_checkpoint(Path(os.environ["INPUT_FILE"]), acc, "Dada2", "1", os.environ["DELIM"])
        except subprocess.CalledProcessError:
            append_with_lock(f"{acc}:DADA2_FAIL", FAILED_FILE, failed_lock)
            return

    update_checkpoint(Path(os.environ["INPUT_FILE"]), acc, "Completed", "1", os.environ["DELIM"])
    logger.info(f"Finished {acc}")

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
        samples = [row for row in reader if row[2]=="16S"]
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_sample, *sample, env_name) for sample in samples]
        for f in as_completed(futures):
            f.result()  # will raise if any failed
    logger.info("All samples finished")

def main():
    args = parse_arguments()
    if args.debug:
        logger.setLevel(logging.DEBUG)

    setup_directories()
    env_name = setup_environment()

    # Read and export CSV delimiter & path
    delim = validate_input_file(Path(args.input_file))
    os.environ["DELIM"] = delim
    os.environ["INPUT_FILE"] = args.input_file
    os.environ["FAILED_FILE"] = str(FAILED_FILE)
    os.environ["OUTPUT_BASE"] = str(OUTPUT_BASE)

    check_rdp_database()
    fastq_dir = Path(args.fastq_dir).resolve()

    # === PRE‑DADA2 STEPS (use ALL CPUs) ===
    pre_workers = os.cpu_count() or DEFAULT_WORKERS

    # 1) Generate checksums if missing
    if not (fastq_dir / "checksums.b3").exists():
        generate_fastq_checksums(fastq_dir, pre_workers)

    # 2) Gather already-validated accessions
    validated_accessions = build_validated_set(Path(args.input_file), delim)

    # 3) Repair FASTQs in parallel (if enabled)
    if args.repair_fastqs:
        fastq_dir = repair_fastq_if_needed(fastq_dir, pre_workers, validated_accessions)

    # 4) Update input.csv with Fastq1/Fastq2/Validated using the same pool size
    update_input_with_fastq_paths(
        Path(args.input_file),
        fastq_dir,
        delim,
        pre_workers
    )

    # === DADA2 processing (throttled to --num-workers) ===
    process_samples(Path(args.input_file), delim, args.num_workers, env_name)
    final_validation_and_merge(Path(args.input_file), delim, env_name)

    logger.info("Pipeline complete.")

if __name__ == "__main__":
    if "--process-sample" in sys.argv:
        logger.error("process-sample mode not fully supported in Python version")
        sys.exit(1)
    main()
