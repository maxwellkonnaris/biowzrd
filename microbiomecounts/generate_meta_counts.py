#!/usr/bin/env python3
# counts_16s.py

import argparse
import csv
import logging
import os
import sys
import time
import subprocess
from datetime import datetime
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

from utils_counts import (
    setup_logging,
    compute_worker_counts,
    submit_slurm_job,
    setup_directories,
    validate_input_file,
    check_database,
    generate_fastq_checksums,
    build_validated_set,
    repair_fastq_if_needed,
    update_input_with_fastq_paths,
    regenerate_flags,
    ensure_required_columns,
    append_with_flock,
    update_checkpoint,
    run_command,
    logger,
    pd
)

# SLURM configuration
SLURM_CONFIG = {
    "job_name": "counts_16s",
    "output": "slurm_16s-%j.out",
    "error": "slurm_16s-%j.err",
    "time": "48:00:00",
    "nodes": 1,
    "ntasks": 1,
    "cpus_per_task": 32,
    "mem": "512G",
    "account": "one",
    "mail_user": "mak6930@psu.edu",
}

# Default paths and parameters
DEFAULT_FASTQ_DIR = "fastq_data/fastq_biologicaldata"
OUTPUT_BASE = Path("dada2_16s")
LOCK_DIR = Path("locks")
FAILED_FILE = Path("failed_16s.log")
RDP_DATABASE = Path("rdp_19_toGenus_trainset.fa.gz")
CONDA_ENV = "dada2_env"
REQUIRED_COLUMNS = [
    "Bioproject", "RunAccession", "SequencingType",
    "Fastq1", "Fastq2", "Validated", "Dada2", "Completed"
]
ZERO_DEFAULTS = {"Validated", "Dada2", "Completed"}
OUTPUT_COLS = ["Fastq1", "Fastq2", "Validated"]
SEQ_TYPE = "16S"

# Regeneration specs for DADA2
DADA2_SPECS = {
    "dada2": {
        "merged": "{biop}_merged.rds",
        "partial": ("asv_{acc}.rds", 100),
        "columns": ["Validated", "Dada2", "Completed"],
        "complete": ["1", "1", "1"],
        "partial_cols": ["Validated", "Dada2"],
        "partial_vals": ["1", "1"],
    }
}

def parse_args():
    parser = argparse.ArgumentParser(description="16S-rRNA FASTQ processing pipeline")
    parser.add_argument("-i", "--input-file", required=True, help="Path to your input CSV file")
    parser.add_argument("-d", "--fastq-dir", default=DEFAULT_FASTQ_DIR,
                        help="Directory containing your FASTQ files")
    parser.add_argument("--env", default="dada2",
                        help="Micromamba/conda environment name to activate")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("--no-repair-fastqs", dest="repair_fastqs", action="store_false",
                        help="Skip the FASTQ repair step if desired")
    parser.add_argument("--submit-slurm", action="store_true",
                        help="Submit this script as a SLURM job and exit")
    parser.add_argument("--process-sample", action="store_true", help="Process a single sample")
    return parser.parse_args()

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
        logger.info(f"Using activated DADA2 environment: {active_env}")
        return active_env
    except subprocess.CalledProcessError:
        logger.error("micromamba not in PATH")
        sys.exit(1)

def process_sample(biop: str, acc: str, st: str, fq1: str, fq2: str, val: str, d2: str, comp: str, env_name: str):
    """Process a single sample with DADA2."""
    if st != SEQ_TYPE or val == "0":
        return
    if not (fq1 or fq2):
        append_with_flock(f"{acc}:NO_FASTQ", FAILED_FILE)
        return
    if comp == "1":
        return

    outdir = OUTPUT_BASE / biop
    outdir.mkdir(parents=True, exist_ok=True)
    seqtab = outdir / f"asv_{acc}.rds"

    if d2 != "1":
        cmd = f"Rscript {Path.cwd()}/run_dada2_partial.R {fq1} {fq2} {seqtab}" if fq2 else f"Rscript {Path.cwd()}/run_dada2_partial.R {fq1} {seqtab}"
        try:
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = "4"
            env["MKL_NUM_THREADS"] = "4"
            run_command(cmd, f"Running DADA2 for {acc}", env=env)
            if not seqtab.exists():
                append_with_flock(f"{acc}:DADA2_FAIL", FAILED_FILE)
                return
            update_checkpoint(Path(args.input_file), LOCK_DIR, acc, "Dada2", "1", delim)
        except subprocess.CalledProcessError:
            append_with_flock(f"{acc}:DADA2_FAIL", FAILED_FILE)
            return

    update_checkpoint(Path(args.input_file), LOCK_DIR, acc, "Completed", "1", delim)
    logger.info(f"Finished {acc}")

def merge_profiles(biop: str, env_name: str):
    """Merge DADA2 profiles for a bioproject."""
    odir = OUTPUT_BASE / biop
    with Path(args.input_file).open() as f:
        reader = csv.reader(f, delimiter=delim)
        next(reader)
        for row in reader:
            if row[0] == biop and row[2] == SEQ_TYPE and row[-1] != "1":
                return
    tabs = list(odir.glob("asv_*.rds"))
    if tabs:
        cmd = f"Rscript {Path.cwd()}/merge_dada2.R {biop} {' '.join(map(str, tabs))}"
        run_command(cmd, f"Merging {biop}", env_name)

def process_samples(input_file: Path, delim: str, num_workers: int, env_name: str):
    """Process samples in parallel."""
    with input_file.open() as f:
        reader = csv.reader(f, delimiter=delim)
        next(reader)
        samples = [row for row in reader if row[2] == SEQ_TYPE]
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_sample, *sample, env_name) for sample in samples]
        for f in as_completed(futures):
            f.result()
    logger.info("All samples finished")

def final_validation_and_merge(input_file: Path, delim: str, env_name: str):
    """Run final validation and merge profiles."""
    bioprojects = set()
    with input_file.open() as f:
        reader = csv.reader(f, delimiter=delim)
        next(reader)
        for row in reader:
            if row[2] == SEQ_TYPE:
                bioprojects.add(row[0])
    for bp in sorted(bioprojects):
        merge_profiles(bp, env_name)
    logger.info("Final validation done")

def regenerate_input_with_validation(input_file: Path, output_dir: Path, delim: str = ",") -> pd.DataFrame:
    """Regenerate flags for 16S pipeline using Rscript validation."""
    df = pd.read_csv(input_file, sep=delim, dtype=str)
    df["Validated"] = "0"
    df["Dada2"] = "0"
    df["Completed"] = "0"

    for idx, row in df.iterrows():
        biop = row["Bioproject"]
        acc = row["RunAccession"]
        biop_dir = output_dir / biop
        asv_file = biop_dir / f"asv_{acc}.rds"
        merged_file = biop_dir / f"{biop}_merged.rds"

        if merged_file.exists():
            df.at[idx, "Validated"] = "1"
            df.at[idx, "Dada2"] = "1"
            df.at[idx, "Completed"] = "1"
            continue

        if not asv_file.exists() or asv_file.stat().st_size < 100:
            continue

        try:
            rscript = f"""
            args <- commandArgs(trailingOnly=TRUE)
            x <- readRDS(args[1])
            if ((is.matrix(x) || is.data.frame(x)) && ncol(x) > 0) {{
              quit(status=0)
            }} else {{
              quit(status=1)
            }}
            """
            subprocess.run(
                ["Rscript", "-e", rscript, str(asv_file)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True
            )
            df.at[idx, "Validated"] = "1"
            df.at[idx, "Dada2"] = "1"
        except subprocess.CalledProcessError:
            pass

    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    backup_path = input_file.with_suffix(f".bak.{timestamp}.csv")
    df.to_csv(backup_path, sep=delim, index=False)
    logger.info(f"[BACKUP] Saved backup to {backup_path}")

    df.to_csv(input_file, sep=delim, index=False)
    logger.info(f"[REGEN] Updated input file based on existing outputs → {input_file}")
    return df

def main():
    global args, delim
    args = parse_args()
    setup_logging(logging.DEBUG if args.debug else logging.INFO)

    if args.process_sample:
        logger.error("process-sample mode not fully supported in Python version")
        sys.exit(1)

    if "SLURM_JOB_ID" not in os.environ and not args.submit_slurm:
        logger.info("Not in a SLURM job; submitting one...")
        submit_slurm_job(SLURM_CONFIG, args.env, args)
        sys.exit(0)

    pre_workers = compute_worker_counts(processing_step=False, seq_type=SEQ_TYPE)
    dada2_workers = compute_worker_counts(processing_step=True, seq_type=SEQ_TYPE)
    logger.info(f"Using {pre_workers} workers for I/O and validation steps ({pre_workers} CPUs)")
    logger.info(f"Using {dada2_workers} workers for DADA2 processing ({dada2_workers * 4} CPUs)")

    setup_directories([Path(args.fastq_dir), LOCK_DIR, OUTPUT_BASE])
    env_name = setup_environment(args.env)
    delim = validate_input_file(Path(args.input_file))

    check_database(RDP_DATABASE)
    fastq_dir = Path(args.fastq_dir).resolve()

    start = time.time()
    if not (fastq_dir / "checksums.b3").exists():
        generate_fastq_checksums(fastq_dir, pre_workers)
    logger.info(f"Checksum generation took {time.time() - start:.2f} seconds")

    start = time.time()
    validated = build_validated_set(Path(args.input_file), delim)
    logger.info(f"Validated set building took {time.time() - start:.2f} seconds")

    start = time.time()
    if args.repair_fastqs:
        fastq_dir = repair_fastq_if_needed(fastq_dir, pre_workers, validated)
    logger.info(f"FASTQ repair took {time.time() - start:.2f} seconds")

    start = time.time()
    update_input_with_fastq_paths(
        Path(args.input_file),
        fastq_dir,
        delim,
        pre_workers,
        seq_type=SEQ_TYPE,
        output_cols=OUTPUT_COLS
    )
    logger.info(f"FASTQ path update took {time.time() - start:.2f} seconds")

    start = time.time()
    if OUTPUT_BASE.exists() and any(OUTPUT_BASE.iterdir()):
        df = regenerate_input_with_validation(Path(args.input_file), OUTPUT_BASE, delim)
        df.to_csv(Path(args.input_file), sep=delim, index=False)
        logger.info("Checkpoint fields restored from existing DADA2 output")
    logger.info(f"Flag regeneration took {time.time() - start:.2f} seconds")

    start = time.time()
    process_samples(Path(args.input_file), delim, dada2_workers, env_name)
    logger.info(f"DADA2 processing took {time.time() - start:.2f} seconds")

    start = time.time()
    final_validation_and_merge(Path(args.input_file), delim, env_name)
    logger.info(f"Final validation took {time.time() - start:.2f} seconds")

    logger.info("Pipeline complete.")

if __name__ == "__main__":
    main()
