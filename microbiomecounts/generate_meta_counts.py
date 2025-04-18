#!/usr/bin/env python3

import argparse
import csv
import gzip
import logging
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from threading import Lock
from typing import Dict, List, Optional, Tuple

# SLURM-like configuration
SLURM_CONFIG = {
    "job_name": "counts_meta",
    "output": "slurm_meta-%j.out",
    "error": "slurm_meta-%j.err",
    "time": "48:00:00",
    "nodes": 1,
    "ntasks": 1,
    "cpus_per_task": 32,
    "mem": "512G",
    "account": "open",
    "mail_user": "mak6930@psu.edu",
}

# Default variables and directories
DEFAULT_DIR = "fastq_data/fastq_biologicaldata"
OUTPUT_BASE = Path("metagenome")
LOCK_DIR = Path("locks_meta")
FAILED_FILE = Path("failed_meta.log")
FAILED_LOCK = LOCK_DIR / "failed.lock"
INPUT_LOCK = LOCK_DIR / "input.lock"
DEFAULT_WORKERS = 8
DEFAULT_MOTUS_TAX_LEVEL = "mOTU"
METAPHLAN_DB_FALLBACK = "/storage/work/mak6930/applicationstorage/micromamba/envs/mpa/lib/python3.7/site-packages/metaphlan/metaphlan_databases"

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# File-locking for thread-safe operations
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
    parser = argparse.ArgumentParser(description="Metagenomic FASTQ processing pipeline")
    parser.add_argument("-i", "--input-file", required=True, help="Input CSV file")
    parser.add_argument("-d", "--fastq-dir", default=DEFAULT_DIR, help="FASTQ directory")
    parser.add_argument("-k", "--motus-tax-level", default=DEFAULT_MOTUS_TAX_LEVEL, help="mOTUs taxonomy level")
    parser.add_argument("-w", "--num-workers", type=int, default=DEFAULT_WORKERS, help="Number of workers")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("--no-repair-fastqs", action="store_false", dest="repair_fastqs", help="Disable FASTQ repair")
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
        logger.info(f"Using activated MPA environment: {active_env}")
        return active_env
    except subprocess.CalledProcessError:
        logger.error("micromamba not in PATH")
        sys.exit(1)

def get_conda_prefix(env_name: str) -> Path:
    """Get the micromamba environment prefix."""
    try:
        result = subprocess.run(
            f"micromamba env list | grep {env_name}",
            shell=True, capture_output=True, text=True, check=True
        )
        prefix = Path(result.stdout.split()[-1])
        return prefix
    except subprocess.CalledProcessError:
        logger.error(f"Could not determine prefix for environment {env_name}")
        sys.exit(1)

def run_command(cmd: str, description: str, env_name: str = None) -> None:
    logger.info(description)
    subprocess.run(cmd, shell=True, check=True, text=True)


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

def check_metaphlan_database(env_name: str) -> Path:
    """Detect or install MetaPhlAn database with fallback."""
    try:
        result = subprocess.run(
            f"micromamba run -n {env_name} python -c 'import metaphlan, pathlib; print(pathlib.Path(metaphlan.__file__).parent / \"metaphlan_databases\")'",
            shell=True, capture_output=True, text=True, check=True
        )
        db_path = Path(result.stdout.strip())
    except subprocess.CalledProcessError:
        logger.debug(f"MetaPhlAn DB auto-detection failed; using fallback: {METAPHLAN_DB_FALLBACK}")
        db_path = Path(METAPHLAN_DB_FALLBACK)

    if not db_path.is_dir():
        logger.error(f"MetaPhlAn DB not found at {db_path}")
        sys.exit(1)

    if not list(db_path.glob("*.bt2l")):
        logger.info("MetaPhlAn DB empty – installing...")
        run_command(f"metaphlan --install --bowtie2db \"{db_path}\"", "Installing MetaPhlAn DB", env_name)
    logger.info(f"MetaPhlAn DB ok: {db_path}")
    return db_path

def build_validated_set(input_file: Path, delim: str) -> List[str]:
    """Read the input CSV and return all RunAccession values where Validated == '1'."""
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
    """Validate a FASTQ file."""
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

    failed_files = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(verify, h, f) for h, f in to_check]
        for future in as_completed(futures):
            res = future.result()
            if res:
                failed_files.append(res)

    failed_log.write_text("\n".join(failed_files) + "\n")
    logger.info(f"{len(failed_files)} files failed checksum")

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

    for fq in fastq_dir.glob("*.fastq.gz"):
        target = repaired_dir / fq.name
        if not target.exists():
            target.symlink_to(fq.resolve())

    generate_fastq_checksums(repaired_dir, num_workers)
    logger.info(f"Repair complete → FASTQ_DIR updated to {repaired_dir}")
    return repaired_dir

def validate_and_check(fq: Path, acc: str) -> Tuple[str, Optional[Path], Optional[str]]:
    """Validate a FASTQ file and determine biological validity."""
    if validate_fastq(fq):
        if is_biological_fastq(fq):
            return acc, fq, None
        else:
            return acc, None, "TECHNICAL"
    else:
        return acc, None, "CORRUPT"

def update_input_with_fastq_paths(input_file: Path, fastq_dir: Path, delim: str, num_workers: int):
    """Update the input CSV with Fastq1, Fastq2, and Validated columns."""
    logger.info("Updating FASTQ paths")

    fastq_map: Dict[str, List[Path]] = {}
    for fq in fastq_dir.glob("*.fastq.gz"):
        acc = re.sub(r"(_[1-4])?\.fastq(\.gz)?$", "", fq.name)
        fastq_map.setdefault(acc, []).append(fq)

    total_acc = len(fastq_map)
    logger.info(f"Found {total_acc} unique accessions")

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

    tmp_output = input_file.with_suffix(".tmp")
    logger.info(f"Writing updated input CSV to {tmp_output}")

    with tmp_output.open("w") as out:
        writer = csv.writer(out, delimiter=delim, lineterminator="\n")
        writer.writerow([
            "Bioproject", "RunAccession", "SequencingType",
            "Fastq1", "Fastq2", "Validated", "Motus", "Metaphlan", "Completed"
        ])
        with input_file.open("r") as inp:
            reader = csv.reader(inp, delimiter=delim)
            next(reader)
            for row in reader:
                biop, acc, st = row[:3]
                val = row[5] if len(row) > 5 else "0"
                mt = row[6] if len(row) > 6 else "0"
                mp = row[7] if len(row) > 7 else "0"
                comp = row[8] if len(row) > 8 else "0"

                if st != "meta":
                    f1 = row[3] if len(row) > 3 else ""
                    f2 = row[4] if len(row) > 4 else ""
                    writer.writerow([biop, acc, st, f1, f2, val, mt, mp, comp])
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
                    mt, mp, comp
                ])

    shutil.move(str(tmp_output), str(input_file))
    logger.info(f"Replaced original input CSV: {input_file}")

    logger.info("Preview of updated CSV:")
    with input_file.open() as f:
        for i, line in enumerate(f):
            if i < 6:
                logger.info("  " + line.strip())
            else:
                break

def update_checkpoint(input_file: Path, accession: str, field: str, value: str, delim: str):
    """Update a checkpoint field in the input CSV."""
    idx_map = {"Validated": 5, "Motus": 6, "Metaphlan": 7, "Completed": 8}
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
            if row and row[1] == accession:
                row[idx] = value
            writer.writerow(row)
    shutil.move(str(tmp), str(input_file))

def convert_metaphlan_to_counts(log_file: Path, profile: Path, output: Path) -> bool:
    """Convert MetaPhlAn profile to counts."""
    if not profile.is_file():
        return False
    with log_file.open() as f:
        for line in f:
            if "Total number of reads mapped" in line:
                mapped = line.split()[5].split("(")[0]
                break
        else:
            return False
    with profile.open() as inf, output.open("w") as outf:
        outf.write("#clade\tabundance\treads\n")
        next(inf)
        for line in inf:
            clade, abundance = line.strip().split("\t")[:2]
            reads = float(mapped) * float(abundance) / 100
            outf.write(f"{clade}\t{abundance}\t{reads}\n")
    return True

def process_sample(biop: str, acc: str, st: str, fq1: str, fq2: str, val: str, mt: str, mp: str, comp: str, env_name: str, metaphlan_db: Path, motus_tax_level: str):
    """Process a single metagenomic sample."""
    if st != "meta" or val == "0":
        return
    if not (fq1 or fq2):
        append_with_lock(f"{acc}:NO_FASTQ", FAILED_FILE, failed_lock)
        return
    if comp == "1":
        return

    outdir = OUTPUT_BASE / biop
    outdir.mkdir(parents=True, exist_ok=True)
    mpl_log = outdir / f"{acc}_metaphlan_log.txt"
    mpl_bt2 = outdir / f"{acc}_bt2.txt"
    mpl_out = outdir / f"{acc}_metaphlan4.txt"
    mpl_counts = outdir / f"{acc}_metaphlan4_counts.txt"
    motus_out = outdir / f"{acc}_motus.txt"

    if mp != "1":
        infq = f"{fq1},{fq2}" if fq2 else fq1
        cmd = f"metaphlan \"{infq}\" --input_type fastq --unclassified_estimation --nproc {SLURM_CONFIG['cpus_per_task']} --bowtie2db \"{metaphlan_db}\" --bowtie2out \"{mpl_bt2}\" -o \"{mpl_out}\" 2> \"{mpl_log}\""
        try:
            run_command(cmd, f"MetaPhlAn on {acc}", env_name)
            if not convert_metaphlan_to_counts(mpl_log, mpl_out, mpl_counts):
                append_with_lock(f"{acc}:MPL_convert_fail", FAILED_FILE, failed_lock)
                return
            update_checkpoint(Path(os.environ["INPUT_FILE"]), acc, "Metaphlan", "1", os.environ["DELIM"])
        except subprocess.CalledProcessError:
            append_with_lock(f"{acc}:MPL_fail", FAILED_FILE, failed_lock)
            return

    if mt != "1":
        if fq2:
            cmd = f"motus profile -f \"{fq1}\" -r \"{fq2}\" -o \"{motus_out}\" -t {SLURM_CONFIG['cpus_per_task']} -c -k {motus_tax_level}"
        else:
            cmd = f"motus profile -s \"{fq1}\" -o \"{motus_out}\" -t {SLURM_CONFIG['cpus_per_task']} -c -k {motus_tax_level}"
        try:
            run_command(cmd, f"mOTUs {'paired' if fq2 else 'single-end'} on {acc}", env_name)
            if not motus_out.exists():
                append_with_lock(f"{acc}:mOTUs_fail", FAILED_FILE, failed_lock)
                return
            update_checkpoint(Path(os.environ["INPUT_FILE"]), acc, "Motus", "1", os.environ["DELIM"])
        except subprocess.CalledProcessError:
            append_with_lock(f"{acc}:mOTUs_fail", FAILED_FILE, failed_lock)
            return

    update_checkpoint(Path(os.environ["INPUT_FILE"]), acc, "Completed", "1", os.environ["DELIM"])
    logger.info(f"Finished {acc}")

def merge_profiles(biop: str, tool: str, env_name: str, conda_prefix: Path):
    """Merge MetaPhlAn or mOTUs profiles for a bioproject."""
    odir = OUTPUT_BASE / biop
    if tool == "metaphlan":
        files = list(odir.glob("*_metaphlan4_counts.txt"))
        if not files:
            return
        out = odir / f"{biop}_metaphlan_merged.txt"
        merge_script = conda_prefix / "lib" / "python3.7" / "site-packages" / "metaphlan" / "utils" / "merge_metaphlan_tables.py"
        cmd = f"python \"{merge_script}\" {' '.join(map(str, files))} > \"{out}\""
    else:  # motus
        files = list(odir.glob("*_motus.txt"))
        if not files:
            return
        out = odir / f"{biop}_motus_merged.txt"
        file_list = ",".join(map(str, files))
        cmd = f"motus merge -i \"{file_list}\" -o \"{out}\""
    run_command(cmd, f"Merging {biop} {tool}", env_name)

def final_validation_and_merge(input_file: Path, delim: str, env_name: str, conda_prefix: Path):
    """Run final validation and merge profiles."""
    bioprojects = set()
    with input_file.open() as f:
        reader = csv.reader(f, delimiter=delim)
        next(reader)
        for row in reader:
            if row[2] == "meta":
                bioprojects.add(row[0])
    for bp in sorted(bioprojects):
        merge_profiles(bp, "metaphlan", env_name, conda_prefix)
        merge_profiles(bp, "motus", env_name, conda_prefix)
    logger.info("Final merge done")

def process_samples(input_file: Path, delim: str, num_workers: int, env_name: str, metaphlan_db: Path, motus_tax_level: str, conda_prefix: Path):
    """Process samples in parallel."""
    with input_file.open() as f:
        reader = csv.reader(f, delimiter=delim)
        next(reader)
        samples = [row for row in reader if row[2] == "meta"]
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(process_sample, *sample, env_name, metaphlan_db, motus_tax_level)
            for sample in samples
        ]
        for f in as_completed(futures):
            f.result()
    logger.info("All samples finished")

def main():
    args = parse_arguments()
    if args.debug:
        logger.setLevel(logging.DEBUG)

    setup_directories()
    env_name = setup_environment()
    conda_prefix = get_conda_prefix(env_name)

    delim = validate_input_file(Path(args.input_file))
    os.environ["DELIM"] = delim
    os.environ["INPUT_FILE"] = args.input_file
    os.environ["FAILED_FILE"] = str(FAILED_FILE)
    os.environ["OUTPUT_BASE"] = str(OUTPUT_BASE)

    metaphlan_db = check_metaphlan_database(env_name)
    fastq_dir = Path(args.fastq_dir).resolve()

    pre_workers = os.cpu_count() or DEFAULT_WORKERS

    if not (fastq_dir / "checksums.b3").exists():
        generate_fastq_checksums(fastq_dir, pre_workers)

    validated_accessions = build_validated_set(Path(args.input_file), delim)

    if args.repair_fastqs:
        fastq_dir = repair_fastq_if_needed(fastq_dir, pre_workers, validated_accessions)

    update_input_with_fastq_paths(
        Path(args.input_file),
        fastq_dir,
        delim,
        pre_workers
    )

    process_samples(
        Path(args.input_file),
        delim,
        args.num_workers,
        env_name,
        metaphlan_db,
        args.motus_tax_level,
        conda_prefix
    )
    final_validation_and_merge(Path(args.input_file), delim, env_name, conda_prefix)

    logger.info("Pipeline complete.")

if __name__ == "__main__":
    main()
