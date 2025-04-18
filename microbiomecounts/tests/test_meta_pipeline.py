import gzip
import pytest
from pathlib import Path

# adjust this import to your actual meta‐script name
import generate_meta_counts as pm

def make_csv(tmp_path, header, rows, delim=","):
    f = tmp_path / "input_meta.csv"
    with f.open("w") as fh:
        fh.write(header + "\n")
        for row in rows:
            fh.write(delim.join(row) + "\n")
    return f

def make_fastq(tmp_path, name, sequences, qualities):
    fq = tmp_path / name
    with gzip.open(fq, "wt") as fh:
        for seq, qual in zip(sequences, qualities):
            fh.write("@id\n")
            fh.write(seq + "\n")
            fh.write("+\n")
            fh.write(qual + "\n")
    return fq

def test_validate_input_file(tmp_path):
    f = make_csv(tmp_path, "Bioproject,RunAccession,SequencingType", [])
    assert pm.validate_input_file(f) == ","

def test_build_validated_set(tmp_path):
    hdr = "Bioproject,RunAccession,SequencingType,Fastq1,Fastq2,Validated"
    rows = [
        ["BPX","RUN1","meta","","","1"],
        ["BPX","RUN2","meta","","","0"],
        ["BPY","RUN3","meta","","","1"]
    ]
    f = make_csv(tmp_path, hdr, rows)
    vs = pm.build_validated_set(f, ",")
    assert sorted(vs) == ["RUN1","RUN3"]

def test_validate_fastq_and_is_biological(tmp_path):
    seqs = ["TTTT","TTTT"]
    quals = ["IIII","IIII"]
    fq1 = make_fastq(tmp_path, "m1.fastq.gz", seqs, quals)
    assert pm.validate_fastq(fq1)
    assert not pm.is_biological_fastq(fq1)
    seqs2 = ["AAAAAAA","AAAAAA"]  # differing lengths
    quals2 = ["IIIIIII","IIIIII"]
    fq2 = make_fastq(tmp_path, "m2.fastq.gz", seqs2, quals2)
    assert pm.validate_fastq(fq2)
    assert pm.is_biological_fastq(fq2)

def test_validate_fastq_error(tmp_path):
    bad = tmp_path / "bad_meta.fastq.gz"
    bad.write_text("oops")
    assert not pm.validate_fastq(bad)
