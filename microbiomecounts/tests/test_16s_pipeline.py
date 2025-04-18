import gzip
import os
import csv
import pytest
from pathlib import Path

# adjust this import to your actual script name
import generate_16s_counts as p16

def make_csv(tmp_path, header, rows, delim=","):
    f = tmp_path / "input.csv"
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

def test_validate_input_file_comma(tmp_path):
    f = make_csv(tmp_path, "Bioproject,RunAccession,SequencingType", [])
    assert p16.validate_input_file(f) == ","

def test_validate_input_file_tab(tmp_path):
    f = make_csv(tmp_path, "Bioproject\tRunAccession\tSequencingType", [], delim="\t")
    assert p16.validate_input_file(f) == "\t"

def test_validate_input_file_bad(tmp_path):
    f = make_csv(tmp_path, "Foo,Bar,Baz", [])
    with pytest.raises(SystemExit):
        p16.validate_input_file(f)

def test_build_validated_set(tmp_path):
    hdr = "Bioproject,RunAccession,SequencingType,Fastq1,Fastq2,Validated"
    rows = [
        ["BP1","ACC1","16S","","","1"],
        ["BP2","ACC2","16S","","","0"],
        ["BP3","ACC3","16S","","","1"]
    ]
    f = make_csv(tmp_path, hdr, rows)
    vs = p16.build_validated_set(f, ",")
    assert sorted(vs) == ["ACC1","ACC3"]

def test_validate_fastq_and_is_biological(tmp_path):
    # create one “technical” (identical short seqs) FASTQ
    seqs = ["A"*10, "A"*10]
    quals = ["I"*10, "I"*10]
    fq1 = make_fastq(tmp_path, "tech.fastq.gz", seqs, quals)
    assert p16.validate_fastq(fq1)
    assert not p16.is_biological_fastq(fq1)

    # and one “biological” (varying lengths)
    seqs2 = ["A"*10, "A"*12]
    quals2 = ["I"*10, "I"*12]
    fq2 = make_fastq(tmp_path, "bio.fastq.gz", seqs2, quals2)
    assert p16.validate_fastq(fq2)
    assert p16.is_biological_fastq(fq2)

def test_validate_fastq_gzip_error(tmp_path):
    bad = tmp_path / "bad.fastq.gz"
    bad.write_text("not a gzipped file")
    assert not p16.validate_fastq(bad)
