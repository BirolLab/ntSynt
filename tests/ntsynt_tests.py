#!/usr/bin/env python3
'''
Tests for ntSynt
'''
import shlex
import subprocess

def launch_ntSynt(*genomes, prefix="ntSynt", k=24, w=1000, **kwargs):
    "Launch ntSynt"
    genome_list = " ".join(genomes)
    more_params = " ".join(f"--{name} {val}" for name, val in kwargs.items())
    cmd = f"ntSynt --force {genome_list} -k{k} -w {w} -d 0.5 --prefix {prefix} {more_params}"
    cmd = shlex.split(cmd)
    ret_code = subprocess.call(cmd)
    assert ret_code == 0

def launch_ntSynt_fof(genome_file_list, prefix="ntSynt", k=24, w=1000, **kwargs):
    "Launch ntSynt"
    more_params = " ".join(f"--{name} {val}" for name, val in kwargs.items())
    cmd = f"ntSynt --force --fastas_list {genome_file_list} -k{k} -w {w} -d 0.5 --prefix {prefix} {more_params}"
    cmd = shlex.split(cmd)
    ret_code = subprocess.call(cmd)
    assert ret_code == 0

def are_expected_blocks(block_file1, block_file2):
    "Check that the given synteny blocks are as expected"
    with open(block_file1, 'r', encoding="utf-8") as fin1:
        with open(block_file2, 'r', encoding="utf-8") as fin2:
            while (line1 := fin1.readline()) and (line2 := fin2.readline()):
                assert line1 == line2

def test_prep_files():
    "Prep the test files"
    infiles = "celegans-chrII-III.fa.gz", "celegans-chrII-III.A.fa.gz", "celegans-chrII-III.B.fa.gz"
    for infile in infiles:
        cmd = shlex.split(f"unpigz {infile}")
        ret_code = subprocess.call(cmd)
        assert ret_code == 0

def test_2_genomes():
    "Testing ntSynt with two input genomes"
    genome1, genome2 = "celegans-chrII-III.fa", "celegans-chrII-III.A.fa"
    launch_ntSynt(genome1, genome2, k=24, prefix="celegans-A-ntSynt", indel=500, merge=3000,
                  hashes=1)
    are_expected_blocks("celegans-A-ntSynt.synteny_blocks.tsv",
                        "expected_result/celegans-A-ntSynt.synteny_blocks.tsv")

def test_3_genomes():
    "Testing ntSynt with three input genomes"
    genome1, genome2, genome3 = "celegans-chrII-III.fa", "celegans-chrII-III.A.fa", "celegans-chrII-III.B.fa"
    launch_ntSynt(genome1, genome2, genome3, k=20, prefix="celegans-A-B-ntSynt", indel=500, merge=3000,
                  hashes=3)
    are_expected_blocks("celegans-A-B-ntSynt.synteny_blocks.tsv",
                        "expected_result/celegans-A-B-ntSynt.synteny_blocks.tsv")

def test_3_genomes_fof():
    "Testing ntSynt with three input genomes, specified using a file of files"
    genome_list = "genome_file_list.tsv"
    launch_ntSynt_fof(genome_list, k=20, prefix="celegans-A-B-ntSynt-fof", indel=500, merge=3000,
                      hashes=3)
    are_expected_blocks("celegans-A-B-ntSynt-fof.synteny_blocks.tsv",
                        "expected_result/celegans-A-B-ntSynt.synteny_blocks.tsv")

def test_cleanup():
    "Re-compress the files"
    infiles = "celegans-chrII-III.fa", "celegans-chrII-III.A.fa", "celegans-chrII-III.B.fa"
    for infile in infiles:
        cmd = shlex.split(f"pigz {infile}")
        ret_code = subprocess.call(cmd)
        assert ret_code == 0
