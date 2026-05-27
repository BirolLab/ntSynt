#!/usr/bin/bash

#Run ntSynt demo to test installation

echo "Running ntSynt Demo..."
echo "This demo assumes that you have the 'ntSynt' script in your PATH"

set -eux -o pipefail

echo "Decompressing input files.."
if [ -e celegans-chrII-III.fa.gz ]; then
	unpigz celegans-chrII-III.fa.gz celegans-chrII-III.A.fa.gz celegans-chrII-III.B.fa.gz
fi

echo "Running ntSynt with 2 input genomes"
ntSynt celegans-chrII-III.fa celegans-chrII-III.A.fa --prefix celegans-A-ntSynt -d 0.5 --merge 3000 --indel 500 --hashes 1 --force

echo "Running ntSynt with 3 input genomes"
ntSynt celegans-chrII-III.fa celegans-chrII-III.A.fa celegans-chrII-III.B.fa -k20 --prefix celegans-A-B-ntSynt -d 0.5 --merge 3000 --indel 500 --force

echo "Compressing input files.."
if [ -e celegans-chrII-III.fa ]; then
	pigz celegans-chrII-III.fa celegans-chrII-III.A.fa celegans-chrII-III.B.fa
fi

echo "DONE! Compare your output files with those in 'expected_result' to ensure that your installation is generating the expected results"
