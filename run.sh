#!/bin/bash

# Download data
# DATA_DIR="data/"
source config.ini


usage()
{
cat << EOF
usage: $0 options


OPTIONS:
   -h      Show this message
   -f      download relevant files
   -p      preprocess
   -c      compress
   -d      Decompress       
 		
EOF
}


# TODO: Make this more general later
download()
{
	echo "*** FASTQ Sequences being downloaded ***"
	mkdir -p data/$basename
	mkdir -p data/$basename/output
	wget -O data/$basename/file_1.fastq.gz $URL_1
	wget -O data/$basename/file_2.fastq.gz $URL_2

	gunzip data/$basename/file_1.fastq.gz data/$basename/file_2.fastq.gz
	cat data/$basename/file_1.fastq data/$basename/file_2.fastq > data/$basename/input.fastq
}

# Removes quality values and N values from the dna
preprocess()
{
	echo "*** Preprocessing ***"
	sed -n '2~4p' data/$basename/input.fastq > data/$basename/input.dna 
	python util/remove_N.py data/$basename/input.dna
}

compress()
{
	g++ src/cpp/noisy/matchsort7_v6.cpp -std=c++11 -o src/reorder_noisy.out
	mkdir -p data/$basename/output 
	./src/reorder_noisy.out data/$basename
	python src/encodernoisy.py data/$basename
}

decompress()
{
	echo "nothing right now"
}

#Process the arguments
while getopts hfpcd opt
do
   case "$opt" in
	h) usage; exit 1;;
	f) download;; 
	p) preprocess;;
	c) compress;;
	d) decompress;;
	?) usage; exit;;
   esac
done