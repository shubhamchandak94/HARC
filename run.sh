#!/bin/bash

# Download data
# DATA_DIR="data/"
source config.py


usage()
{
cat << EOF
usage: $0 options


OPTIONS:
   -h      Show this message
   -f      download relevant files
   -p      preprocess
   -g      generateConfig
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

# Removes quality values and N values from the dna and writes reads with N to a separate file
preprocess()
{
	echo "*** Preprocessing ***"
	sed -n '2~4p' data/$basename/input.fastq > data/$basename/input.dna 
	python util/remove_N.py data/$basename/input.dna
	sort -o data/$basename/input_N.dna data/$basename/input_N.dna
}

generateConfig()
{
	echo "#define maxmatch $maxmatch" > src/cpp/noisy/config.h
	echo "#define thresh $thresh" >> src/cpp/noisy/config.h
	echo "#define numdict $numdict" >> src/cpp/noisy/config.h
	echo "#define dict1_start $dict1start" >> src/cpp/noisy/config.h
	echo "#define dict1_end $dict1end" >> src/cpp/noisy/config.h
	echo "#define dict2_start $dict2start" >> src/cpp/noisy/config.h
	echo "#define dict2_end $dict2end" >> src/cpp/noisy/config.h
	readlen="$(wc -L < data/$basename/input_clean.dna)"
	echo "#define readlen $readlen" >> src/cpp/noisy/config.h
}
compress()
{
	g++ src/cpp/noisy/matchsort7_v8.cpp -Isrc/cpp/noisy/sparsepp/ -std=c++11 -o src/reorder_noisy.out
	mkdir -p data/$basename/output 
	./src/reorder_noisy.out data/$basename
	split -a 4 -d -l $chunksize data/$basename/output/temp.dna data/$basename/output/temp.dna.
	split -a 4 -d -b $chunksize data/$basename/output/tempflag.txt data/$basename/output/tempflag.txt.
	python src/encodernoisy_parallel.py data/$basename
	cat data/$basename/output/read_seq.txt.* > data/$basename/output/read_seq.txt
	cat data/$basename/output/read_pos.txt.* > data/$basename/output/read_pos.txt
	cat data/$basename/output/read_noise.txt.* > data/$basename/output/read_noise.txt
	cat data/$basename/output/read_noisepos.txt.* > data/$basename/output/read_noisepos.txt
	cp data/$basename/input_N.dna data/$basename/output/input_N.dna
       # remove unwanted files
	rm data/$basename/output/temp.dna*
	rm data/$basename/output/tempflag.txt*
	rm data/$basename/output/read*txt.*
       #create tarball
	tar -cf data/$basename/output.tar data/$basename/output
	xz -f data/$basename/output.tar
	rm -r data/$basename/output/
}

decompress()
{
	echo "Decompression ..."
	xz -dkf data/$basename/output.tar.xz
	tar -xf data/$basename/output.tar
	rm data/$basename/output.tar
	python src/decodernoisy.py data/$basename
}

#Process the arguments
while getopts hfpcdg opt
do
   case "$opt" in
	h) usage; exit 1;;
	f) download;; 
	p) preprocess;;
	g) generateConfig;;
	c) compress;;
	d) decompress;;
	?) usage; exit;;
   esac
done
