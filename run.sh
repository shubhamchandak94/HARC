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
   -e      Compute Entropy Bound       
 		
EOF
}


# TODO: Make this more general later
download()
{
   
	echo "*** FASTQ Sequences being downloaded ***"
	mkdir -p data/$basename
	mkdir -p data/$basename/output
	if [[ -z "${URL_2}" ]]; then
		wget -O data/$basename/file_1.fastq.gz $URL_1
		gunzip data/$basename/file_1.fastq.gz 
		mv data/$basename/file_1.fastq data/$basename/input.fastq
	else
		echo "Combining 2 FASTQ files"
		wget -O data/$basename/file_1.fastq.gz $URL_1
		wget -O data/$basename/file_2.fastq.gz $URL_2
		gunzip data/$basename/file_1.fastq.gz data/$basename/file_2.fastq.gz
		cat data/$basename/file_1.fastq data/$basename/file_2.fastq > data/$basename/input.fastq
	fi
}

# Removes quality values and N values from the dna and writes reads with N to a separate file
preprocess()
{
	echo "*** Preprocessing ***"
	sed -n '2~4p' data/$basename/input.fastq > data/$basename/input.dna
	sed -n '4~4p' data/$basename/input.fastq > data/$basename/input.quality  
	python util/remove_N.py data/$basename/input.dna
	sort -o data/$basename/input_N.dna data/$basename/input_N.dna
}

generateConfig()
{
	echo "#define maxmatch $maxmatch" > src/cpp/noisy/config.h
	echo "#define thresh $thresh" >> src/cpp/noisy/config.h
	echo "#define numdict $numdict" >> src/cpp/noisy/config.h
	if [  ${dict1start+x} ]; then echo "#define dict1_start $dict1start" >> src/cpp/noisy/config.h; fi
	if [  ${dict1end+x} ]; then echo "#define dict1_end $dict1end" >> src/cpp/noisy/config.h; fi
	if [  ${dict2start+x} ]; then echo "#define dict2_start $dict2start" >> src/cpp/noisy/config.h; fi
	if [  ${dict2end+x} ]; then echo "#define dict2_end $dict2end" >> src/cpp/noisy/config.h; fi
	if [  ${dict3start+x} ]; then echo "#define dict3_start $dict3start" >> src/cpp/noisy/config.h; fi
	if [  ${dict3end+x} ]; then echo "#define dict3_end $dict3end" >> src/cpp/noisy/config.h; fi
	if [  ${dict4start+x} ]; then echo "#define dict4_start $dict4start" >> src/cpp/noisy/config.h; fi
	if [  ${dict4end+x} ]; then echo "#define dict4_end $dict4end" >> src/cpp/noisy/config.h; fi
	#readlen="$(wc -L < data/$basename/input_clean.dna)"
	readlen="$(head data/$basename/input_clean.dna | wc -L)"
	echo "#define readlen $readlen" >> src/cpp/noisy/config.h
	echo "#define num_thr $num_thr" >> src/cpp/noisy/config.h
}
compress()
{
	g++ src/cpp/noisy/matchsort7_v13.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o src/reorder_noisy.out
	mkdir -p data/$basename/output 
	./src/reorder_noisy.out data/$basename
	cp data/$basename/input_N.dna data/$basename/output/input_N.dna
	g++ src/cpp/noisy/encoder.cpp -march=native -O3 -fopenmp -std=c++11 -o src/encoder.out
	./src/encoder.out data/$basename

	# add > to top of read_seq.txt (needed for MFCompress)	
	#echo ">" > data/$basename/output/seq_header.txt
	#cat data/$basename/output/seq_header.txt data/$basename/output/read_seq.txt > data/$basename/output/read_seq.txt.1
	#mv data/$basename/output/read_seq.txt.1 data/$basename/output/read_seq.txt
	
	#remove temporary files
	rm data/$basename/output/temp.dna
	rm data/$basename/output/tempflag.txt
	rm data/$basename/output/temppos.txt
       	
	#compress and create tarball
	7z a data/$basename/output/read_pos.txt.7z data/$basename/output/read_pos.txt
	7z a data/$basename/output/read_noise.txt.7z data/$basename/output/read_noise.txt
	7z a data/$basename/output/read_noisepos.txt.7z data/$basename/output/read_noisepos.txt
	7z a data/$basename/output/input_N.dna.7z data/$basename/output/input_N.dna
	7z a data/$basename/output/read_meta.txt.7z data/$basename/output/read_meta.txt
	7z a data/$basename/output/read_rev.txt.7z data/$basename/output/read_rev.txt
	7z a -m0=PPMd -mo=12 -mmem=31 data/$basename/output/read_seq.txt.7z data/$basename/output/read_seq.txt
	#./util/MFCompress/MFCompressC data/$basename/output/read_seq.txt
	rm data/$basename/output/*.txt data/$basename/output/*.dna  data/$basename/output/*.bin
	tar -cf data/$basename/output.tar data/$basename/output
	rm -r data/$basename/output/
}

decompress()
{
	echo "Decompression ..."
	tar -xf data/$basename/output.tar
	7z e data/$basename/output/read_pos.txt.7z -odata/$basename/output/
	7z e data/$basename/output/read_noise.txt.7z -odata/$basename/output/
	7z e data/$basename/output/read_noisepos.txt.7z -odata/$basename/output/
	7z e data/$basename/output/input_N.dna.7z -odata/$basename/output/
	7z e data/$basename/output/read_meta.txt.7z -odata/$basename/output/
	7z e data/$basename/output/read_rev.txt.7z -odata/$basename/output/
	7z e data/$basename/output/read_seq.txt.7z -odata/$basename/output/
	#./util/MFCompress/MFCompressD data/$basename/output/read_seq.txt.mfc
	#tr -d '\r\n>' < data/$basename/output/read_seq.txt.mfc.d > data/$basename/output/read_seq.txt
	#rm data/$basename/output/read_seq.txt.mfc.d
#	python src/decodernoisy.py data/$basename
	./src/decoder.out data/$basename
}

compute_entropy()
{
    mkdir -p logs
	echo "Downloading the FASTA File"
	wget -O data/$basename/genome_fasta.fa.gz $URL_genome
	gunzip data/$basename/genome_fasta.fa.gz
	
    echo "Computing FASTA File entropy"
    ./util/MFCompress/MFCompressC -3 data/$basename/genome_fasta.fa

	echo "computing Noise entropy"
        cp config.py logs/"$basename"_entropy_computation.log
	python util/compute_entropy.py data/$basename/input.quality data/$basename/genome_fasta.fa.mfc data/$basename/genome_fasta.fa data/$basename/output.tar | tee -a logs/"$basename"_entropy_computation.log

}
#Process the arguments
while getopts hfpcdge opt
do
   case "$opt" in
	h) usage; exit 1;;
	f) download;; 
	p) preprocess;;
	g) generateConfig;;
	c) compress;;
	d) decompress;;
	e) compute_entropy;;
	?) usage; exit;;
   esac
done
