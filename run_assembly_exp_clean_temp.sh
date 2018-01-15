#!/bin/bash

### Inputs
prefix=$1
ref_fasta=$2

#### Log folders for the experiment
log_folder=$prefix.log_folder
mkdir $log_folder
tmp_folder=$prefix.temp_folder_2
mkdir $tmp_folder
denoising_stats_log=$log_folder/overall_denoising_stats.log_file

##### TOOLS to use
spades_tool=../SPAdes-3.11.0-Linux/bin/spades.py
quast_tool=../quast-4.5/quast.py


fastq_file=$prefix.fastq
harc_compressed_file=$prefix".harc"
HARC_fastq_1=$prefix"_1.HARC.fastq"
HARC_fastq_2=$prefix"_2.HARC.fastq"


## Run On Clean file
#python util/combine_into_fastq_clean.py $prefix"_1" 1
#python util/combine_into_fastq_clean.py $prefix"_2" 2
python $spades_tool -t 16 --only-assembler -1 $prefix"_1.clean.fastq" -2 $prefix"_2.clean.fastq" -o $tmp_folder
quast_output=$log_folder/quastreport_clean
python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels clean

