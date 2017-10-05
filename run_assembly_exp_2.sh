#!/bin/bash

### Inputs
prefix=$1
ref_fasta=$2

#### Log folders for the experiment
log_folder=$prefix.log_folder
mkdir $log_folder
tmp_folder=$prefix.temp_folder
mkdir $tmp_folder
denoising_stats_log=$log_folder/overall_denoising_stats.log_file

##### TOOLS to use
spades_tool=../SPAdes-3.11.0-Linux/bin/spades.py
quast_tool=../quast-4.5/quast.py


fastq_file=$prefix.fastq
harc_compressed_file=$prefix".harc"
HARC_fastq_1=$prefix"_1.HARC.fastq"
HARC_fastq_2=$prefix"_2.HARC.fastq"


#Run on different parameters of HARC
for n in 0.05 0.1
do
    for r in 2 4
    do
        for a in 0 2 4 8
        do
            ./harc -c $fastq_file -pq -r $r -a $a -n $n -t 8
            ./harc -d $harc_compressed_file -p
            param_string=_n$n._r$r._a$a
            python get_denoising_stats.py $prefix $log_folder/$param_string  >> $denoising_stats_log  
        done
    done
done
        
