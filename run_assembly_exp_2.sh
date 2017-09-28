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

## Run Original File
python $spades_tool -t 16 --only-assembler -1 $prefix"_1.fq" -2 $prefix"_2.fq" -o $tmp_folder
quast_output=$log_folder/quastreport_original_only_assembler
python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels original_only_assembler

## Run Original File with spades_denoising
python $spades_tool -t 16  -1 $prefix"_1.fq" -2 $prefix"_2.fq" -o $tmp_folder
quast_output=$log_folder/quastreport_original_with_spades_denoising
python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels original_with_spades_denoising

## Run On Clean file
python util/combine_into_fastq_clean.py $prefix"_1" 1
python util/combine_into_fastq_clean.py $prefix"_2" 2
python $spades_tool -t 16 --only-assembler -1 $prefix"_1.clean.fastq" -2 $prefix"_2.clean.fastq" -o $tmp_folder
quast_output=$log_folder/quastreport_clean
python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels clean

#Run on different parameters of HARC
for n in 0.05 0.1 0.2
do
    for r in 2 4
    do
        for a in 2 4 8
        do
            ./harc -c $fastq_file -pq -r $r -a $a -n $n -t 16
            ./harc -d $harc_compressed_file -p
            param_string=_n$n._r$r._a$a
            python get_denoising_stats.py $prefix $log_folder/$param_string  >> $denoising_stats_log            
           

            python util/combine_into_fastq.py $prefix
            COUNT=$( cat $prefix.output.fastq | wc -l)
            COUNT=$((COUNT/2))
            split -l $COUNT $prefix".output.fastq" $tmp_folder/x  
            mv $tmp_folder/xaa $HARC_fastq_1
            mv $tmp_folder/xab $HARC_fastq_2

            python $spades_tool -t 16 --only-assembler -1 $HARC_fastq_1 -2 $HARC_fastq_2 -o $tmp_folder

            quast_output=$log_folder/quastreport_$param_string
            python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $param_string
        done
    done
done
        
