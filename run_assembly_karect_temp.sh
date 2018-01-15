#!/bin/bash

### Inputs
#prefix=$1
#ref_fasta=$2
#denoising_stats_log=$log_folder/overall_denoising_stats.log_file

##### TOOLS to use
spades_tool=../SPAdes-3.11.0-Linux/bin/spades.py
quast_tool=../quast-4.5/quast.py
reckoner_tool=../RECKONER/bin/reckoner
karect_tool=../karect/karect

dirname=assembly_data/datasets
for dataset in ERR204808 
do
	prefix=$dirname"/"$dataset
	log_folder=$dirname.log_folder
	assembly_output_folder=$prefix.assembly_output	

	ref_fasta=$prefix".genome.fasta"
	## Run Reckoner denoised file
	#suffix="using_reckoner"
	#/usr/bin/time -v $reckoner_tool -prefix $dirname -threads 16 $prefix"_1.fastq" $prefix"_2.fastq" | tee $log_folder/$suffix".reckoner_log"
	#mv $prefix"_1.corrected.fastq" $prefix"_1.reckoner.fastq"
	#mv $prefix"_2.corrected.fastq" $prefix"_2.reckoner.fastq"	
	#tmp_folder=$assembly_output_folder/$suffix
        #/usr/bin/time -v python $spades_tool -t 16 --only-assembler -1 $prefix"_1.reckoner.fastq" -2 $prefix"_2.reckoner.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
	#quast_output=$log_folder/$suffix
	#python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
	
	## Run Karect Denoised File
	suffix="using_Karect_matchtype_edit"
	#/usr/bin/time -v $karect_tool -correct -inputfile=$prefix"_1.fastq" -inputfile=$prefix"_2.fastq" -resultdir=$dirname -celltype=diploid -matchtype=edit -threads=16 | tee $log_folder/$suffix".karect_log"
	#mv $dirname"/karect_"$dataset"_1.fastq" $prefix"_1.karect.fastq"
	#mv $dirname"/karect_"$dataset"_2.fastq" $prefix"_2.karect.fastq"	
	tmp_folder=$assembly_output_folder/$suffix
	/usr/bin/time -v python $spades_tool -t 16 --only-assembler -1 $prefix"_1.karect.fastq" -2 $prefix"_2.karect.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
	quast_output=$log_folder/$suffix
	python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
done
