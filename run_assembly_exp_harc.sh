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
for dataset in SRR823377
do
	
	prefix=$dirname"/"$dataset
	#### Log folders for the experiment
	log_folder=$prefix.log_folder
	mkdir $log_folder
	assembly_output_folder=$prefix.assembly_output
	mkdir $assembly_output_folder
	
	ref_fasta=$dirname"/"$dataset".genome.fasta"

#	## Run Original File
#	suffix="original_only_assembly"
#	tmp_folder=$assembly_output_folder/$suffix
#	/usr/bin/time -v python $spades_tool -t 32 --only-assembler -1 $prefix"_1.fastq" -2 $prefix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#	quast_output=$log_folder/$suffix
#	python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#	
#	### Run on Clean dataset
#	#suffix="clean_assembly_only"
#	#tmp_folder=$assembly_output_folder/$suffix
#	#/usr/bin/time -v python $spades_tool -t 32 --only-assembler -1 $prefix"_1.clean.fastq" -2 $prefix"_2.clean.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#	#quast_output=$log_folder/$suffix
#	#python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#	
#	## Run Original File with inbuilt denoising
#	suffix="original_bayeshammer"
#	tmp_folder=$assembly_output_folder/$suffix
#	/usr/bin/time -v python $spades_tool -t 32 -1 $prefix"_1.fastq" -2 $prefix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#	quast_output=$log_folder/$suffix
#	python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#
#	
#	## Run Reckoner denoised file
#	suffix="using_reckoner"
#	/usr/bin/time -v $reckoner_tool -prefix $dirname -threads 32 $prefix"_1.fastq" $prefix"_2.fastq" | tee $log_folder/$suffix".reckoner_log"
#	mv $prefix"_1.corrected.fastq" $prefix"_1.reckoner.fastq"
#	mv $prefix"_2.corrected.fastq" $prefix"_2.reckoner.fastq"	
#	tmp_folder=$assembly_output_folder/$suffix
#        /usr/bin/time -v python $spades_tool -t 32 --only-assembler -1 $prefix"_1.reckoner.fastq" -2 $prefix"_2.reckoner.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#	quast_output=$log_folder/$suffix
#	python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#	
#	## Run Karect Denoised File
#	suffix="using_Karect"
#	/usr/bin/time -v $karect_tool -correct -inputfile=$prefix"_1.fastq" -inputfile=$prefix"_2.fastq" -resultdir=$dirname -celltype=diploid -matchtype=hamming -threads=32 | tee $log_folder/$suffix".karect_log"
#	mv $dirname"/karect_"$dataset"_1.fastq" $prefix"_1.karect.fastq"
#	mv $dirname"/karect_"$dataset"_2.fastq" $prefix"_2.karect.fastq"	
#	tmp_folder=$assembly_output_folder/$suffix
#	/usr/bin/time -v python $spades_tool -t 32 --only-assembler -1 $prefix"_1.karect.fastq" -2 $prefix"_2.karect.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#	quast_output=$log_folder/$suffix
#	python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix

	
	## HARC Experiments
#	denoise_type="adaptive_denoising"
#	for n in 50 100 150
#	do
#		for r in 2 4
#		do
#			for a in 4 8
#			do
#				suffix="harc_"$denoise_type"_n_"$n"_r_"$r"_a_"$a
#				./harc -c $prefix".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $log_folder/$suffix
#				./harc -d $prefix".harc" -p
#				python util/combine_into_fastq.py $prefix $suffix 
#				tmp_folder=$assembly_output_folder/$suffix
#				/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#				quast_output=$log_folder/$suffix".quast_output"
#				python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#			done
#		done
#	done
#
	denoise_type="adaptive_denoising"
	for n in 1 2 3 4 5 6 7 8 9 
	do
		for s in 1
		do
		for r in 4
		do
			for a in 8
			do
				suffix="harc_"$denoise_type"_n_"$n"_r_"$r"_a_"$a
				./harc -c $prefix".fastq" -p -u $denoise_type -n $n -s $s -r $r -a $a -t 16 | tee $log_folder/$suffix
				./harc -d $prefix".harc" -p
				python util/combine_into_fastq.py $prefix $suffix 
				tmp_folder=$assembly_output_folder/$suffix
				/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
				quast_output=$log_folder/$suffix".quast_output"
				python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
			done
		done
		done
	done

#	denoise_type="counts_and_quality"
#	for n in 0.1
#	do
#		for s in 38
#		do
#		for r in 4
#		do
#			for a in 8
#			do
#				suffix="harc_"$denoise_type"_n_"$n"_s_"$s"_r_"$r"_a_"$a
#				./harc -c $prefix".fastq" -p -u $denoise_type -n $n -s $s -r $r -a $a -t 16 | tee $log_folder/$suffix
#				./harc -d $prefix".harc" -p
#				python util/combine_into_fastq.py $prefix $suffix 
#				tmp_folder=$assembly_output_folder/$suffix
#				/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#				quast_output=$log_folder/$suffix".quast_output"
#				python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#			done
#		done
#		done
#	done
#	denoise_type="counts_and_quality"
#	for n in 0.1
#	do
#		for s in 30
#		do
#		for r in 2
#		do
#			for a in 4 8
#			do
#				suffix="harc_"$denoise_type"_n_"$n"_s_"$s"_r_"$r"_a_"$a
#				./harc -c $prefix".fastq" -p -u $denoise_type -n $n -s $s -r $r -a $a -t 16 | tee $log_folder/$suffix
#				./harc -d $prefix".harc" -p
#				python util/combine_into_fastq.py $prefix $suffix 
#				tmp_folder=$assembly_output_folder/$suffix
#				/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#				quast_output=$log_folder/$suffix".quast_output"
#				python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#			done
#		done
#		done
#	done
#	denoise_type="only_using_counts"
#	for n in 0.1 0.2 0.5
#	do
#		for r in 2 4
#		do
#			for a in 8
#			do
#				suffix="harc_"$denoise_type"_n_"$n"_r_"$r"_a_"$a
#				./harc -c $prefix".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $log_folder/$suffix
#				./harc -d $prefix".harc" -p
#				python util/combine_into_fastq.py $prefix $suffix 
#				tmp_folder=$assembly_output_folder/$suffix
#				/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#				quast_output=$log_folder/$suffix".quast_output"
#				python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#			done
#		done
#	done
#
#	denoise_type="quality_threshold"
#	for n in 30 35
#	do
#		for r in 2 4
#		do
#			for a in 8
#			do
#				suffix="harc_"$denoise_type"_n_"$n"_r_"$r"_a_"$a
#				./harc -c $prefix".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $log_folder/$suffix
#				./harc -d $prefix".harc" -p
#				python util/combine_into_fastq.py $prefix $suffix 
#				tmp_folder=$assembly_output_folder/$suffix
#				/usr/bin/time -v python $spades_tool --only-assembler -t 32 -1 $prefix"."$suffix"_1.fastq" -2 $prefix"."$suffix"_2.fastq" -o $tmp_folder | tee $log_folder/$suffix".spades_log"
#				quast_output=$log_folder/$suffix".quast_output"
#				python $quast_tool $tmp_folder/contigs.fasta -R $ref_fasta -o $quast_output --plots-format pdf --labels $suffix
#			done
#		done
#	done
#
done
