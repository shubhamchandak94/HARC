#!/bin/bash

#mkdir $1.log_folder
#
denoise_type="adaptive_denoising"
denoising_stats_log=$1.log_folder/overall_denoising_corrected.$denoise_type".log_file"
#
for n in 50.0 80.0 100.0 150.0 200.0 500.0 
do
    for r in 2 4 
    do
        for a in 0 2 4 8 
        do
            ./harc -c $1".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n._r$r._a$a.log_file
            ./harc -d $1.harc -p
            python get_denoising_stats.py $1 $1.log_folder/$denoise_type"_n"$n._r$r._a$a >> $denoising_stats_log
        done 
    done
done

#denoise_type="counts_and_quality"
#denoising_stats_log=$1.log_folder/overall_denoising_harc20dict_corrected.$denoise_type".log_file"

#for n in 0.1 0.2
#do
#    for s in 20 30
#    do	
#	for r in 2 4
#    	do
#        	for a in 0
#       		do
#            		./harc -c $1".fastq" -p -u $denoise_type -n $n -s $s -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n."_s"$s._r$r._a$a.log_file
#            		./harc -d $1.harc -p
#            		python get_denoising_stats.py $1 $1.log_folder/$denoise_type"_n"$n."_s"$s._r$r._a$a >> $denoising_stats_log
#		done
#        done 
#    done
#done

#denoise_type="only_using_counts"
#denoising_stats_log=$1.log_folder/overall_denoising.$denoise_type".log_file"
#
#for n in 0.01 0.05 0.1 0.2 0.5 
#do
#    for r in 2 4 8 
#    do
#        for a in 0 2 4 8
#        do
#            ./harc -c $1".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n._r$r._a$a.log_file
#            ./harc -d $1.harc -p
#            python get_denoising_stats.py $1 $1.log_folder/$denoise_type"_n"$n._r$r._a$a >> $denoising_stats_log
#        done 
#    done
#done
