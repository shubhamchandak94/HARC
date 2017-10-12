#!/bin/bash

#mkdir $1.log_folder
##
#denoise_type="adaptive_denoising"
#denoising_stats_log=$1.log_folder/overall_denoising_new_buildcontig.$denoise_type".log_file"
##
#for n in 100.0 
#do
#    for r in 2
#    do
#        for a in 8
#        do
#            ./harc -c $1".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n._r$r._a$a.log_file
#            ./harc -d $1.harc -p
#            ./util/get_denoising_stats.out $1 $1.log_folder/$denoise_type"_n"$n._r$r._a$a >> $denoising_stats_log
#        done 
#    done
#done

#denoise_type="counts_and_quality"
#denoising_stats_log=$1.log_folder/overall_denoising_new_buildcontig.$denoise_type".log_file"
#
#for n in 0.2
#do
#    for s in 30
#    do	
#	for r in 4
#    	do
#        	for a in 8
#       		do
#            		./harc -c $1".fastq" -p -u $denoise_type -n $n -s $s -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n."_s"$s._r$r._a$a.log_file
#            		./harc -d $1.harc -p
#		        ./util/get_denoising_stats.out $1 $1.log_folder/$denoise_type"_n"$n."_s"$s._r$r._a$a >> $denoising_stats_log
#		done
#        done 
#    done
#done

denoise_type="quality_threshold"
denoising_stats_log=$1.log_folder/overall_denoising_new_buildcontig.$denoise_type".log_file"

for n in 30
do
    for r in 2 
    do
        for a in 8
        do
            ./harc -c $1".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n._r$r._a$a.log_file
            ./harc -d $1.harc -p
	   ./util/get_denoising_stats.out $1 $1.log_folder/$denoise_type"_n"$n._r$r._a$a >> $denoising_stats_log
        done 
    done
done

#denoise_type="only_using_counts"
#denoising_stats_log=$1.log_folder/overall_denoising_new_buildcontig.$denoise_type".log_file"
#
#for n in 0.1
#do
#    for r in 4 
#    do
#        for a in 8
#        do
#            ./harc -c $1".fastq" -p -u $denoise_type -n $n -r $r -a $a -t 16 | tee $1.log_folder/$denoise_type"_n"$n._r$r._a$a.log_file
#            ./harc -d $1.harc -p
#	   ./util/get_denoising_stats.out $1 $1.log_folder/$denoise_type"_n"$n._r$r._a$a >> $denoising_stats_log
#        done 
#    done


#done
