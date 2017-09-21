#!/bin/bash

mkdir $1.log_folder
denoising_stats_log=$1.log_folder/overall_denoising_stats.log_file

for n in 0.05 0.1 0.2 0.5 0.8
do
    for r in 2 4 8
    do
        for a in 0 2 4 8 
        do
            ./harc -c $1 -pq -n $n -r $r -a $a | tee $1.log_folder/n$n._r$r._a$a.log_file
            ./harc -d $1.harc -p
            python get_denoising_stats.py $1 $1.log_folder/_n$n._r$r._a$a >> $denoising_stats_log
        done 
    done
done

