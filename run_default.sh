#!/bin/bash
set -e
# Download data
# DATA_DIR="data/"


usage()
{
cat << EOF
HARC compression tool for genomic reads. Works on fixed length reads of length less than 256.

Usage: 
Compression - compresses FASTQ reads. Output written to .tar file
./run_default.sh -c PATH_TO_FASTQ [-p] [-t NUM_THREADS] [-q]
-p = Preserve order of reads
-t NUM_THREADS - default 8
-q = Write quality values to .quality file. Quality values are appropriately reordered if -p is not specified. 

Decompression - decompresses reads. Output written to .dna.d file
./run_default.sh -d PATH_TO_TAR [-p] [-t NUM_THREADS] [-m max_memory]
-p = Get reads in original order (slower). Only applicable if -p was used during compression.
-t NUM_THREADS - default 8
-m max_memory - controls memory-time tradeoff for decompression with -p. Specify max memory in GB (default 7 GB). Minimum of 3 GB memory is always used. Example: -m 10 for 10 GB maximum memory.

Help (this message)
./run_default.sh -h

EOF
exit 0
}

compress()
{
	pathname=$(dirname $filename)
	echo "*** Preprocessing ***"
	echo $filename
	if [[ $preserve_quality == "True" ]];then
		./src/preprocess_quality.out $filename $pathname
	else
		./src/preprocess.out $filename $pathname
	fi
	readlen="$(head $pathname/input_clean.dna | wc -L)"
	if (($readlen > 256));then
		echo "Maximum read length exceeded" 
		exit 1
	fi
	echo "#define maxmatch $((readlen/2))" > src/config.h
	echo "#define thresh 4" >> src/config.h
	echo "#define numdict 2" >> src/config.h
	echo "#define maxsearch 1000" >> src/config.h
	echo "#define dict1_start $(( readlen > 100 ? readlen/2-32 : readlen/2-readlen*32/100 ))" >> src/config.h
	echo "#define dict1_end $((readlen/2-1))" >> src/config.h
	echo "#define dict2_start $((readlen/2))" >> src/config.h
	echo "#define dict2_end $(( readlen > 100 ? readlen/2-1+32 : readlen/2-1+readlen*32/100 ))" >> src/config.h

	echo "#define readlen $readlen" >> src/config.h
	echo "#define num_thr $num_thr" >> src/config.h

	g++ src/reorder.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o src/reorder.out
	mkdir -p $pathname/output/ 
	./src/reorder.out $pathname
	mv $pathname/input_N.dna $pathname/output/input_N.dna
	mv $pathname/input_clean.dna $pathname/output/input_clean.dna
	mv $pathname/read_order_N.bin $pathname/output/read_order_N.bin
	if [[ $preserve_quality == "True" ]];then
		mv $pathname/input_N.quality $pathname/output/input_N.quality
		mv $pathname/input_clean.quality $pathname/output/input_clean.quality
	fi
	g++ src/encoder.cpp -march=native -O3 -fopenmp -std=c++11 -o src/encoder.out
	./src/encoder.out $pathname
	
	#remove temporary files
	rm $pathname/output/temp.dna
	rm $pathname/output/tempflag.txt
	rm $pathname/output/temppos.txt
	
	#compress and create tarball
	7z a $pathname/output/read_pos.txt.7z $pathname/output/read_pos.txt -mmt=$num_thr
	7z a $pathname/output/read_noise.txt.7z $pathname/output/read_noise.txt -mmt=$num_thr
	7z a $pathname/output/read_noisepos.txt.7z $pathname/output/read_noisepos.txt -mmt=$numt_thr
	7z a $pathname/output/input_N.dna.7z $pathname/output/input_N.dna -mmt=$numt_thr
	7z a $pathname/output/read_meta.txt.7z $pathname/output/read_meta.txt -mmt=$numt_thr
	7z a $pathname/output/read_rev.txt.7z $pathname/output/read_rev.txt -mmt=$num_thr
	./src/libbsc/bsc e $pathname/output/read_seq.txt $pathname/output/read_seq.txt.bsc -b512p -tT #-tT for single thread - uses too much memory in multi-threaded
	rm $pathname/output/*.txt $pathname/output/*.dna  
	if [[ $preserve_order == "True" ]];then
		7z a $pathname/output/read_order.bin.7z $pathname/output/read_order.bin -mmt=$num_thr
		7z a $pathname/output/read_order_N.bin.7z $pathname/output/read_order_N.bin -mmt=$num_thr
		if [[ $preserve_quality == "True" ]];then
			./src/merge_quality.out $pathname
			mv $pathname/output/output.quality $pathname/$(basename "$filename" .fastq).quality
		fi
	else
		if [[ $preserve_quality == "True" ]];then
			echo "Reordering quality values"
			g++ src/reorder_quality.cpp -march=native -O3 -std=c++11 -o src/reorder_quality.out
			./src/reorder_quality.out $pathname
			mv $pathname/output/output.quality $pathname/$(basename "$filename" .fastq).quality
		fi	
	fi
	rm $pathname/output/*.bin $pathname/output/*.quality
	tar -cf $pathname/$(basename "$filename" .fastq).tar -C $pathname/output .
	rm -r $pathname/output/
}

decompress()
{
	echo "Decompression ..."
	pathname=$(dirname $filename)
	mkdir -p $pathname/output
	tar -xf $filename -C $pathname/output
	if [[ $preserve_order == "True" ]];then
		if [ ! -f $pathname/output/read_order.bin.7z ];then
			echo "Not compressed using -p flag. Order cannot be restored"
			usage
			exit 1
		fi
	fi
	7z e $pathname/output/read_pos.txt.7z -o$pathname/output/
	7z e $pathname/output/read_noise.txt.7z -o$pathname/output/
	7z e $pathname/output/read_noisepos.txt.7z -o$pathname/output/
	7z e $pathname/output/input_N.dna.7z -o$pathname/output/
	7z e $pathname/output/read_meta.txt.7z -o$pathname/output/
	7z e $pathname/output/read_rev.txt.7z -o$pathname/output/
	./src/libbsc/bsc d $pathname/output/read_seq.txt.bsc $pathname/output/read_seq.txt -tT
	if [[ $preserve_order == "True" ]];then
		readlen=$( cat $pathname/output/read_meta.txt )
		echo "#define readlen $readlen" > src/config.h
		echo "#define MAX_BIN_SIZE $memory" >> src/config.h
		g++ src/decoder_preserve.cpp -O3 -march=native -std=c++11 -o src/decoder_preserve.out
		7z e $pathname/output/read_order.bin.7z -o$pathname/output/
		7z e $pathname/output/read_order_N.bin.7z -o$pathname/output/
		./src/decoder_preserve.out $pathname
		./src/merge_N.out $pathname
		echo "Done!"
	else
		./src/decoder.out $pathname
	fi
	rm -r $pathname/output/
	mv $pathname/output.dna $pathname/$(basename "$filename" .tar).dna.d
}

#Initialize variables to default values.
num_thr=8

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
 usage
fi

mode=''
preserve_order="False"
preserve_quality="False"
memory=''

while getopts ':c:d:t:m:pqh' opt; do
  case "$opt" in
    c) [[ -n "$mode" ]] && usage || mode='c' && filename=$OPTARG;;
    d) [[ -n "$mode" ]] && usage || mode='d' && filename=$OPTARG;;
    t) num_thr=$OPTARG;;
    m) memory=$OPTARG;;
    p) preserve_order="True";;
    q) preserve_quality="True";;	
    h) usage ;;
    \?) usage ;;
    *) usage ;;
  esac
done

if [[ $mode == 'c' ]];then
compress
elif [[ $mode == 'd' ]];then
decompress
else
echo "Either -c or -d required"
usage
exit 1
fi;
