#!/bin/bash
set -e
# Download data
# DATA_DIR="data/"


usage()
{
cat << EOF
usage: 
Compression
./run_default.sh -c PATH_TO_FASTQ -t NUM_THREADS[=8]

Decompression
./run_default.sh -d PATH_TO_TAR -t NUM_THREADS[=8]

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
	./src/preprocess.out $filename $pathname

	readlen="$(head $pathname/input_clean.dna | wc -L)"
	if (($readlen > 256));then
		echo "Maximum read length exceeded" 
		exit 1
	fi
	echo "#define maxmatch $((readlen/2))" > src/cpp/noisy/config.h
	echo "#define thresh 4" >> src/cpp/noisy/config.h
	echo "#define numdict 2" >> src/cpp/noisy/config.h
	echo "#define maxsearch 1000" >> src/cpp/noisy/config.h
	echo "#define dict1_start $(( readlen > 100 ? readlen/2-32 : readlen/2-readlen*32/100 ))" >> src/cpp/noisy/config.h
	echo "#define dict1_end $((readlen/2-1))" >> src/cpp/noisy/config.h
	echo "#define dict2_start $((readlen/2))" >> src/cpp/noisy/config.h
	echo "#define dict2_end $(( readlen > 100 ? readlen/2-1+32 : readlen/2-1+readlen*32/100 ))" >> src/cpp/noisy/config.h

	echo "#define readlen $readlen" >> src/cpp/noisy/config.h
	echo "#define num_thr $num_thr" >> src/cpp/noisy/config.h

	g++ src/cpp/noisy/matchsort7_v14.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o src/reorder_noisy.out
	mkdir -p $pathname/output/ 
	./src/reorder_noisy.out $pathname
	mv $pathname/input_N.dna $pathname/output/input_N.dna
	mv $pathname/read_order_N.bin $pathname/output/read_order_N.bin
	g++ src/cpp/noisy/encoder.cpp -march=native -O3 -fopenmp -std=c++11 -o src/encoder.out
	./src/encoder.out $pathname
	
	#remove temporary files
	rm $pathname/output/temp.dna
	rm $pathname/output/tempflag.txt
	rm $pathname/output/temppos.txt
       	rm $pathname/input.dna
	
	#compress and create tarball
	7z a $pathname/output/read_pos.txt.7z $pathname/output/read_pos.txt -mmt=$num_thr
	7z a $pathname/output/read_noise.txt.7z $pathname/output/read_noise.txt -mmt=$num_thr
	7z a $pathname/output/read_noisepos.txt.7z $pathname/output/read_noisepos.txt -mmt=$numt_thr
	7z a $pathname/output/input_N.dna.7z $pathname/output/input_N.dna -mmt=$numt_thr
	7z a $pathname/output/read_meta.txt.7z $pathname/output/read_meta.txt -mmt=$numt_thr
	7z a $pathname/output/read_rev.txt.7z $pathname/output/read_rev.txt -mmt=$num_thr
	./src/libbsc/bsc e $pathname/output/read_seq.txt $pathname/output/read_seq.txt.bsc -b512p -tT #-tT for single thread - uses too much memory in multi-threaded
	rm $pathname/output/*.txt $pathname/output/*.dna  $pathname/output/*.bin
	tar -cf $pathname/$(basename "$filename" .fastq).tar $pathname/output
	rm -r $pathname/output/
}

decompress()
{
	echo "Decompression ..."
	pathname=$(dirname $filename)
	tar -xf $filename -C $pathname
	7z e $pathname/output/read_pos.txt.7z -o$pathname/output/
	7z e $pathname/output/read_noise.txt.7z -o$pathname/output/
	7z e $pathname/output/read_noisepos.txt.7z -o$pathname/output/
	7z e $pathname/output/input_N.dna.7z -o$pathname/output/
	7z e $pathname/output/read_meta.txt.7z -o$pathname/output/
	7z e $pathname/output/read_rev.txt.7z -o$pathname/output/
	./src/libbsc/bsc d $pathname/output/read_seq.txt.bsc $pathname/output/read_seq.txt -tT
	./src/decoder.out $pathname
	rm -r $pathname/output/
}

#Initialize variables to default values.
num_thr=8

#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
 usage
fi

mode=''

while getopts ':c:d:t:' opt; do
  case "$opt" in
    c) [[ -n "$mode" ]] && usage || mode='c' && filename=$OPTARG;;
    d) [[ -n "$mode" ]] && usage || mode='d' && filename=$OPTARG;;
    t) num_thr=$OPTARG;;
    h) usage ;;
    \?) usage ;;
    *) usage ;;
  esac
done

if [[ $mode == 'c' ]];then
compress
else
decompress
fi;
