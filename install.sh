#!/bin/bash

# Create dirs
mkdir -p data

#libbsc for read_seq compression
rm -rf src/libbsc
git clone https://github.com/shubhamchandak94/libbsc.git src/libbsc
(cd src/libbsc && make)
cp src/libbsc/bsc bin/

#pip install --user distance biopython joblib tqdm

#Compilation of some files
g++ src/preprocess.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/preprocess.out
g++ src/decoder.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/decoder.out
g++ src/pack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pack_order.out
g++ src/pe_encode.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pe_encode.out
g++ src/unpack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/unpack_order.out

for (( i=1; i<=8; i++))
do
	j=$(($i*32))
	mkdir -p bin/len_$j
	echo "#define MAX_READ_LEN $j" > src/config.h
	g++ src/reorder.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/len_$j/reorder.out
	g++ src/encoder.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/len_$j/encoder.out
	g++ src/decoder_preserve.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/len_$j/decoder_preserve.out
	g++ src/decoder_pe.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/len_$j/decoder_pe.out
done
rm src/config.h
