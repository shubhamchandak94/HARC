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
g++ src/preprocess_pe.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/preprocess_pe.out
g++ src/decoder.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/decoder.out
g++ src/pack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pack_order.out
g++ src/pe_encode.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pe_encode.out
g++ src/unpack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/unpack_order.out
g++ src/reorder_quality.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/reorder_quality.out
g++ src/reorder_compress_quality_id_pe.cpp src/ID_compression/src/*.c src/qvz/src/*.c -O3 -march=native -fopenmp -Isrc/ID_compression/include -Isrc/qvz/include -DLINUX -std=c++11 -o bin/reorder_compress_quality_id_pe.out
g++ src/decompress_quality_id_pe.cpp src/ID_compression/src/*.c src/qvz/src/*.c -O3 -march=native -fopenmp -Isrc/ID_compression/include -Isrc/qvz/include -std=c++11 -DLINUX -o bin/decompress_quality_id_pe.out
g++ src/pe_decode.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pe_decode.out

mkdir -p bin/reorder
for bitset_size in 64 128 192 256 320 384 448 512
do
	max_read_len=$(($bitset_size/2))
	echo "#define MAX_READ_LEN $max_read_len" > src/config.h
	g++ src/reorder.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/reorder/reorder_$bitset_size".out"
done

mkdir -p bin/encoder
mkdir -p bin/decoder_preserve
mkdir -p bin/decoder_pe
for bitset_size in 64 128 192 256 320 384 448 512 576 640 704 768
do
	max_read_len=$(($bitset_size/3))
	echo "#define MAX_READ_LEN $max_read_len" > src/config.h
	g++ src/encoder.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/encoder/encoder_$bitset_size".out"
	g++ src/decoder_preserve.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/decoder_preserve/decoder_preserve_$bitset_size".out"
	g++ src/decoder_pe.cpp -march=native -O3 -fopenmp -lpthread -std=c++11 -o bin/decoder_pe/decoder_pe_$bitset_size".out"
done
rm src/config.h
