#!/bin/bash

# Create dirs
mkdir -p data

#libbsc for read_seq compression
rm -rf src/libbsc
git clone https://github.com/shubhamchandak94/libbsc.git src/libbsc
(cd src/libbsc && make)

#pip install --user distance biopython joblib tqdm

#Compilation of some files
g++ src/preprocess.cpp -O3 -march=native -fopenmp -std=c++11 -o src/preprocess.out
g++ src/decoder.cpp -O3 -march=native -fopenmp -std=c++11 -o src/decoder.out
g++ src/pack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o src/pack_order.out
g++ src/pe_encode.cpp -O3 -march=native -fopenmp -std=c++11 -o src/pe_encode.out
g++ src/unpack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o src/unpack_order.out
g++ src/merge_N.cpp -O3 -march=native -std=c++11 -o src/merge_N.out
