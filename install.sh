#!/bin/bash

# Create dirs
mkdir -p data

#libbsc for read_seq compression
rm -rf src/libbsc
git clone https://github.com/shubhamchandak94/libbsc.git src/libbsc
(cd src/libbsc && make)

#pip install --user distance biopython joblib tqdm

#Compilation of some files
g++ src/preprocess.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/preprocess.out
g++ src/decoder.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/decoder.out
g++ src/pack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pack_order.out
g++ src/pe_encode.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pe_encode.out
g++ src/pe_decode.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/pe_decode.out
g++ src/unpack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o bin/unpack_order.out


