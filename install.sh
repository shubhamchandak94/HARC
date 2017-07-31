#!/bin/bash

# Create dirs
mkdir -p data

#libbsc for read_seq compression
git clone https://github.com/shubhamchandak94/libbsc.git src/libbsc
(cd src/libbsc && make)

#pip install --user distance biopython joblib tqdm



#Compilation of some files
g++ src/preprocess.cpp -O3 -march=native -std=c++11 -o src/preprocess.out
g++ src/preprocess_quality.cpp -O3 -march=native -std=c++11 -o src/preprocess_quality.out
g++ src/decoder.cpp -O3 -march=native -std=c++11 -o src/decoder.out
g++ src/merge_N.cpp -O3 -march=native -std=c++11 -o src/merge_N.out
g++ src/merge_quality_N.cpp -O3 -march=native -std=c++11 -o src/merge_quality_N.out
