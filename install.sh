#!/bin/bash

# Create dirs
mkdir -p data

#libbsc for read_seq compression
git clone https://github.com/shubhamchandak94/libbsc.git src/libbsc
(cd src/libbsc && make)

#pip install --user distance biopython joblib tqdm



#Compiling preprocessor
g++ src/cpp/noisy/preprocess.cpp -O3 -march=native -std=c++11 -o src/preprocess.out
g++ src/cpp/noisy/preprocess_final.cpp -O3 -march=native -std=c++11 -o src/preprocess_final.out

#Compiling decoder
g++ src/cpp/noisy/decoder.cpp -O3 -march=native -std=c++11 -o src/decoder.out
g++ src/cpp/noisy/decoder_preserve.cpp -O3 -march=native -std=c++11 -o src/decoder_preserve.out
g++ src/cpp/noisy/merge_N.cpp -O3 -march=native -std=c++11 -o src/merge_N.out

