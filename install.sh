#!/bin/bash

# Create dirs
mkdir -p data

#libbsc for read_seq compression
git clone https://github.com/shubhamchandak94/libbsc.git src/libbsc
(cd src/libbsc && make)

#pip install --user distance biopython joblib tqdm

# Install MFCompress
wget  http://sweet.ua.pt/ap/software/mfcompress/MFCompress-linux64-1.01.tgz
tar -zxvf MFCompress-linux64-1.01.tgz
mkdir -p util/MFCompress
mv MFCompress-linux64-1.01/* util/MFCompress
rm -r MFCompress-linux64-1.01*

#Compiling preprocessor
g++ src/cpp/noisy/preprocess.cpp -O3 -march=native -std=c++11 -o src/preprocess.out

#Compiling decoder
g++ src/cpp/noisy/decoder.cpp -O3 -march=native -std=c++11 -o src/decoder.out
