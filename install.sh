#!/bin/bash

# Create dirs
mkdir -p data

# Install Sparsepp
git clone https://github.com/shubhamchandak94/sparsepp src/cpp/noisy/sparsepp

pip install distance biopython joblib tqdm

# Install MFCompress
wget  http://sweet.ua.pt/ap/software/mfcompress/MFCompress-linux64-1.01.tgz
tar -zxvf MFCompress-linux64-1.01.tgz
mkdir -p util/MFCompress
mv MFCompress-linux64-1.01/* util/MFCompress
chmod 741 ./util/MFCompress/MFCompressC
chmod 741 ./util/MFCompress/MFCompressD 
rm -r MFCompress-linux64-1.01*

#Compiling encoder and decoder
g++ src/cpp/noisy/encoder.cpp -O3 -march=native -std=c++11 -o src/encoder.out
g++ src/cpp/noisy/decoder.cpp -O3 -march=native -std=c++11 -o src/decoder.out
