#!/bin/bash

# Create dirs
mkdir -p data

pip install distance biopython joblib tqdm

# Install MFCompress
wget  http://sweet.ua.pt/ap/software/mfcompress/MFCompress-linux64-1.01.tgz
tar -zxvf MFCompress-linux64-1.01.tgz
mkdir -p util/MFCompress
mv MFCompress-linux64-1.01/* util/MFCompress 
rm -r MFCompress-linux64-1.01*
