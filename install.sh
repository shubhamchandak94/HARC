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
#<<<<<<< HEAD
#g++ src/preprocess_quality_id.cpp -O3 -march=native -fopenmp -std=c++11 -o src/preprocess_quality.out ##### Temporary change by KEDAR
#=======
#>>>>>>> 5a6551d005a1a95c49abd259bcadb8addc4039fc
g++ src/decoder.cpp -O3 -march=native -fopenmp -std=c++11 -o src/decoder.out
g++ src/pack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o src/pack_order.out
g++ src/unpack_order.cpp -O3 -march=native -fopenmp -std=c++11 -o src/unpack_order.out
g++ src/merge_N.cpp -O3 -march=native -std=c++11 -o src/merge_N.out
