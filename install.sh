#!/bin/bash

# Create dirs
mkdir -p data

pip install distance biopython

# Compile the reordering code
g++ src/cpp/noisy/matchsort7_v6.cpp -std=c++11 -o src/reorder_noisy.out 

