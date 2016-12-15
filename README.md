# readcompression

##### Downloading datasets
###### Usual reads
```
wget -b ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_1.fastq.gz
wget -b ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_2.fastq.gz
gunzip SRR065390_1.fastq.gz SRR065390_2.fastq.gz
cat SRR065390_1.fastq SRR065390_2.fastq > SRR065390.fastq
```

For some datasets (e.g. SRR327342 and SRR870667), the two fastq files may have reads of different lengths

###### Metagenomics
B
```
wget -b http://public.genomics.org.cn/BGI/gutmeta/High_quality_reads/MH0001/081026/MH0001_081026_clean.1.fq.gz
wget -b http://public.genomics.org.cn/BGI/gutmeta/High_quality_reads/MH0001/081026/MH0001_081026_clean.2.fq.gz
gunzip MH0001_081026_clean.1.fq.gz MH0001_081026_clean.2.fq.gz
cat MH0001_081026_clean.1.fq MH0001_081026_clean.2.fq
```

###### Human genome
```
wget -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
```

Inside this file different chromosomes are demarcated.

##### Typical fastq format
```
@seq id
read
+
quality score
```

##### Running orcom
```
git clone https://github.com/lrog/orcom.git
cd orcom
make boost
cd bin
./orcom_bin e -iPATH/SRR065390.fastq -oPATH/SRR065390.bin
./orcom_pack e -iPATH/SRR065390.bin -oPATH/SRR065390.orcom
```

##### Getting orcom ordered file
```
./orcom_pack d -iPATH/SRR065390.orcom -oPATH/SRR065390.dna
```

The dna file is of the form:
```
read1
read2
..
```

##### Converting Fastq file to dna file (needed for our code)
```
sed -n '2~4p' SRR065390.fastq > SRR065390.dna &
```

##### Removing reads containing 'N' (python code)
```python
fout = open('SRR065390_clean.dna','w')
with open('SRR065390.dna','r') as f:
  for line in f:
     if 'N' not in line:
          fout.write(line)
```

##### Converting dna file to fastq file (with fake seq id and quality scores) (python code)
```python
fout = open('SRR065390_clean.fastq','w')
with open('SRR065390_clean.dna','r') as f:
  for line in f:
    fout.write('@\n'+line+'+\n'+line)
```

##### Compiling and running C++ codes
```
g++ matchsort7.cpp -std=c++11 -o a.out
./a.out
```
##### Using Google sparsehashmap
###### Installing
```
git clone https://github.com/sparsehash/sparsehash.git
cd sparsehash
./configure
sudo make install
```
###### Using
In the code replace 
```cpp
#include<unordered_map>
```
by
```cpp
#include<sparsehash/sparse_hash_map>
```
and replace the data type
```cpp
std::unordered_map<>
```
by
```cpp
google::sparse_hash_map<>
```
There should be no need to change anything else. However note that sparsehashmap does not seem to work with bitset as the index. One way to get around this (if the index length is not too large) is to use the bitset::to_ulong function.

##### Some python packages
```
sudo pip install distance
sudo pip install biopython
sudo pip install joblib
```

##### Running Python code
```
python matchsortnoisy8.py
python packernoisy2.py
```

##### Calculating Number of Hard reads
```
grep 0 -o read_flag.txt | wc -l
```

##### Calculating Number of Singleton reads (i.e. 00's in the flag file) 
```
tr -cs 0 '\012' < read_flag.txt | awk '/00/{n += length - 1}; END {print n+0}'
```

##### Generating reads using gen_fastq (from orcom repo)
###### Error-free reads (with reverse complementation) - 35M reads of length 100 from chrom 22
```
cd gen_fastq
make
./gen_fastq 35000000 100 PATH/chrom22clean.fasta PATH/chrom22_reads.fastq
```

###### Reads with 1% uniform substitution rate (each base is equally likely to be changed to any of the 4 possibilities e.g. A - C,T,G,N) (with reverse complementation) - 35M reads of length 100 from chrom 22
```
./gen_fastq 35000000 100 PATH/chrom22clean.fasta PATH/chrom22_reads.fastq -e
```

B
To get reads without reverse complementation, replace gen_fastq by gen_fastq_noRC above.
