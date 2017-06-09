# readcompression

Compression of genomic reads in FASTQ format. Compresses only the read sequences. Achieves near-optimal compression ratios and fast decompression. Supports upto 4.29 Billion fixed-length reads with lengths at most 256. Requires around 50 Bytes of RAM/read for read length 100. The algorithm requires C++11 and g++ compiler and works on Linux.

#### Installation
```bash
git clone https://github.com/shubhamchandak94/readcompression.git
cd readcompression
./install.sh
```

#### Usage
Compression - compresses FASTQ reads. Output written to .tar file
```bash
./run_default.sh -c PATH_TO_FASTQ [-p] [-t NUM_THREADS]
```
-p = Preserve order of reads (compression ratio 2-4x worse if order preserved)

-t NUM_THREADS - default 8

Decompression - decompresses reads. Output written to .dna.d file
```bash
./run_default.sh -d PATH_TO_TAR [-p] [-t NUM_THREADS]
```
-p = Get reads in original order (slower). Only applicable if -p was used during compression.

-t NUM_THREADS - default 8

Help (this message)
```bash
./run_default.sh -h
```
##### Downloading datasets
###### Usual reads
```bash
wget -b ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_1.fastq.gz
wget -b ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_2.fastq.gz
gunzip SRR065390_1.fastq.gz SRR065390_2.fastq.gz
cat SRR065390_1.fastq SRR065390_2.fastq > SRR065390.fastq
```

For some datasets (e.g. SRR327342 and SRR870667), the two fastq files may have reads of different lengths

###### Metagenomics data
```bash
wget -b http://public.genomics.org.cn/BGI/gutmeta/High_quality_reads/MH0001/081026/MH0001_081026_clean.1.fq.gz
wget -b http://public.genomics.org.cn/BGI/gutmeta/High_quality_reads/MH0001/081026/MH0001_081026_clean.2.fq.gz
gunzip MH0001_081026_clean.1.fq.gz MH0001_081026_clean.2.fq.gz
cat MH0001_081026_clean.1.fq MH0001_081026_clean.2.fq
```

###### Human genome (hg19 - for generating simulated reads)
```bash
wget -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
```

Inside this file different chromosomes are demarcated.

##### Generating reads using gen_fastq (from orcom repo)
###### Error-free reads (without reverse complementation) - 35M reads of length 100 from chrom 22
```bash
cd util/gen_fastq_noRC
make
./gen_fastq_noRC 35000000 100 PATH/chrom22.fasta PATH/chrom22_reads.fastq
```

###### Reads with 1% uniform substitution rate (each base is equally likely to be changed to any of the 4 possibilities e.g. A - C,T,G,N) (without reverse complementation) - 35M reads of length 100 from chrom 22
```bash
./gen_fastq_noRC 35000000 100 PATH/chrom22.fasta PATH/chrom22_reads.fastq -e
```

##### Typical fastq format
```
@seq id
read
+
quality score
```

#### Other compressors (for evaluation)

##### Installing & Running orcom (boost should be installed)
```bash
git clone https://github.com/lrog/orcom.git
cd orcom
make boost
cd bin
./orcom_bin e -iPATH/SRR065390_clean.fastq -oPATH/SRR065390_clean.bin
./orcom_pack e -iPATH/SRR065390_clean.bin -oPATH/SRR065390_clean.orcom
```

##### Getting orcom ordered file
```bash
./orcom_pack d -iPATH/SRR065390_clean.orcom -oPATH/SRR065390_clean.dna
```

The dna file is of the form:
```
read1
read2
..
```

##### Installing and running Leon (seq-only mode)
```bash
git clone https://github.com/GATB/leon.git
cd leon
sh INSTALL
./leon -file PATH_TO_FASTQ -c -seq-only -nb-cores NUM_THR
```

##### Leon Decompression (order preserved)
```bash
./leon -file PATH_TO_LEON -d -nb-cores NUM_THR
```

The fasta.d file is of the form:
```
>1
read1
>2
read2
..
```

<!---

##### Some python packages
```
sudo pip install distance
sudo pip install biopython
sudo pip install joblib
```

##### Running Proposed tool (noiseless)
First set the parameters at the top of the files.
```
g++ reordernoiseless.cpp -std=c++11 -o a.out
./a.out
python encodernoiseless.py
xz -k outfiles
```

##### Running Proposed tool (noisy)
First set the parameters at the top of the files.
```
g++ reordernoisy.cpp -std=c++11 -o a.out
./a.out
python encodernoisy.py
xz -k outfiles
```

##### Decompressing (noisy)
First set the parameters at the top of the files.
```
python decodernoisy.py
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
-->

<!--

##### Running Python code
```
python matchsortnoisy8.py
python packernoisy2.py
```
##### Calculating Number of Hard reads
```
grep 0 -o read_flag.txt | wc -l
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

##### Calculating Number of Singleton reads (i.e. 00's in the flag file) 
```
tr -cs 0 '\012' < read_flag.txt | awk '/00/{n += length - 1}; END {print n+0}'
```
##### Generating reads using gen_fastq (from orcom repo)
###### Error-free reads (without reverse complementation) - 35M reads of length 100 from chrom 22
```
cd gen_fastq_noRC
make
./gen_fastq_noRC 35000000 100 PATH/chrom22.fasta PATH/chrom22_reads.fastq
```

###### Reads with 1% uniform substitution rate (each base is equally likely to be changed to any of the 4 possibilities e.g. A - C,T,G,N) (without reverse complementation) - 35M reads of length 100 from chrom 22
```
./gen_fastq_noRC 35000000 100 PATH/chrom22.fasta PATH/chrom22_reads.fastq -e
```


##### C++ files
These are the more important files in the C++ folder. For the other files, see comments on top of those files. Note that the noisy files are currently unable to handle reads with N. Also read length is assumed to be constant for all codes.
###### Noiseless and no RC
1. matchsort2.cpp - reordering, parameters - matchlen, maxmatch. Generates outfile with reordered reads.

###### Noisy
1. matchsort3.cpp - v1 reordering described in the report. Tries to find mathces to the current read (no clean reference). Parameters - numdict, dictionary indices, maxmatch, thresh. Generates three files: outfile which has the reordered reads (some of which reverse complemented), outfileRC which has flags (0/1) to tell if the read has been reverse complemented, outfileflag which has the flags (0/1) to tell if the read is matched or not. This information helps packernoisy2_noN.py and packernoisy4_noN.py
2. matchsort6.cpp - has a recovery step after the reordering which tries to place the singleton reads before a matching read. thresh2 is the threshold for this second stage process. Other parameters and output files are same as matchsort3.
3. matchsort7.cpp - v2 reordering described in the report. Uses majority-based reference read for the reordering. Parameters and output files are same as matchsort3

##### Python files
These are the more important files in the python folder. For the other files, see comments on top of those files. Note that the some of the noisy files can handle reads with N. Also read length is assumed to be constant for all codes.
###### Noiseless and no RC
1. matchsort.py - reordering, parameters - matchlen, maxmatch. Generates outfile with reordered reads.
2. packer.py - encoding, takes infile with reordered reads and generates outfile_seq and outfile_flag. Parameter - maxmatch, neg (whether we want to look for matches with opposite shift as well.

###### Noisy for real data
1. matchpacknoisyRC.py - v2 reordering described in the report along with encoding. Uses majority-based reference read for the reordering. Parameters - numdict, dictionary indices, maxmatch, thresh. Directly generates encoding in the form of five files - outfile_seq, outfile_flag, outfile_rev, outfile_noise and outfile_noisepos. Works with reads containing N.
3. packernoisy2.py - encoding for reordered reads as described in report (uses reference). Input - infile containing reordered reads. Parameters - thresh, maxmatch. The reordering should generate the reverse complement flags separately, this does not consider RC. Produces 4 files - outfile_seq, outfile_flag, outfile_noise and outfile_noisepos. The current implementation is slow (due to findmajority function), see packernoisy2_noN.py for faster implementation.
4. packernoisy2_parallel.py - Parallel implementation of packernoisy2.py using joblib library - see comments on top of file. Leaves a blank line in each output file after each thread. Extra parameter - numthreads. Still slower than packernoisy2_noN.py.
5. packernoisy4.py - Similar to packernoisy2.py except to encoding of noise. Instead of storing the noisy base, we store cyclic shift using 1,2,3 or 4 (note that N is also a possibility). See function encodenoise for the exact encoding.

###### Encoders for C++ generated reordering (no reads with N)
1. packernoisy2_noN.py - Similar to packernoisy2.py but can't handle reads with N. Also an extra input file is needed - infile_flag which contains the flags (0/1 - unmatched/matched) produced by the C++ codes (matchsort3,6,7.cpp). Much faster than packernoisy2.py due to the flag file and the better written findmajority function.
2. packernoisy4_noN.py - Similar to packernoisy4.py but can't handle reads with N. Also an extra input file is needed - infile_flag which contains the flags (0/1 - unmatched/matched) produced by the C++ codes (matchsort3,6,7.cpp). Much faster than packernoisy4.py due to the flag file and the better written findmajority function.

###### Reordering for noisy simulated data without RC
1. matchsortnoisy2.py - Similar to v1 reordering described in the report. Tries to find matches to the current read (no clean reference). 4 dictionaries. Parameters - matchlen, maxmatch, thresh. Produces outfile with reordered reads. Works with reads containing N.
2. matchsortnoisy8.py - Similar to v2 reordering described in the report. Tries to find matches to the clean reference. 5 dictionaries. Parameters - matchlen, maxmatch, thresh. Produces outfile with reordered reads. Works with reads containing N.
-->
