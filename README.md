# readcompression

##### Downloading datasets:
```
wget -b ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_1.fastq.gz
wget -b ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR065/SRR065390/SRR065390_2.fastq.gz
gunzip SRR065390_1.fastq.gz SRR065390_2.fastq.gz
cat SRR065390_1.fastq SRR065390_2.fastq > SRR065390.fastq
```

##### Typical fastq format:
```
@seq id
read
+
quality score
```

##### Running orcom:
```
git clone https://github.com/lrog/orcom.git
cd orcom
make boost
cd bin
./orcom_bin e -iPATH/SRR065390.fastq -oPATH/SRR065390.bin
./orcom_pack e -iPATH/SRR065390.bin -oPATH/SRR065390.orcom
```

##### Getting orcom ordered file:
```
./orcom_pack d -iPATH/SRR065390.orcom -oPATH/SRR065390.dna
```

The dna file is of the form:
```
read1
read2
..
```

##### Converting Fastq file to dna file (needed for our code):
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

##### Converting dna file to fastq file (with fake seq id and quality scores) (python code):
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
