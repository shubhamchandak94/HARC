import sys
import os
import numpy as np
from tqdm import tqdm
from itertools import izip

basename = sys.argv[1]
pair=sys.argv[2]
infile_reads = basename + ".clean"
outfile_fastq = basename + ".clean.fastq"
readlen = 100

f_reads = open(infile_reads,'r')
f_out_fastq = open(outfile_fastq,'w')

separator = "+\n"

for i,read in enumerate(f_reads):
	readid = "@read_" + str(i)+"/" + pair + "\n"
	quality= "I"*(len(read)-1) + "\n"
	f_out_fastq.write(readid)
	f_out_fastq.write(read)
	f_out_fastq.write(separator)
	f_out_fastq.write(quality)
