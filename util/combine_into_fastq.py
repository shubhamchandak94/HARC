import sys
import os
import numpy as np
from tqdm import tqdm
from itertools import izip

basename = sys.argv[1]
infile_reads = basename + ".dna.d"
infile_quality = basename + ".quality"
infile_readid = basename + ".readid"
outfile_fastq = basename + ".output.fastq"

f_reads = open(infile_reads,'r')
f_quality = open(infile_quality, 'r')
f_readid = open(infile_readid, 'r')
f_out_fastq = open(outfile_fastq,'w')

separator = "+\n"


for readid,read, quality in izip(f_readid, f_reads, f_quality):
	f_out_fastq.write(readid) 
	f_out_fastq.write(read)
	f_out_fastq.write(separator)
	f_out_fastq.write(quality)