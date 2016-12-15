import operator
from itertools import imap

#hamming2 function from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming2(str1, str2):
    #assert len(str1) == len(str2)
    #ne = str.__ne__  ## this is surprisingly slow
    ne = operator.ne
    return sum(imap(ne, str1, str2))

infile = "tempte.dna"
outfile_seq = "read_seq16.txt"
outfile_pos = "read_pos16.txt"
outfile_noise = "read_noise16.txt"
outfile_noisepos = "read_noisepos16.txt"

readlen = 100
minmatch = 28
thresh = 20 # maximum number of mismatches allowed 

f_seq = open(outfile_seq,'w')
f_pos = open(outfile_pos,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
k = 0
prev = 'A'*readlen
with open(infile,'r') as f:
	for line in f:
		k = k + 1
		if k%1000000 == 0:
			print str(k//1000000)+'M done'
		current = line.rstrip('\n')
		flag = 0
		for i in range(minmatch):
			if(hamming2(current[:(readlen-i)],prev[i:])<=thresh):
				f_pos.write('+')
				f_seq.write(current[(readlen-i+1):]+'\n')
				for j in range(readlen-i):
					if current[j]!=prev[i+j]:
						f_noise.write(current[j])
						f_noisepos.write("%02d"%j)
				prev = current
				f_noise.write('\n')
				flag = 1
				break
		
		if flag == 0:
			f_pos.write('0')
			f_seq.write(current+'\n')
			prev = current
f_seq.close()
f_pos.close()	
f_noise.close()
f_noisepos.close()

