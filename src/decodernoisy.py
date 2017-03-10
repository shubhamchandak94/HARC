#Decode the reordered reads. readlen is the length of read (assumed constant). the infile's are the 5 encoded files and
#outfile is the file to which the reordered reads will be written.

from Bio.Seq import Seq
import sys
import os
<<<<<<< HEAD
#sys.path.append("..")
#import config.ini

basename = sys.argv[1]
basedir = os.path.join(basename,"output") 
outfile = os.path.join(basedir,"output.dna")
infile_meta = os.path.join(basedir,"read_meta.txt")
infile_seq = os.path.join(basedir,"read_seq140.txt")
infile_pos = os.path.join(basedir,"read_pos140.txt")
infile_noise = os.path.join(basedir,"read_noise140.txt")
infile_noisepos = os.path.join(basedir,"read_noisepos140.txt")
infile_rev = os.path.join(basedir,"tempRC0.txt")

f_meta = open(infile_meta,'r')
f_pos = open(infile_pos,'r')
f_noise = open(infile_noise,'r')
f_noisepos = open(infile_noisepos,'r')
f_rev = open(infile_rev,'r')
f_out = open(outfile,'w')

readlen = 0 
for line in f_meta:
	current = line.rstrip('\n')
	readlen = len(current)
	break

asciitoint = {'a':0,'b':1,'c':2,'d':3,'e':4,'f':5,'g':6,'h':7,'i':8,'j':9,'k':10,'l':11,'m':12,'n':13,'o':14,'p':15,'q':16,'r':17,'s':18,'t':19,'u':20,'w':21,'x':22,'y':23,'z':24,'A':25,'B':26,'C':27,'D':28,'E':29,'F':30,'G':31,'H':32,'I':33,'J':34,'K':35,'L':36,'M':37,'N':38,'O':39,'P':40,'v':readlen}


with open(infile_seq,'r') as f_seq:
	for line in f_seq:
		ref = line.rstrip('\n')
		pos = f_pos.readline().rstrip('\n')
		prevpos = 0
		for i in pos:
			currentread = ref[prevpos+asciitoint[i]:prevpos+asciitoint[i]+readlen]
			prevpos = prevpos + asciitoint[i]
			noise = f_noise.readline().rstrip('\n')
			prevnoisepos = 0
			for n in noise:
				noisepos = int(f_noisepos.read(2))+prevnoisepos
				currentread = currentread[:noisepos]+n+currentread[noisepos+1:]
				prevnoisepos = noisepos
			rev = f_rev.read(1)
			if rev == 'd':
				f_out.write(currentread+'\n')	
			else:
				f_out.write(str(Seq(currentread).reverse_complement())+'\n')

f_out.close()
f_pos.close()	
f_noise.close()
f_noisepos.close()
f_rev.close()
