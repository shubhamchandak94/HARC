import operator
from itertools import imap

#hamming2 function from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming2(str1, str2):
    #assert len(str1) == len(str2)
    #ne = str.__ne__  ## this is surprisingly slow
    ne = operator.ne
    return sum(imap(ne, str1, str2))

def char2index(c):
	if c == 'A':
		return 0		
	if c == 'C':
		return 1
	if c == 'G':
		return 2
	if c == 'T':
		return 3
	if c == 'N':
		return 4
def findmajority(count):
	l = []
	for i in range(len(count[0])):
		s = [count[j][i] for j in range(5)]
		maxcount = max(s[0:4])
		if maxcount == 0: #only N's seen so far
			l.append('A')
			continue
		if s[0] == maxcount:
			l.append('A')
			continue
		if s[1] == maxcount:
			l.append('C')
			continue
		if s[2] == maxcount:
			l.append('G')
			continue
		if s[3] == maxcount:
			l.append('T')
			continue
	return ''.join(l)



infile = "tempte1.dna"
outfile_seq = "read_seq22.txt"
outfile_pos = "read_pos22.txt"
outfile_noise = "read_noise22.txt"
outfile_noisepos = "read_noisepos22.txt"

readlen = 100
minmatch = 20
thresh = 20 # maximum number of mismatches allowed 

f_seq = open(outfile_seq,'w')
f_pos = open(outfile_pos,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
k = 0
ref = 'A'*readlen # ref is the reference which is constantly updated (introduced because matching a read to previous read leads to double noise than actual)
prev = ['A'*readlen]*8
count = [[1]*100,[0]*100,[0]*100,[0]*100,[0]*100] #number of A's,C's,T's,G's and N's seen at each position in ref
#Note: N is never considered in the ref - we arbitrarily place an A if only N's are seen at some position
with open(infile,'r') as f:
	for line in f:
		k = k + 1
		if k%1000000 == 0:
			print str(k//1000000)+'M done'
		current = line.rstrip('\n')
		flag = 0
		for i in range(minmatch):
			if(hamming2(current[:(readlen-i)],ref[i:])<=thresh):
				if i == 0:
					bestmatch = hamming2(current,ref)
					bestmatchpos = -1
					for j in range(8):
						if(hamming2(current,prev[j])<bestmatch):
							bestmatch = hamming2(current,prev[j])
							bestmatchpos = j
					if bestmatchpos == -1:
						f_pos.write('r')
						f_seq.write('\n')
						prevj = 0;
						for j in range(readlen):
							count[char2index(current[j])][j] += 1		
							if current[j]!=ref[j]:
								f_noise.write(current[j])
								f_noisepos.write("%02d"%(j-prevj))#delta encoding
								prevj = j	
					else:
						f_pos.write(str(bestmatchpos+1))
						f_seq.write('\n')
						prevj = 0;
						for j in range(readlen):
							count[char2index(current[j])][j] += 1		
							if current[j]!=prev[bestmatchpos][j]:
								f_noise.write(current[j])
								f_noisepos.write("%02d"%(j-prevj))#delta encoding
								prevj = j	
				else:			
					if(hamming2(current[:(readlen-i)],ref[i:])<=hamming2(current[:(readlen-i)],prev[0][i:])):
						f_pos.write('r')
						f_seq.write(current[(readlen-i+1):]+'\n')
						prevj = 0;
						for j in range(readlen-i):
							count[char2index(current[j])][i+j] += 1		
							if current[j]!=ref[i+j]:
								f_noise.write(current[j])
								f_noisepos.write("%02d"%(j-prevj))#delta encoding
								prevj = j	
					else:
						f_pos.write('1')
						f_seq.write(current[(readlen-i+1):]+'\n')
						prevj = 0;
						for j in range(readlen-i):
							count[char2index(current[j])][i+j] += 1		
							if current[j]!=prev[0][i+j]:
								f_noise.write(current[j])
								f_noisepos.write("%02d"%(j-prevj))#delta encoding
								prevj = j	
					
				count = [count[j][i:]+[0]*i for j in range(5)]
				for j in range(readlen-i,readlen):
					count[char2index(current[j])][j] = 1
				
				ref = findmajority(count)	
				#ref = current#ref[i:]+current[readlen-i:]
				f_noise.write('\n')
				flag = 1
				break
		
		if flag == 0:
			f_pos.write('0')
			f_seq.write(current+'\n')
			count = [[0]*100 for j in range(5)]
			for j in range(readlen):
				count[char2index(current[j])][j] = 1
			ref = findmajority(count)
		prev = [current]+prev[:7]						
f_seq.close()
f_pos.close()	
f_noise.close()
f_noisepos.close()

