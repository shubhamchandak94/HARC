#find best match with ref or prev
from joblib import Parallel, delayed
from distance import hamming

infile = "temp.dna"
outfile_seq = "read_seq31.txt"
outfile_flag = "read_flag31.txt"
outfile_noise = "read_noise31.txt"
outfile_noisepos = "read_noisepos31.txt"

no_reads = 5372832
readlen = 98
maxmatch = 18
thresh = 20 # maximum number of mismatches allowed 
numcores = 12

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

def pack(startline,endline,s_seq,s_flag,s_noise,s_noisepos):
  f = open(infile,'r')
  f.seek(startline*(readlen+1))
  ref = 'A'*readlen # ref is the reference which is constantly updated (introduced because matching a read to previous read leads to double noise than actual)
  prev = 'A'*readlen
  count = [[1]*readlen,[0]*readlen,[0]*readlen,[0]*readlen,[0]*readlen] #number of A's,C's,T's,G's and N's seen at each position in ref
  #Note: N is never considered in the ref - we arbitrarily place an A if only N's are seen at some position
  for k in range(startline,endline+1):
      current = f.readline().rstrip('\n')
      flag = 0
      for i in range(maxmatch):
        if(hamming(current[:(readlen-i)],ref[i:])<=thresh):
          if(hamming(current[:(readlen-i)],ref[i:])<=hamming(current[:(readlen-i)],prev[i:])):
            s_flag.append('r')
            s_seq.append(current[(readlen-i):]+'\n')
            prevj = 0;
            for j in range(readlen-i):
              count[char2index(current[j])][i+j] += 1		
              if current[j]!=ref[i+j]:
                s_noise.append(current[j])
                s_noisepos.append("%02d"%(j-prevj))#delta encoding
                prevj = j	
          else:
            s_flag.append('p')
            s_seq.append(current[(readlen-i):]+'\n')
            prevj = 0;
            for j in range(readlen-i):
              count[char2index(current[j])][i+j] += 1		
              if current[j]!=prev[i+j]:
                s_noise.append(current[j])
                s_noisepos.append("%02d"%(j-prevj))#delta encoding
                prevj = j	

          count = [count[j][i:]+[0]*i for j in range(5)]
          for j in range(readlen-i,readlen):
            count[char2index(current[j])][j] = 1

          ref = findmajority(count)	
          #ref = current#ref[i:]+current[readlen-i:]
          s_noise.append('\n')
          flag = 1
          break

      if flag == 0:
        s_flag.append('0')
        s_seq.append(current+'\n')
        count = [[0]*readlen for j in range(5)]
        for j in range(readlen):
          count[char2index(current[j])][j] = 1
        ref = findmajority(count)
      prev = current						

f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')

s_seq = ['']*numcores
s_flag = ['']*numcores
s_noise = ['']*numcores
s_noisepos = ['']*numcores

i = no_reads/numcores
startline = [i*j for j in range(numcores)]
endline = [startline[j+1]-1 for j in range(numcores-1)]
endline.append(no_reads)

Parallel(n_jobs=numcores)(delayed(pack)(startline[i],endline[i],s_seq[i],s_flag[i],s_noise[i],s_noisepos[i]) for i in range(numdict))

for i in range(numcores):
  f_seq.write(s_seq[i]+'\n')
  f_flag.write(s_flag[i]+'\n')
  f_noise.write(s_noise[i]+'\n')
  f_noisepos.write(s_noisepos[i]+'\n')

f_seq.close()
f_flag.close()	
f_noise.close()
f_noisepos.close()
