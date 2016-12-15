from distance import hamming
from Bio.Seq import Seq
from bitstring import BitArray

infile = "SRR065390.dna"
outfile = "temp.dna"
readlen = 100
no_reads = 67617092
matchlen = 80 
maxmatch = 20
num_dict = 1 # should divide matchlen

#ind = [[i for i in range(j,matchlen,num_dict)] for j in range(num_dict)]#
ind = [[i for i in range(40,60)]]

char2base = {'A':[0,0],'C':[1,0],'G':[0,1],'T':[1,1]}

def convert2bitarray(l):
	return BitArray([i for j in l for i in char2base[j]])

listN = []
lines = []
print "Reading file"
f = open(infile,'r')
for i in range(no_reads):
	line = f.readline().rstrip('\n')
	if 'N' in line:
		listN.append(line)
	else:
		lines.append(convert2bitarray(line))
	
f.close()

no_cleanreads = len(lines)

print "Generating reverse complements"
revlines = [~lines[i] for i in range(no_cleanreads)]

print "Constructing dictionaries"
d1 = [{} for i in range(num_dict)]
d2 = [{} for i in range(num_dict)] #reverse complement dictionary
for i in range(no_reads):
	l = [[lines[i][j] for j in ind[k]] for k in range(num_dict)]
	s = [''.join(l[k]) for k in range(num_dict)]
	for k in range(num_dict):
		if s[k] in d1[k]:
			d1[k][s[k]].append(i)
		else:
			d1[k][s[k]] = [i]

	l = [[revlines[i][j] for j in ind[k]] for k in range(num_dict)]
	s = [''.join(l[k]) for k in range(num_dict)]
	for k in range(num_dict):
		if s[k] in d2[k]:
			d2[k][s[k]].append(i)
		else:
			d2[k][s[k]] = [i]

print "Ordering reads and writing to file"
f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
f_noise = open(outfile_noise,'w')
f_noisepos = open(outfile_noisepos,'w')
f_rev = open(outfile_rev,'w')

remainingreads = set([i for i in range(no_reads)])
current = 0
prev = lines[current]
count = [[0]*readlen for j in range(5)]
for j in range(readlen):
	count[char2index(lines[current][j])][j] = 1
f_flag.write('0')
f_seq.write(lines[current]+'\n')

ref = findmajority(count)

while True:
	flag = 0
	if len(remainingreads)%1000000 == 0:
		print str(len(remainingreads)//1000000)+'M reads remain'
	remainingreads.remove(current)
	if(len(remainingreads)==0):
		break
	for i in range(maxmatch):
		l = [[ref[j+i] for j in ind[k]] for k in range(num_dict)]
		s = [''.join(l[k]) for k in range(num_dict)]
		inter = set()
		for k in range(num_dict):	
			if s[k] in d1[k]:
				inter = inter.union(set(d1[k][s[k]]))
		inter = inter.intersection(remainingreads)
		if len(inter)>0:
			for j in inter:
			#	if(hamming(ref[i:],lines[j][:readlen-i])<=thresh):
					f_rev.write('d') # d for direct
					current = j
					currentseq = lines[j]
					flag = 1
					break

		if flag == 0:
			inter = set()
			for k in range(num_dict):	
				if s[k] in d2[k]:
					inter = inter.union(set(d2[k][s[k]]))
			inter = inter.intersection(remainingreads)
			if len(inter)>0:
				for j in inter:
				#	if(hamming(ref[i:],revlines[j][:readlen-i])<=thresh):
						f_rev.write('r') # r for reverse
						current = j
						currentseq = revlines[j]
						flag = 1
						break
		if flag == 1:
			if(hamming(currentseq[:(readlen-i)],ref[i:])<=hamming(currentseq[:(readlen-i)],prev[i:])):
				f_flag.write('r') #r for ref
				f_seq.write(currentseq[(readlen-i):]+'\n')
				prevj = 0;
				for j in range(readlen-i):
					count[char2index(currentseq[j])][i+j] += 1		
					if currentseq[j]!=ref[i+j]:
						f_noise.write(currentseq[j])
						f_noisepos.write("%02d"%(j-prevj))#delta encoding
						prevj = j	
			else:
				f_flag.write('p') #p for prev
				f_seq.write(currentseq[(readlen-i):]+'\n')
				prevj = 0;
				for j in range(readlen-i):
					count[char2index(currentseq[j])][i+j] += 1		
					if currentseq[j]!=prev[i+j]:
						f_noise.write(currentseq[j])
						f_noisepos.write("%02d"%(j-prevj))#delta encoding
						prevj = j	
			
			count = [count[j][i:]+[0]*i for j in range(5)]
			for j in range(readlen-i,readlen):
				count[char2index(currentseq[j])][j] = 1
			
			ref = findmajority(count)	
			f_noise.write('\n')
			break

	if flag == 1:
		prev = currentseq
		continue
	
	current = remainingreads.pop()
	remainingreads.add(current)
	prev = lines[current]
	count = [[0]*readlen for j in range(5)]
	for j in range(readlen):
		count[char2index(lines[current][j])][j] = 1
	f_flag.write('0')
	f_seq.write(lines[current]+'\n')

	ref = findmajority(count)

print "Done"
f_seq.close()
f_flag.close()	
f_rev.close()
f_noise.close()
f_noisepos.close()

