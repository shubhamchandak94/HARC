from itertools import imap
import operator

# hamming2 function from http://code.activestate.com/recipes/499304-hamming-distance/
def hamming2(str1, str2):
    #assert len(str1) == len(str2)
    #ne = str.__ne__  ## this is surprisingly slow
    ne = operator.ne
    return sum(imap(ne, str1, str2))
#	return sum([str1[i]!=str2[i] for i in range(len(str1))])	

infile = "chrom22_50x_noRC_noisy.dna"
outfile = "tempte5.dna"
readlen = 100
no_reads = 17500000
matchlen = 80
maxmatch = 20
thresh = [4,5,6,7]+[8]*16 #maximum number of mismatches allowed
ind = [[i for i in range(j,80,5)] for j in range(5)]

print "Reading file"
f = open(infile,'r')
lines = [f.readline().rstrip('\n') for i in range(no_reads)]
f.close()

print "Constructing dictionaries"
d = [{} for i in range(5)]
for i in range(no_reads):
	l = [[lines[i][j] for j in ind[k]] for k in range(5)]
	s = [''.join(l[k]) for k in range(5)]
	for k in range(5):
		if s[k] in d[k]:
			d[k][s[k]].append(i)
		else:
			d[k][s[k]] = [i]

print "Ordering reads and writing to file"
remainingreads = set([i for i in range(no_reads)])
current = 0
fout = open(outfile,'w')
while True:
	flag = 0
	if len(remainingreads)%1000000 == 0:
		print str(len(remainingreads)//1000000)+'M reads remain'
	fout.write(lines[current]+'\n')
	remainingreads.remove(current)
	if(len(remainingreads)==0):
		break
	for i in range(maxmatch):
		l = [[lines[current][j+i] for j in ind[k]] for k in range(5)]
		s = [''.join(l[k]) for k in range(5)]
		inter = set()
		for k in range(5):
			if s[k] in d[k]:
				inter = set(d[k][s[k]]).union(inter)
		inter = inter.intersection(remainingreads)	
		if len(inter)>0:
			for j in inter:
				if(hamming2(lines[current][i:],lines[j][:readlen-i])<=thresh[i]):
					current = j
					flag = 1
					break
		if flag == 1:
			break		

#if len(d[lines[current][0:matchlen]]) == 0:
#		del d[lines[current][0:matchlen]]
#	else:
#		for i in d[lines[current][0:matchlen]]:
#			if hamming2(lines[current][matchlen:],lines[i][matchlen:]) <= thresh:
#				current = i
#				flag = 1
#				break
#	if flag == 1:
#		continue
#	for j in range(1,maxmatch):
#		if lines[current][j:j+matchlen] in d:
#			for i in d[lines[current][j:j+matchlen]]:
#				if hamming2(lines[current][j+matchlen:],lines[i][matchlen:readlen-j]) <= thresh:
#					current = i
#					flag = 1
#					break
#			if flag == 1:
#					break
	if flag == 1:
		continue
	current = remainingreads.pop()
	remainingreads.add(current)

print "Done"
fout.close()