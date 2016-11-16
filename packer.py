infile = "temp2.dna"
outfile_seq = "read_seq4.txt"
outfile_flag = "read_flag4.txt"
readlen = 98
maxmatch = 18
neg = False #whether we want to match the current read with the previous read shifted rightward 
f_seq = open(outfile_seq,'w')
f_flag = open(outfile_flag,'w')
prev = 'A'*readlen
k = 0
with open(infile,'r') as f:
	for line in f:
		k = k + 1
		if (k%1000000 == 0):
			print str(k//1000000)+"M reads done"
		current = line.rstrip('\n')
		flag = 0
		for i in range(maxmatch):
			if(current[:readlen-i]==prev[i:]):
				f_flag.write('+')
				f_seq.write(current[(readlen-i):]+'\n')
				prev = current
				flag = 1
				break
		if flag == 0 and neg == True:
			for i in range(maxmatch):
				if(current[i:]==prev[:(readlen-i)]):
					f_flag.write('-')
					f_seq.write(current[:i]+'\n')
					prev = current
					flag = 1
					break
		
		if flag == 0:
			f_flag.write('0')
			f_seq.write(current+'\n')
			prev = current
f_seq.close()
f_flag.close()			
