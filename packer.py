infile = "temp2.dna"
outfile_seq = "11111.txt"
outfile_pos = "11112.txt"
readlen = 100
minmatch = 20
neg = False #whether we want to match the current read with the previous read shifted rightward 
f_seq = open(outfile_seq,'w')
f_pos = open(outfile_pos,'w')
prev = 'A'*readlen
with open(infile,'r') as f:
	for line in f:
		current = line.rstrip('\n')
		flag = 0
		for i in range(minmatch):
			if(current[:readlen-i]==prev[i:]):
				f_pos.write('+')
				f_seq.write(current[(readlen-i+1):]+'\n')
				prev = current
				flag = 1
				break
		if flag == 0 and neg == True:
			for i in range(minmatch):
				if(current[i:]==prev[:(readlen-i)]):
					f_pos.write('-')
					f_seq.write(current[:i]+'\n')
					prev = current
					flag = 1
					break
		
		if flag == 0:
			f_pos.write('0')
			f_seq.write(current+'\n')
			prev = current
f_seq.close()
f_pos.close()	
		

