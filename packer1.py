
#Instead of \n in read_seq we have a separate file readpos which stores the offset from previous read
infile = "temp2.dna"
outfile_seq = "readseq.txt"
outfile_match = "readmatch.txt"
outfile_pos = "readpos.txt"

readlen = 100
minmatch = 20
neg = False #whether we want to match the current read with the previous read shifted rightward 
f_seq = open(outfile_seq,'w')
f_match = open(outfile_match,'w')
f_pos = open(outfile_pos,'w')
prev = 'A'*readlen
with open(infile,'r') as f:
	for line in f:
		current = line.rstrip('\n')
		flag = 0
		for i in range(minmatch):
			if(current[:readlen-i]==prev[i:]):
				f_match.write('+')
				f_seq.write(current[(readlen-i+1):])
				f_pos.write("%02d"%i)
				prev = current
				flag = 1
				break
		if flag == 0 and neg == True:
			for i in range(minmatch):
				if(current[i:]==prev[:(readlen-i)]):
					f_match.write('-')
					f_seq.write(current[:i])
					f_pos.write("%02d"%i)
					prev = current
					flag = 1
					break
		
		if flag == 0:
			f_match.write('0')
			f_seq.write(current)
			prev = current
f_seq.close()
f_pos.close()	
f_match.close()		
