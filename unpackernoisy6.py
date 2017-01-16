#find best match with ref or prev - for files with no Ns

outfile = "temptemp.dna"
infile_seq = "read_seq107.txt"
infile_pos = "read_pos107.txt"
infile_noise = "read_noise107.txt"
infile_noisepos = "read_noisepos107.txt"

readlen = 100
asciitoint = {'a':0,'b':1,'c':2,'d':3,'e':4,'f':5,'g':6,'h':7,'i':8,'j':9,'k':10,'l':11,'m':12,'n':13,'o':14,'p':15,'q':16,'r':17,'s':18,'t':19,'u':20,'v':readlen}
inttoascii = {0:'a',1:'b',2:'c',3:'d',4:'e',5:'f',6:'g',7:'h',8:'i',9:'j',10:'k',11:'l',12:'m',13:'n',14:'o',15:'p',16:'q',17:'r',18:'s',19:'t',20:'u',readlen:'v'}


f_pos = open(infile_pos,'r')
f_noise = open(infile_noise,'r')
f_noisepos = open(infile_noisepos,'r')
f_out = open(outfile,'w')

with open(infile_seq,'r') as f_seq:
	for line in f_seq:
		ref = line.rstrip('\n')
		pos = f_pos.readline().rstrip('\n')
		prevpos = 0
		for i in pos:
			read = ref[prevpos+asciitoint[i]:prevpos+asciitoint[i]+readlen]
			prevpos = prevpos + asciitoint[i]
			noise = f_noise.readline().rstrip('\n')
			prevnoisepos = 0
			for n in noise:
				noisepos = int(f_noisepos.read(2))+prevnoisepos
				read = read[:noisepos]+n+read[noisepos+1:]
				prevnoisepos = noisepos
			f_out.write(read+'\n')	

f_out.close()
f_pos.close()	
f_noise.close()
f_noisepos.close()

