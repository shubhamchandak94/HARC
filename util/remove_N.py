import sys
import os

inFile = sys.argv[1]
print("infile: ", inFile)

outFile = os.path.splitext(inFile)[0] + "_clean" + os.path.splitext(inFile)[1]
outFile_N = os.path.splitext(inFile)[0] + "_N" + os.path.splitext(inFile)[1]

print("outFile clean: ", outFile)
print("outfile N:", outFile_N)
fout = open(outFile,'w')
fout_N = open(outFile_N,'w')
with open(inFile,'r') as f:
  for line in f:
     if 'N' not in line:
          fout.write(line)
     else:
	  fout_N.write(line)
