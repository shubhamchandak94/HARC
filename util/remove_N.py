import sys
import os

inFile = sys.argv[1]
print("infile: ", inFile)

outFile = os.path.splitext(inFile)[0] + "_clean" + os.path.splitext(inFile)[1]
print("outFile: ", outFile)
fout = open(outFile,'w')
with open(inFile,'r') as f:
  for line in f:
     if 'N' not in line:
          fout.write(line)