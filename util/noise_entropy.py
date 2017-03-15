#!/usr/bin/python2
## computes noise_entropy by assuming per position

import sys
import os
import numpy as np
from tqdm import tqdm


###################################
def quality_to_prob(qual_string):
	_q = [ord(c) for c in qual_string]

	_q = (33.0 - np.array(_q))/10.0
	prob = np.power(10.0,_q)
	return prob


###################################
def get_M_N(inFile):
	M = 0
	N = 0
	with open(inFile,'r') as f:
	  for line in f:
	  	M += 1
	  	N = len(line) -1
	return (M,N)

###################################
def compute_entropy(probs):
	_e = -( probs*np.log(probs) + (1-probs)*np.log(1-probs) + probs*np.log(3.0) ) 
	_e = _e/np.log(2.0)
	entropy = np.sum(_e)
	return entropy,_e 

###################################
def main():
	inFile = sys.argv[1]
	print("infile: ", inFile)
	M,N = get_M_N(inFile)
	print "Read Length: ", N
	probs_sum = np.zeros(N)
	with open(inFile,'r') as f:
		for line in tqdm(f,total=M,desc="computing read probs ...",ascii=True):
			line = line.rstrip('\n')
			_p = quality_to_prob(line)
			probs_sum += _p
			# print line
			# print _p

	probs = probs_sum/(1.0*M)
	entropy,_e = compute_entropy(probs)

	print "Error Probability: "
	print probs
	print "Entropy per position: "
	print _e

	print "Overall Entropy:", entropy

###################################
if __name__ == '__main__':
    main()

