
# coding: utf-8

# In[4]:

#!/usr/bin/python
import sys
import numpy as np
import editdistance
import distance
from itertools import izip
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from subprocess import call


# ## Files to be read

prefix = sys.argv[1]
input_fastq = prefix + ".fq"
clean_file = prefix + ".clean"
denoised_file = prefix + ".fq.dna.d"
noisy_file = prefix + ".noisy"
image_prefix = sys.argv[2] 


def get_M_readlen(inFile):
    M = 0
    N = 0
    with open(inFile,'r') as f:
      for line in f:
        M += 1
        N = len(line) -1
    return (M,N)



M,readlen = get_M_readlen(clean_file)

per_position_original_errors = np.zeros(readlen+1)
per_position_denoised_errors = np.zeros(readlen+1)
per_position_undetected_errors = np.zeros(readlen+1)
per_position_new_errors = np.zeros(readlen+1)
new_errors_hist = np.zeros(readlen+1)

for clean_read, denoised_read, noisy_read in tqdm( izip(open(clean_file), open(denoised_file), open(noisy_file)), total=M,desc="computing distance ...", ascii=True):

    true_errors = (np.fromstring(clean_read,dtype=np.uint8)!=np.fromstring(noisy_read,dtype=np.uint8))
    after_denoising_errors = (np.fromstring(clean_read,dtype=np.uint8)!=np.fromstring(denoised_read,dtype=np.uint8))
    undetected_errors = np.logical_and(true_errors,after_denoising_errors)
    new_errors = np.logical_and(np.logical_not(true_errors), after_denoising_errors)
    
    num_new_errors = np.sum(new_errors)
    new_errors_hist[num_new_errors] += 1
    per_position_original_errors += true_errors
    per_position_denoised_errors += after_denoising_errors
    per_position_undetected_errors += undetected_errors
    per_position_new_errors += new_errors

total_errors = np.sum(per_position_original_errors) 
total_after_denoising_errors = np.sum(per_position_denoised_errors) 
total_undetected_errors = np.sum(per_position_undetected_errors) 
total_new_errors = np.sum(per_position_new_errors)

print sys.argv[2]
print "Original Errors: ", total_errors
print "After Denoising Errors: ", total_after_denoising_errors
print "After Undetected Errors: ", total_undetected_errors
print "False Positives: ", total_new_errors
print "\n"

plt.plot(per_position_original_errors)
plt.ylabel('per_position_errors')
plt.savefig(image_prefix + ".per_pos_errors.png")

plt.plot(per_position_denoised_errors)
plt.ylabel('per_position_after_denoising_errors')
plt.savefig(image_prefix + ".per_pos_efter_denoising_errors.png")

plt.plot(per_position_undetected_errors)
plt.ylabel('per_position_undetected_errors')
plt.savefig(image_prefix + ".per_pos_undetected_errors.png")

plt.plot(per_position_new_errors)
plt.ylabel('new_errors')
plt.savefig(image_prefix + ".new_errors.png")


plt.plot(np.log10(new_errors_hist))
plt.ylabel('new_errors histogram')
plt.savefig(image_prefix + ".new_errors_histogram.png")
    
