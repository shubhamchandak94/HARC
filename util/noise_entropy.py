import sys
import os
import numpy as np
import struct

readlen = 63


###################################
def quality_to_prob(qual_string):
	_q = [ord(c) for c in qual_string]

	_q = (33.0 - np.array(_q))/10.0
	prob = np.power(10.0,_q)
	return prob

###################################
def get_M_readlen(inFile):
	M = 0
	N = 0
	with open(inFile,'r') as f:
	  for line in f:
	  	M += 1
	  	N = len(line) -1
	return (M,N)


###################################
def qv_to_prob():
	_q  = -np.arange(42)/10.0 
	prob = np.reshape(np.power(10.0,_q),[1,42])
	return prob


def quality_value_stats(inFile):
  
    qv_counts_0order = np.ones((42,readlen))
    qv_counts_1order = np.ones((42,42,readlen))
    qv_counts_2order = np.ones((42,42,42,readlen))
    
    qv_counts_0order_1 = np.ones((42,readlen))
    qv_counts_1order_1 = np.ones((42,42,readlen))
    qv_counts_2order_1 = np.ones((42,42,42,readlen))
    
    qv_counts_0order_2 = np.ones((42,readlen))
    qv_counts_1order_2 = np.ones((42,42,readlen))
    qv_counts_2order_2 = np.ones((42,42,42,readlen))

    
    f = open(inFile,'r')
    s = f.read(8)
    num_reads_1 = struct.unpack('Q',s)[0]
    for i in range(readlen):
	for j in range(42):
		s = f.read(8)
		qv_counts_0order_1[j,i] = struct.unpack('Q',s)[0] 
    for i in range(readlen-1):
	for j in range(42):
		for k in range(42):
			s = f.read(8)
			qv_counts_1order_1[j,k,i] = struct.unpack('Q',s)[0] 
    for i in range(readlen-2):
	for j in range(42):
		for k in range(42):
			for l in range(42):
				s = f.read(8)
				qv_counts_2order_1[j,k,l,i] = struct.unpack('Q',s)[0] 

    s = f.read(8)
    num_reads_2 = struct.unpack('Q',s)[0]
    for i in range(readlen):
	for j in range(42):
		s = f.read(8)
		qv_counts_0order_2[j,i] = struct.unpack('Q',s)[0]
    for i in range(readlen-1):
	for j in range(42):
		for k in range(42):
			s = f.read(8)
			qv_counts_1order_2[j,k,i] = struct.unpack('Q',s)[0] 
    for i in range(readlen-2):
	for j in range(42):
		for k in range(42):
			for l in range(42):
				s = f.read(8)
				qv_counts_2order_2[j,k,l,i] = struct.unpack('Q',s)[0] 
    
    num_reads = num_reads_1 + num_reads_2	
    qv_prob_0order = (qv_counts_0order_1+qv_counts_0order_2)/(1.0*num_reads)
    qv_prob_1order = (qv_counts_1order_1+qv_counts_1order_2)/(1.0*num_reads)
    qv_prob_2order = (qv_counts_2order_1+qv_counts_2order_2)/(1.0*num_reads)

    qv_prob_0order_1 = qv_counts_0order_1/(1.0*num_reads_1)
    qv_prob_1order_1 = qv_counts_1order_1/(1.0*num_reads_1)
    qv_prob_2order_1 = qv_counts_2order_1/(1.0*num_reads_1)
    
    qv_prob_0order_2 = qv_counts_0order_2/(1.0*num_reads_2)
    qv_prob_1order_2 = qv_counts_1order_2/(1.0*num_reads_2)
    qv_prob_2order_2 = qv_counts_2order_2/(1.0*num_reads_2)
    
    return (qv_prob_0order, qv_prob_1order, qv_prob_2order,
            qv_prob_0order_1, qv_prob_1order_1, qv_prob_2order_1,
            qv_prob_0order_2, qv_prob_1order_2, qv_prob_2order_2,
            num_reads,num_reads_1, num_reads_2)


def avg_quality_value(inFile):
    num_reads,readlen = get_M_readlen(inFile)
    qv_avg = np.zeros(num_reads)
    with open(inFile,'r') as f_qv:
        for i,qv in tqdm(enumerate(f_qv),total=num_reads,desc="computing read probs ...",ascii=True):
            qv_id = qv.rstrip('\n')
            qv_id = [ord(c) - 33.0 for c in qv_id]
            qv_avg[i] = np.mean(qv_id)
    return qv_avg

inFile = "../temp"
ret = quality_value_stats(inFile)


# In[ ]:

(qv_prob_0order, qv_prob_1order, qv_prob_2order,qv_prob_0order_1, qv_prob_1order_1, qv_prob_2order_1,qv_prob_0order_2, qv_prob_1order_2, qv_prob_2order_2,num_reads,num_reads_1, num_reads_2) = ret


# In[ ]:

# print qv_prob
# print "Cluster 1 \n", qv_prob_1
# print "Cluster 2 \n", qv_prob_2

print "Total number of reads"
print num_reads
print "Number of reads in cluster 1"
print num_reads_1
print "Number of reads in cluster 2"
print num_reads_2


# In[ ]:

###################################
def compute_Ni_old_entropy(probs):
	_e = -( probs*np.log(probs) + (1-probs)*np.log(1-probs) ) + probs*np.log(3.0)  
	_e = _e/np.log(2.0)
	entropy = np.sum(_e)
	return entropy,_e 

###################################
def compute_Ni_entropy(probs):
	_e = -( probs*np.log(probs) + (1-probs)*np.log(1-probs)  + probs*(0.3*np.log(0.15)+0.7*np.log(0.7)))
	_e = _e/np.log(2.0)
	entropy = np.sum(_e)
	return entropy,_e 

def compute_Si_prob(qv_prob):
    prob_substitution = np.dot(qv_to_prob(),qv_prob)
    return prob_substitution


# In[ ]:

Si_prob = compute_Si_prob(qv_prob_0order)
# print Si_prob 

Si_prob_1 = compute_Si_prob(qv_prob_0order_1)
# print Si_prob_1 

Si_prob_2 = compute_Si_prob(qv_prob_0order_2)
# print Si_prob_2 


def xlogx(p):
    _e = p*np.log(p)
    _e = _e/np.log(2.0)
    return _e


# In[ ]:


# In[ ]:


def compute_N1N2_old_entropy(qv_N1N2_probs):
    qv_to_prob_0 = 1.0 - qv_to_prob()
    qv_to_prob_1 = qv_to_prob()
    
    temp_prob_00 = np.dot(np.dot(qv_to_prob_0,qv_N1N2_probs),np.transpose(qv_to_prob_0))
    temp_prob_01 = np.dot(np.dot(qv_to_prob_0,qv_N1N2_probs),np.transpose(qv_to_prob_1))/3.0
    temp_prob_10 = np.dot(np.dot(qv_to_prob_1,qv_N1N2_probs),np.transpose(qv_to_prob_0))/3.0
    temp_prob_11 = np.dot(np.dot(qv_to_prob_1,qv_N1N2_probs),np.transpose(qv_to_prob_1))/9.0
    entropy = -(xlogx(temp_prob_00) + 3.0*xlogx(temp_prob_01) + 3.0*xlogx(temp_prob_10) + 9.0*xlogx(temp_prob_11))
    #print entropy
    return entropy

def compute_N1N2N3_entropy(qv_N1N2N3_probs):
    qv_to_prob_0 = 1.0 - qv_to_prob()
    qv_to_prob_1 = qv_to_prob()
    qv_2_prob = np.stack((qv_to_prob_0, qv_to_prob_1))
    #print qv_2_prob.shape
    noise_entropy = -(0.3*np.log(0.15)+0.7*np.log(0.7))/np.log(2.0)

    entropy = 0
    
    def prob_S1S2S3(s1,s2,s3):
        temp_prob = 0;
        for i0 in range(42):
            for i1 in range(42):
                for i2 in range(42):
                    temp_prob += qv_N1N2N3_probs[i0,i1,i2]*qv_2_prob[s1,0,i0]*qv_2_prob[s2,0,i1]*qv_2_prob[s3,0,i2]
        return temp_prob
    
    prob_s1s2s3 = np.zeros((2,2,2))
    for s1 in [0,1]:
        for s2 in [0,1]:
            for s3 in [0,1]:
                prob_s1s2s3[s1,s2,s3] = prob_S1S2S3(s1,s2,s3)
                entropy += (-xlogx(prob_s1s2s3[s1,s2,s3]))
    prob_s1 = np.sum(prob_s1s2s3[1,:,:])
    prob_s2 = np.sum(prob_s1s2s3[:,1,:])
    prob_s3 = np.sum(prob_s1s2s3[:,:,1])
    #print np.sum(prob_s1s2s3)
    entropy += (prob_s1 + prob_s2 + prob_s3)*noise_entropy
    #print entropy
    
    return entropy
    
                

def compute_N1N2_entropy(qv_N1N2_probs):
    qv_to_prob_0 = 1.0 - qv_to_prob()
    qv_to_prob_1 = qv_to_prob()
    noise_entropy = -(0.3*np.log(0.15)+0.7*np.log(0.7))/np.log(2.0)
    
    temp_prob_00 = np.dot(np.dot(qv_to_prob_0,qv_N1N2_probs),np.transpose(qv_to_prob_0))
    temp_prob_01 = np.dot(np.dot(qv_to_prob_0,qv_N1N2_probs),np.transpose(qv_to_prob_1))
    temp_prob_10 = np.dot(np.dot(qv_to_prob_1,qv_N1N2_probs),np.transpose(qv_to_prob_0))
    temp_prob_11 = np.dot(np.dot(qv_to_prob_1,qv_N1N2_probs),np.transpose(qv_to_prob_1))
    p_S1_1 = temp_prob_10 + temp_prob_11
    p_S2_1 = temp_prob_01 + temp_prob_11
    entropy = -(xlogx(temp_prob_00) + xlogx(temp_prob_01) + xlogx(temp_prob_10) + xlogx(temp_prob_11))
    entropy += (p_S1_1 + p_S2_1)*noise_entropy
    #print entropy
    return entropy
    
def compute_Ni_joint_probability(qv_joint_probs,Ni_perpos_entropy,total_Ni_entropy):
    readlen = Ni_perpos_entropy.shape[1]
    #print readlen
    #print total_Ni_entropy
    entropy = Ni_perpos_entropy[0,0] 
    #print Ni_perpos_entropy
    perpos_1order_entropy = np.zeros(readlen-1)
    for i in range(readlen-1):
        qv_N1N2_probs = qv_joint_probs[:,:,i]
        perpos_1order_entropy[i] = compute_N1N2_entropy(qv_N1N2_probs)
        entropy += perpos_1order_entropy[i] - Ni_perpos_entropy[0,i]
    return entropy,perpos_1order_entropy

def compute_2order_entropy(qv_prob_2order, perpos_1order_entropy, total_1order_entropy):
    readlen = perpos_1order_entropy.shape[0] + 1
    #print readlen
    #print total_Ni_entropy
    entropy = perpos_1order_entropy[0]
    #print "initiali entropy:",entropy
    #print Ni_perpos_entropy
    
    for i in range(readlen-2):
        qv_N1N2N3_probs = qv_prob_2order[:,:,:,i]
        entropy += compute_N1N2N3_entropy(qv_N1N2N3_probs) - perpos_1order_entropy[i]
        #print "cumul: ", entropy
    return entropy


print "\n0th order"
total_Ni_entropy,Ni_perpos_entropy = compute_Ni_entropy(Si_prob)
# print Ni_perpos_entropy
print "Overall (without clustering):",total_Ni_entropy

total_Ni_entropy_1,Ni_perpos_entropy_1 = compute_Ni_entropy(Si_prob_1)
# print Ni_perpos_entropy_1
print "Cluster 1:",total_Ni_entropy_1

total_Ni_entropy_2,Ni_perpos_entropy_2 = compute_Ni_entropy(Si_prob_2)
# print Ni_perpos_entropy_2
print "Cluster 2: ",total_Ni_entropy_2


print "Overall (with clustering): ", -xlogx(1.0*num_reads_1/num_reads) - xlogx(1.0*num_reads_2/num_reads) + total_Ni_entropy_1*num_reads_1/num_reads +            total_Ni_entropy_2*num_reads_2/num_reads


# In[ ]:
print "\n1st order entropy"
total_1order_entropy, perpos_1order_entropy = compute_Ni_joint_probability(qv_prob_1order,Ni_perpos_entropy,total_Ni_entropy)
print "Overall (without clustering): ",total_1order_entropy

total_1order_entropy_1, perpos_1order_entropy_1 = compute_Ni_joint_probability(qv_prob_1order_1,Ni_perpos_entropy_1,total_Ni_entropy_1)
print "Cluster 1: ",total_1order_entropy_1

total_1order_entropy_2, perpos_1order_entropy_2 = compute_Ni_joint_probability(qv_prob_1order_2,Ni_perpos_entropy_2,total_Ni_entropy_2)
print "Cluster 2: ",total_1order_entropy_2

print "Overall (with clustering): ", -xlogx(1.0*num_reads_1/num_reads) - xlogx(1.0*num_reads_2/num_reads) + total_1order_entropy_1*num_reads_1/num_reads +            total_1order_entropy_2*num_reads_2/num_reads


# In[ ]:
print "\n2nd order entropy"
entropy_2order = compute_2order_entropy(qv_prob_2order, perpos_1order_entropy, total_1order_entropy)
print "Overall (without clustering): ", entropy_2order
entropy_2order_1 = compute_2order_entropy(qv_prob_2order_1, perpos_1order_entropy_1, total_1order_entropy_1)
print "Cluster 1: ",entropy_2order_1
entropy_2order_2 = compute_2order_entropy(qv_prob_2order_2, perpos_1order_entropy_2, total_1order_entropy_2)
print "Cluster 2: ", entropy_2order_2
print "Overall (with clustering): ",-xlogx(1.0*num_reads_1/num_reads) - xlogx(1.0*num_reads_2/num_reads) + entropy_2order_1*num_reads_1/num_reads +            entropy_2order_2*num_reads_2/num_reads

