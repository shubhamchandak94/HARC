def DUDEBEC(Z,k):

    n = len(Z)
    Y = [0]*n

    #first and last k symbols
    for i in range(k):
      if Z[i] == 2:
        Y[i] = 0
      else: 
        Y[i] = Z[i]
        

    for i in range(n-k,n):
      if Z[i] == 2:
        Y[i] = 0
      else: 
        Y[i] = Z[i]
        

    #generating counts
    d = {}
    for i in range(k,n-k):
        l = Z[i-k:i-1]+Z[i+1:i+k]
        s = ''.join([str(j) for j in l])
        if s in d:
            if Z[i] == 0:
                d[s][0]+=1
            elif Z[i] == 1:
                d[s][1]+=1
        else:
            if Z[i] == 0:
                d[s] = [1,0]
            elif Z[i] == 1:
                d[s] = [0,1]
            else:
                d[s] = [0,0]

    for i in range(k,n-k):
        if Z[i] == 2:
            l = Z[i-k:i-1]+Z[i+1:i+k]
            s = ''.join([str(j) for j in l])
            if d[s][1] > d[s][0]:
                Y[i] = 1
            else:
                Y[i] = 0
        else:
            Y[i] = Z[i]
    return Y       

def BER(X,Y):
    n = len(X)
    return sum([X[i]!=Y[i] for i in range(n)])/n
            
import math
import random

n = 1000000
alpha = 0.1
e = 0.2

X = [0]*n;
Z = [0]*n;
k = math.floor(0.5*math.log(n)/math.log(3));
D = alpha*e/(1-(1-2*alpha)*e*e);

#Generate X
X[0] = round(random.random());
for i in range(1,n):
  r = random.random();
  if r < alpha:
    X[i] = 1 - X[i-1];
  else:
    X[i] = X[i-1];

#Generate Z
for i in range(n):
  r = random.random()
  if r < e:
    Z[i] = 2;
  else:
      Z[i] = X[i];

Y = DUDEBEC(Z,k);
DUDE = BER(X,Y);
DUDE_list = [BER(X,DUDEBEC(Z,i)) for i in range(10)]
DUDE_opt = min(DUDE_list)

