import numpy as np
import scipy as sp

#Calculate the conductance matrix 
def ConductanceMatrix(N):
	cond = np.matrix([np.zeros(N-1),]*(N-1))
	inds = np.array(range(N-1))+2
	#First calculate element-by-element conductance from min(i,j)+2*max(i,j)
	#Elements are negative because we'll subtract Vj/Rij from Vi/Rij
	for i in range(N-1):
		for j in range(N-1):
			iv = inds[i]
			jv = inds[j]
			cond[i,j] = 1.0/(min(iv,jv) + 2*max(iv,jv))
	#Make diagonal element resistance sum of all resistance
	#Turns bunch of differences in sumR*Vi - sum(Vj/Rij)
	for i in range(N-1):
		cond[i,i] = 0.0
		#Want to use Python list.remove() to create sequence
		#Forces stupid indexing here
		l = range(N)
		l.remove(i+1)
		l = np.array(l)+1
		for j in l:
			iv = i+2.0
			jv = j 
			cond[i,i] -= 1.0/(min(iv,jv) + 2*max(iv,jv))
	return cond

N = 2  	#Number of nodes
cmat = ConductanceMatrix(N)
currents = np.zeros(N-1) 
currents[0] = 1 #Net current is zero at all nodes except 2 and 1
currents = np.array(currents)

#Now we have conductances (A) and currents (b), solve Ax = b for voltages (x)
v = np.linalg.solve(cmat, currents)

#R = V/I, I = -1 between node 2 and 1, so resistance is same as -voltage at node 2
print v[0]



		

