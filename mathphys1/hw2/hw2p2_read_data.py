import numpy as np
import scipy as sp

#Open file and read in x, y, sigma
dat = open('hw2.2.dat', 'r')
l = [line for line in dat.readlines()]
nrows = len(l)
idx = 0
x = np.array(np.zeros(nrows))
y = np.array(np.zeros(nrows))
sig = np.array(np.zeros(nrows)) 
for line in l:
	j = line.split()
	x[idx] = j[0]
	y[idx] = j[1]
	sig[idx] = j[2]
	idx += 1

#data is now stored in x, y and sig

