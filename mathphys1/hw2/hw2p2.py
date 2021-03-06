import numpy as np
import scipy as sp
import scipy.optimize as optimization
from matplotlib import pyplot as pp


def chi2(_N, _M, _alpha, _x, _y, _sig):
	c2 = 0.0
	for i in range(_N):
		msum = _y[i]
		for j in range(_M):
			msum -= _alpha[j]*_x[i]**j
		msum /= _sig[i]
		c2 += msum**2
	return c2
	
def approx(_N, _M, _x, _y, _sig):
	#Form x^0 column
	x0 = np.array(np.ones(_N)) / _sig

	#Create coefficient matrix
	A = np.matrix([x0]*_M)
	A = A.transpose()

	#Fill in coefficient matrix
	for j in range(_M-1):
		tx = (_x**(j+1))/_sig
		A[:,j+1] = np.matrix(tx).transpose()

	#SVD on A (note that "V" is already transposed
	U, S, V = np.linalg.svd(A, False, True)
	print "Smallest W:"
	print min(S)

	#Check SVD decomp by reconstructing A 
	A2 = U*np.diag(S)*V
	SI = np.diag(1/S) #W "inverse"
	AI = V.transpose()*SI*U.transpose()

	ysig = np.matrix(_y)/_sig;
	alpha = AI*ysig.transpose()
	return alpha


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

am = 2

alpha = approx(nrows, am, x, y, sig)
print 'Solution: '
print alpha

print 'Chi2: '
print chi2(nrows, am, alpha, x, y, sig)

pp.figure(10)
pp.plot(x, y, 'x')
xv = np.linspace(-1,1,1000)
fv = 0
for i in range(am):
	fv += float(alpha[i])*xv**i
pp.plot(xv,fv)
pp.show()

