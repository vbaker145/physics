#Math Phys HW5 problem 1
import numpy as np
import scipy as sp
from matplotlib import pyplot as pp


#Fourier series expansion of sawtooth function
#(2/pi)*sum[ (1/n) sin pi*n*x ]
def fourierSaw(dx, nterms):
	nsteps = int(1./dx)
	xv = np.array(range(nsteps)) / float(nsteps)
	n = np.array(range(nterms)) + 1.
	ninv = 1. / n
	fx = np.array(np.zeros(nsteps))
	for jj in range(nsteps):
		fx[jj] = (2./np.pi)*np.sum(ninv*np.sin(np.pi*n*xv[jj])	)
	
	return fx

#Main function, numerical integration of two-body system
dx = 0.0005
nt = np.array([1,5,10,20,50,100,500])
nsteps = int(1./dx)
xv = np.array(range(nsteps)) / float(nsteps)

for nidx in range(len(nt)):
	f = fourierSaw(dx, nt[nidx])
	os = round(np.max(f),4)
	osidx = np.argmax(f)
	pp.figure(nidx+1)
	pp.plot(xv, f, '-x')
	tstr = str(nt[nidx]) + ' terms, max overshoot ' + str(os) + ' at ' + str(osidx)
	pp.title(tstr)

pp.show()










