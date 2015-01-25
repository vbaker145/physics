#Math Phys HW5 problem 4
import numpy as np
import scipy as sp
from matplotlib import pyplot as pp

#Runge-Kutta 4 numerical integration
#    y - Function value at start (vector)
#    h - step size
#    dFunc - derivatives of y's (physics of the problem)
#         dFunc(x, y) - Returns derivatives of y's at x's 
def RK4(y, h, dFunc):
	yout = np.array(np.zeros(len(y)))
	h2 = h/2.0
	h6 = h/6.0
	#First step
	dy = dFunc(y)
	yt = y + h2*dy 
	#Second step
	dyt = dFunc(yt)
	yt  = y + h2*dyt
	#Third step
	dym = dFunc(yt)
 	yt = y + h*dym
	dym = dym + dyt
	#Fourth step
	dyt = dFunc(yt)
	yout = y + h6*(dy+dyt+2.0*dym)
	return yout

#P4 system
#y - y values
#
#Returns derivative values
def p4Law(y):
	dy = np.array(np.zeros(len(y)))
	dy[0] = y[1]
	dy[1] = -y[1]-50*y[0]
	return dy

#Integrates the P4 system with a given dy(0)
# plot - True/False, plot y(x) from 0-1
def integrateP4(dy, plot):
	y = np.array(np.zeros(2))
	y[0] = 2.0 #Fixed initial condition
	y[1] = dy

	dx = 0.01
	xmax = 1.0
	nsteps = int(xmax/dx)
	x = np.linspace(0,xmax, nsteps)

	if plot:
		yv = np.zeros(nsteps+1)
		yv[0] = 2.0

	for jj in range(nsteps):
		y       = RK4(y,dx,p4Law)
		if plot:
			yv[jj+1] = y[0]

	if plot:
		lstr = str(dy)
		xv = np.array(range(nsteps+1)) / float(nsteps+1)
		pp.plot(xv, yv, '-', label=lstr)

	return y[0] + 2

#Main function
#Finds dy(0) using bisection
dy1 = -1000.0
dy2 = 1000.0
em = 1.0

e1 = integrateP4(dy1, False)
e2 = integrateP4(dy2, False)
plotCnt = 5

pp.figure(1)
pp.hold(True)

while abs(em) > 0.0001:
	dym = (dy2+dy1)/2.0
	if plotCnt == 5:
		em = integrateP4(dym, True)
		plotCnt = 0
	else:
		em = integrateP4(dym, False)
		plotCnt = plotCnt + 1
	print 'dy1 dym dy2'
	print [dy1, dym, dy2]
	print 'E1 EM E2'
	print [e1, em, e2]
	if em*e1 < 0:
		#Root in first interval	
		dy2 = dym	
		e2  = em
	else:
		#Root in second interval
		dy1 = dym
		e1  = em

#Plot final solution
em = integrateP4(dym, True)
pp.legend()
pp.xlabel('X')
pp.ylabel('Y')
pp.show()










