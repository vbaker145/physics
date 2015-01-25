#Math Phys HW5 problem 2
import numpy as np
import scipy as sp
from matplotlib import pyplot as pp

#Runge-Kutta 1 (Euler method) numerical integration
#    y - Function value at start (vector)
#    h - step size
#    dFunc - derivatives of y's (physics of the problem)
#         dFunc(x, y) - Returns derivatives of y's at x's 
def RK1(y, h, dFunc):
	yout = np.array(np.zeros(len(y)))
	dy = dFunc(y);
	yout = y + h*dy;
	return yout
		

#Runge-Kutta 2 (midpoint method) numerical integration
#    y - Function value at start (vector)
#    h - step size
#    dFunc - derivatives of y's (physics of the problem)
#         dFunc(x, y) - Returns derivatives of y's at x's 
def RK2(y, h, dFunc):
	yout = np.array(np.zeros(len(y)))
	dy = dFunc(y)
	ymid = y + (h/2.0)*dy
	dym = dFunc(ymid)
	yout = y + h*dym
	return yout

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

#Predictor-corrector numerical integration
#    y - Function value at start (vector)
#    h - step size
#    dFunc - derivatives of y's (physics of the problem)
#         dFunc(y) - Returns derivatives of y's at x's 
def PredCorr(y, h, dFunc):
	yout = np.array(np.zeros(len(y)))
	h2 = h*h
	dy = dFunc(y)
	#Predictor
	yout[0] = y[0] + dy[0]*h + 0.5*dy[2]*h2;
	yout[1] = y[1] + dy[1]*h + 0.5*dy[3]*h2;
	yout[2] = y[2] + dy[2]*h - 0.5*dy[2]*h;
	yout[3] = y[3] + dy[3]*h - 0.5*dy[3]*h;
	#Corrector
	dyend = dFunc(yout)
	yout[2] = yout[2] + 0.5*dyend[2]*h
	yout[3] = yout[3] + 0.5*dyend[3]*h
	return yout

#1/r2 central force law
#y - y values (vector of {x1,x2,v1,v2})
#
#Returns derivative values (dx/dt, dv/dt)
def r2law(y):
	dy = np.array(np.zeros(4))
	d = np.sqrt(y[0]*y[0]+y[1]*y[1])
	dy[0] = y[2]; #dx/dt is v
	dy[1] = y[3];
	dy[2] = -y[0]/pow(d,3.0) #dv/dt given by 1/r2 force
	dy[3] = -y[1]/pow(d,3.0)
	return dy

#1/r2 energy function
def energyR2(y):
	v2 = y[2]*y[2] + y[3]*y[3]
	r  = np.sqrt(y[0]*y[0] + y[1]*y[1])
	return (0.5*v2 - 1.0/r)

#Main function, numerical integration of two-body system
y = np.array(np.zeros(4))
y[0] = 1
y[1] = 0
y[2] = 0
y[3] = .75

#Fixed time step
dt = 0.00125
tmax = 100.0
nsteps = int(tmax/dt)
print 'nSteps '
print nsteps

#Positions
p = np.matrix([np.zeros(nsteps), np.zeros(nsteps)])

#Energy
e = np.array(np.zeros(nsteps))
e0   = energyR2(y)
print e0
eMax = 0.0

for jj in range(nsteps):
	p[0,jj] = y[0]
	p[1,jj] = y[1]
	e[jj]   = energyR2(y)
	y       = RK2(y,dt,r2law)
	e2      = energyR2(y) - e0
	if abs(e2) > eMax:
		eMax = abs(e2)

print 'eMax '
print eMax

pp.figure(1)
pp.subplot(2,1,1)
pp.plot(p[0,:], p[1,:],'-x')
pp.axis('equal')
tstr = 'Orbit plot, t=' + str(dt)
pp.title(tstr)
pp.subplot(2,1,2)
pp.plot(e[:], 'x')
round(eMax,4)
tstr = 'Max energy error ' + str(eMax)
pp.title(tstr)
pp.show()










