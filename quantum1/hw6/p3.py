#Quantum 1 HW6 problem 3
import numpy as np
import scipy as sp
from matplotlib import pyplot as pp

#Runge-Kutta 4 numerical integration
#    y - Function value at start (vector)
#    h - step size
#    dFunc - derivatives of y's (physics of the problem)
#         dFunc(x, y) - Returns derivatives of y's at x's 
def RK4(y, t, h, dFunc):
	yout = np.array(np.zeros(len(y)))
	h2 = h/2.0
	h6 = h/6.0
	#First step
	dy = dFunc(y,t)
	yt = y + h2*dy 
	#Second step
	dyt = dFunc(yt,t)
	yt  = y + h2*dyt
	#Third step
	dym = dFunc(yt,t)
 	yt = y + h*dym
	dym = dym + dyt
	#Fourth step
	dyt = dFunc(yt,t)
	yout = y + h6*(dy+dyt+2.0*dym)
	return yout

#Rabi problem
#y - y values (vector of {a,b,c,d})
#
#Returns derivative values (da,db,dc,dd)
def rabi(y,t):
	eps   = 1
	gamma = 1
	omega = 1.5
	dy = np.array(np.zeros(len(y)))
	dy[0] = eps*y[1]+gamma*np.cos(omega*t)*y[3]
	dy[1] = -eps*y[0]-gamma*np.cos(omega*t)*y[2]
	dy[2] = gamma*np.cos(omega*t)*y[1]-eps*y[3]
	dy[3] = -gamma*np.cos(omega*t)*y[0]+eps*y[2]
	return dy

#Main function, numerical integration of Lorenz system
y = np.array(np.zeros(4))
y[0] = 0
y[1] = 0
y[2] = 1
y[3] = 0
dt = 0.001
tmax = 10.0
nsteps = int(tmax/dt)
t = np.linspace(0,tmax, nsteps)

#Values
p = np.matrix([np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps)])
a = np.array(np.zeros(nsteps))
b = np.array(np.zeros(nsteps))
c = np.array(np.zeros(nsteps))
d = np.array(np.zeros(nsteps))


for jj in range(nsteps):
	a[jj] = y[0]
	b[jj] = y[1]
	c[jj] = y[2]
	d[jj] = y[3]
	y     = RK4(y,t[jj],dt,rabi)

pp.figure(1)
pp.hold(True)
pp.plot(t, a, 'b', label='a')
pp.plot(t, b, 'r', label='b')
pp.plot(t, c, 'g', label='c')
pp.plot(t, d, 'y', label='d')
pp.xlabel('Time (sec)')
pp.ylabel('Value')
pp.title('ABCD evolution')

pp.figure(2, figsize=(6,4))
p1=np.array(np.zeros(nsteps))
p2=np.array(np.zeros(nsteps))
for jj in range(nsteps):
	p1[jj] = a[jj]*a[jj]+b[jj]*b[jj]
	p2[jj] = c[jj]*c[jj]+d[jj]*d[jj]
pp.plot(t, p1, 'b', label='P10')
pp.plot(t, p2, 'r', label='P01')
pp.xlabel('time')
pp.ylabel('probability')
pp.title('Probability')
pp.legend()

pp.show()










