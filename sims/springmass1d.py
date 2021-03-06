#Math Phys HW5 problem 3
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

#Lorenz problem
#y - y values (vector of {x,y,z})
#
#Returns derivative values (dx, dy, dz)
def lorenzLaw(y):
	sigma = 10
	b     = 8./3.
	r     = 30
	dy = np.array(np.zeros(len(y)))
	dy[0] = -sigma*y[0] + sigma*y[1]
	dy[1] = r*y[0] - y[1] - y[0]*y[2]
	dy[2] = -b*y[2] + y[0]*y[1]
	return dy


#Main function, numerical integration of Lorenz system
y = np.array(np.zeros(3))
y[0] = 5
y[1] = 6
y[2] = 10
dt = 0.01
tmax = 100.0
nsteps = int(tmax/dt)
t = np.linspace(0,tmax, nsteps)

#Positions
p = np.matrix([np.zeros(nsteps), np.zeros(nsteps), np.zeros(nsteps)])
x1 = np.array(np.zeros(nsteps))
x2 = np.array(np.zeros(nsteps))
x3 = np.array(np.zeros(nsteps))


for jj in range(nsteps):
	x1[jj] = y[0]
	x2[jj] = y[1]
	x3[jj] = y[2]
	y       = RK4(y,dt,lorenzLaw)

pp.figure(1)
pp.hold(True)
pp.plot(t, x1, 'b', label='x')
pp.plot(t, x2, 'r', label='y')
pp.plot(t, x3, 'g', label='z')
pp.xlabel('Time (sec)')
pp.ylabel('Value')
pp.title('XYZ evolution for Lorenz system')
pp.legend()
pp.figure(2)
pp.plot(x1, x3)
pp.xlabel('X')
pp.ylabel('Z')
pp.title('X-Z cut of Lorenz attractor')
pp.show()










