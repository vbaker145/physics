#!/usr/bin/env python

# Solution to PHYS 502 Homework 1, Problems 1b,c.

from math import pi
import numpy as np
import matplotlib.pyplot as plt

def linplot(plot, l, nx, a):
    plot.set_xlim(0,nx)
    plot.set_title(l)
    plot.plot(a)

def logplot(plot, l, nx, a):
    plot.set_xscale('log')
    plot.set_yscale('log')
    plot.set_xlim(1,nx)
    plot.set_title(l)
    plot.plot(np.abs(a))

#----------------------------------------------------------------------
# Create the data.

N = 32*1024
rr = 2*np.random.random(N) - 1

# Specify the window function.

which = 0

window = np.ones(N)
if which == 1:
    window[:N/2] = np.arange(N/2)*2./N
    window[N/2:] = (N/2.-np.arange(N/2))*2./N
elif which == 2:
    window = 0.5*(1-np.cos(2*pi*np.arange(N)/(N-1.)))

wfac = N/np.sum(window**2)

# Compute the transform R and the power spectrum P.

r = rr*window
R = np.fft.fft(r)

P = np.zeros(N/2+1)
P[1:N/2] = np.abs(R[1:N/2])**2 + np.abs(R[N-1:N/2:-1])**2
P[N/2] = abs(R[N/2])**2
P *= wfac

# Average the power spectrum: A.

nav = 64
nav2 = nav/2
A = np.zeros(N/2+1)
for k in range(nav2, N/2-nav2):
    A[k] = np.sum(P[k-nav2:k+nav2])/nav

# Compute the random walk w.

w = np.zeros(N)
for j in range(1,N):
    w[j] = w[j-1] + rr[j]

# Compute the transform W and the power spectrum Q.

w *= window
W = np.fft.fft(w)

Q = np.zeros(N/2+1)
Q[1:N/2] = np.abs(W[1:N/2])**2 + np.abs(W[N-1:N/2:-1])**2
Q[N/2] = abs(W[N/2])**2
Q *= wfac

# Average the power spectrum: B.

nav = 64
nav2 = nav/2
B = np.zeros(N/2+1)
for k in range(nav2, N/2-nav2):
    B[k] = np.sum(Q[k-nav2:k+nav2])/nav

#----------------------------------------------------------------------
# Plot the results.

fig = plt.figure(figsize=(10,5))
fig.canvas.set_window_title('PHYS 502, Homework 1.2')
plt.subplots_adjust(hspace=0.5, wspace=0.25, \
                    left=0.06, right=0.975, bottom=0.1, top=0.925)

plot1 = fig.add_subplot(2,4,1)
linplot(plot1, 'power P(k)', N/2, P)

plot2 = fig.add_subplot(2,4,5)
logplot(plot2, 'power P(k)',  N/2, P)
plot2.plot([1.e0,2.e4], [2.*N/3, 2.*N/3], 'r-')

plot3 = fig.add_subplot(2,4,2)
logplot(plot3, 'averaged power A(k)', N/2, A)
plot3.plot([1.e0,2.e4], [2.*N/3, 2.*N/3], 'r-')

plot4 = fig.add_subplot(2,4,6)
plot4.set_ylim(0.01,1.e6)
logplot(plot4, 'averaged power A(k)', N/2, A)
plot4.plot([1.e0,2.e4], [2.*N/3, 2.*N/3], 'r-')

plot5 = fig.add_subplot(2,4,3)
linplot(plot5, 'random walk w(j)', N, w)

plot6 = fig.add_subplot(2,4,4)
logplot(plot6, 'power Q(k)', N/2, Q)

plot7 = fig.add_subplot(2,4,8)
plot7.set_ylim(0.1,1.e12)
logplot(plot7, 'averaged power B(k)', N/2, B)
plot7.plot([1.e2,1.e4], [1.e7, 1.e3], 'r-')
plot7.text(4.e1, 3.e3, 'slope = -2', color='red')

plt.show()
