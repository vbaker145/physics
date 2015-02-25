#!/usr/bin/env python

# Solution to PHYS 502 Homework 1, Problem 2.

import sys
import numpy as np
import matplotlib.pyplot as plt

def read_data(file):
    c = []
    f = open(file)
    for l in f.readlines():
        c.append(float(l.split()[1]))
    f.close()
    return np.array(c)

def linplot(plot, l, nx, a):
    plot.set_xlim(0,nx)
    plot.set_title(l)
    plot.plot(np.abs(a))

def logplot(plot, l, nx, a):
    plot.set_xscale('log')
    plot.set_yscale('log')
    plot.set_xlim(1,nx)
    plot.set_title(l)
    plot.plot(np.abs(a))

def normalize(a):
    sum = np.sum(np.abs(a))
    return a/sum

#----------------------------------------------------------------------
# Read the data.

c = read_data('../corrupt.dat')
N = len(c)

# Compute the transform.

C = np.fft.fft(c)

# Filter the data.

ncut = 16	# last retained point
if len(sys.argv) > 1: ncut = int(sys.argv[1])
print 'cut at n =', ncut

F = np.array(C)
F[ncut+1:N-ncut] = 0

# Compute the Gaussian kernel, properly normalized and reordered...

width = 2048.
g = np.zeros(N)
g[:N/2+1] = np.exp(-(np.arange(N/2+1)/width)**2)
g[N/2+1:] = g[N/2-1:0:-1]	# careful!
g = normalize(g)

# ...and its transform.

G = np.fft.fft(g)
G = np.where(abs(G) < 1.e-20, 1.e-20, G)

# Deconvolve...

S = np.zeros(N, dtype=complex)
S = np.where(abs(F) > 0, F/G, S)

# ...and transform back.

s = np.fft.ifft(S)

#----------------------------------------------------------------------
# Plot the results.

fig = plt.figure(figsize=(12,5))
fig.canvas.set_window_title('PHYS 502, Homework 1.3')
plt.subplots_adjust(hspace=0.5, wspace=0.25, \
                    left=0.06, right=0.975, bottom=0.1, top=0.925)

plotc = fig.add_subplot(2,4,1)
linplot(plotc, 'corrupted c(k)', N, c)

plotC = fig.add_subplot(2,4,5)
logplot(plotC, 'transformed C(n)',  N/2, C)

plotg = fig.add_subplot(2,4,2)
linplot(plotg, 'response g(k)', N, g)

plotF = fig.add_subplot(2,4,6)
logplot(plotF, 'filtered F(n)', N/2, F)

plotG = fig.add_subplot(2,4,3)
logplot(plotG, 'response G(n)', N/2, G)

plotS = fig.add_subplot(2,4,7)
logplot(plotS, 'recovered S(n)', N/2, S)

plots = fig.add_subplot(2,4,4)
linplot(plots, 'recovered s(k)', N, s)
try:
    o = read_data('../original_uncorrupted.dat')
    linplot(plots, 'recovered s(k), original o(k)', N, o)
    plotu = fig.add_subplot(2,4,8)
    linplot(plotu, 'difference s(k)-o(k)', N, s-o)
except:
    pass

plt.show()
