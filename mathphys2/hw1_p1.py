#Math Phys 2 HW1 problem 1
import numpy as np
import scipy as sp
import random as rnd
from matplotlib import pyplot as pp

def magplot(w, psd, plottype):
	pp.figure(figsize=(10,12))
	pp.subplot(2,1,1)
	pp.ylabel('power (dB)')
	pp.xlabel('frequency (dB)')
	pp.subplot(2,1,2)
	pp.ylabel('average power (dB)')
	pp.xlabel('frequency (dB)')

	if 'log' in plottype:
		pp.subplot(2,1,1)
		pp.plot(np.log10(w), np.log10(psd))
		pp.title('PSD')

	if 'lin' in plottype:
		pp.subplot(2,1,1)
		pp.plot(np.log10(w), psd)
		pp.ylabel('power (linear)')
		pp.subplot(2,1,2)
		pp.plot(np.log10(w), np.log10(psd))
		pp.ylabel('power (dB)')

	if 'avg' in plottype:
		winSize  = 64
		winSize2 = winSize/2
		dLen = psd.size
		psdtmp = np.zeros(dLen+winSize)
		psdtmp[0:winSize2] = psd[0:winSize2]
		psdtmp[winSize2:winSize2+dLen] = psd
		psdtmp[-winSize2:] = psd[-winSize2:]
		psdavg = np.zeros(dLen)
		for i in range(dLen):
			psdavg[i] = np.sum(psdtmp[i+winSize2:i+winSize2+winSize])/dLen
		pp.subplot(2,1,2)
		pp.plot((np.log10(w)), np.log10(psdavg))
		
		pp.title('PSD (averaged '+str(winSize)+')')

def testFFT(npts, dtype):
	data = np.zeros(npts)
	k = np.array(range(npts/2)) + 1
	if 'uniform' in dtype:
		for idx in range(npts):
			data[idx] = rnd.random()*2 - 1.0
	if 'walk' in dtype:
		dsum = 0.0
		for idx in range(npts):
			data[idx] = dsum
			dsum = dsum + rnd.random()*2 - 1.0	
		pp.figure();
		pp.plot(data)
	win = np.hamming(npts)
	data = data * win
	dfft = np.fft.fft(data)
	dfft_c = np.conj(dfft)
	dMag = np.zeros(npts/2)
	for idx in k:
		dMag[idx-1] = abs(dfft[idx]*dfft_c[idx] + dfft[npts-idx]*dfft_c[npts-idx])
 
	print np.mean(dMag)
	print np.mean(dMag)/npts
	print np.var(dMag)
	print np.var(dMag)/(npts*npts)
	magplot(k, dMag, 'log avg')

#Main fxn
#testFFT(512, 'uniform')
#testFFT(2048, 'uniform')
testFFT(32768, 'walk')
pp.show()

