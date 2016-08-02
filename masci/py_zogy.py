############################################################
# Python implementation of ZOGY image subtraction algorithm
# See Zackay, Ofek, and Gal-Yam 2016 for details
# http://arxiv.org/abs/1601.02655
# SBC - 6 July 2016
############################################################

import sys
import numpy as np
import numpy.fft as fft
import astropy.io.fits as fits
from scipy.fftpack import fft2, ifft2

def py_zogy(newimage, refimage, newpsfimage, refpsfimage, newsigimage,
            refsigimage, subimage, subpsfimage, subcorrimage, 
            newgain=1.539, newbkg=303.922, refgain=1.539, 
            refbkg=309.590, newsigma=52.764, refsigma=2.259, Eps=1.0e-8):

	'''Python implementation of ZOGY image subtraction algorithm.
	As per Frank's instructions, will assume images have been aligned,
	background subtracted, and gain-matched.'''
	
	# Load the new and ref images into memory
	newim = fits.open(newimage)[0].data
	refim = fits.open(refimage)[0].data
	
	# Load the variance images into memory
	newsig = fits.open(newsigimage)[0].data
	refsig = fits.open(refsigimage)[0].data
	
	# Pad the images (if necessary)
	[ny, nx] = refim.shape
	if ny % 2 == 0: 
		PadY = 1
	else: 
		PadY = 0
	if nx % 2 == 0:
		PadX = 1
	else:
		PadX = 0
	newim = np.pad(newim, [(0, PadY), (0, PadX)], 'constant', 
	               constant_values=0.0)
	newsig = np.pad(newsig, [(0, PadY), (0, PadX)], 'constant', 
	               constant_values=0.0)
	refim = np.pad(refim, [(0, PadY), (0, PadX)], 'constant', 
	               constant_values=0.0)
	refsig = np.pad(refsig, [(0, PadY), (0, PadX)], 'constant', 
	               constant_values=0.0)
	
	# Load the PSFs into memory
	newpsf = fits.open(newpsfimage)[0].data
	refpsf = fits.open(refpsfimage)[0].data
	
	# Pad the PSFs (if necessary)
	PadImY = int(0.5 * (newim.shape[0] - (refpsf.shape[0]-1) - 1))
	PadImX = int(0.5 * (newim.shape[1] - (refpsf.shape[1]-1) - 1))
	newpsf = np.pad(newpsf, [(PadImY, PadImY), (PadImX, PadImX)],
	                'constant', constant_values=0.0)
	refpsf = np.pad(refpsf, [(PadImY, PadImY), (PadImX, PadImX)],
	                'constant', constant_values=0.0)
	                
	# Shift the PSF to the origin so it will not introduce a shift
	newpsf = fft.fftshift(newpsf)
	refpsf = fft.fftshift(refpsf)
	
	# Take Fourier Transform of all arrays
	Fnewim = fft2(newim.astype(float))
	Frefim = fft2(refim.astype(float))
	Fnewpsf = fft2(newpsf.astype(float))
	Frefpsf = fft2(refpsf.astype(float))
	
	# Numerator of FT of difference image (Equation 13)
	# Need element-by-element multiplication here
	FPrFN = Frefpsf * Fnewim
	FPnFR = Fnewpsf * Frefim
	
	# Denominator of FT of difference image (Equation 13)
	# Need element-by-element multiplication here
	FPn2 = Fnewpsf * np.conjugate(Fnewpsf)
	FPr2 = Frefpsf * np.conjugate(Frefpsf)
	Den = refsigma**2 * FPn2 + newsigma**2 * FPr2 + Eps
	SqrtDen = np.sqrt(Den)
	
	# Subtraction image (and its FT)
	FD = (FPrFN - FPnFR) / SqrtDen
	D = ifft2(FD) 
	
	# Remove padding
	D = D[:-PadY,:-PadX]
	
	# Write out resulting subtraction image
	outim = fits.open(newimage)
	outim[0].data = D.real
	outim.writeto(subimage, output_verify="warn")
	
	# FT of PSF of subtraction image
	F_D = 1.0 / np.sqrt(newsigma**2 + refsigma**2)
	Fsubpsf = Fnewpsf * Frefpsf / SqrtDen / F_D
	
	# PSF of subtraction image
	subpsf = ifft2(Fsubpsf)
	subpsf = fft.ifftshift(subpsf)
	subpsf = subpsf[PadImY:-PadImY,PadImX:-PadImX]
	outpsf = fits.open(newpsfimage)
	outpsf[0].data = subpsf.real
	outpsf.writeto(subpsfimage, output_verify="warn")
	
	# Calculate Score image
	# ***Note this is formally off by a factor of F_r***
	FS = (np.conjugate(Fnewpsf) * FPr2 * Fnewim - 
	      np.conjugate(Frefpsf) * FPn2 * Frefim) / Den
	S = ifft2(FS)
	
	# Alternate method to get score image
	# Use this for now - seem identical (to numerical precision)
	FS2 = F_D * FD * np.conjugate(Fsubpsf)
	S2 = ifft2(FS2)
	
	# Now for the noise image (Snoise)
	# First term is the variance in the new image (in e-)
	Vnew = np.power(newsig,2)
	FVnew = fft2(Vnew)
	Fkn = np.conjugate(Fnewpsf) * FPr2 / Den # Off by a factor of Fn
	FVSn = FVnew * fft2(ifft2(Fkn)**2)
	
	# Second is the variance in the ref image
	Vref = np.power(refsig,2)
	FVref = fft2(Vref)
	Fkr = np.conjugate(Frefpsf) * FPn2 / Den # Off by a factor of Fr
	FVSr = FVref * fft2(ifft2(Fkr)**2)
	
	# Calculate Snoise (leaving out astrometric noise, gradient, ...)
	Snoise = np.sqrt( ifft2(FVSn) + ifft2(FVSr) )
	
	# Scorr: calculate and write out
	Scorr = S2 / Snoise
	#Scorr = fft.ifftshift(Scorr)
	Scorr = Scorr[:-PadY,:-PadX]
	scorrim = fits.open(newimage)
	scorrim[0].data = Scorr.real
	scorrim.writeto(subcorrimage, output_verify="warn")
	
	return
	
if __name__ == "__main__":

	# Usage
	if len(sys.argv) != 10:
		print "Usage: python py_zogy.py <newimage> <refimage> <newpsfimage> <refpsfimage> <newsigmaimage> <refsigmaimage> <subimage> <subpsfimage> <subcorrimage>"
	else:
		py_zogy(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
		        sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9])

