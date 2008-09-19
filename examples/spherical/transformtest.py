#improt system modules
import sys
import os

try:
	import pygsl.sf as sf
except:
	print "Warning: Could not imort pygsl.sf. pygsl is required for certain functions"

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("../.."))

#Load and reload pyprop in order to get recent changes
import pyprop
pyprop = reload(pyprop)

#numpy an pylab for good measure
from numpy import *
from pylab import *


def MapLmIndex(l, m):
	assert(abs(m) <= l)
	return (l + 1) * l + m

def MapIndexLm(index):
	l = int(sqrt(index))
	m = index - (l + 1) * l
	return l,m

class Transform:
	def __init__(self, lmax):
		self.lmax = lmax
		self.trans = self.GetTransform(lmax)
		self.lmsize = (1, self.trans.GetAssociatedLegendrePolynomial().shape[1])
		self.omegasize = (1, self.trans.GetOmegaGrid().shape[0])
		#do one transform to initialize fourier transform
		lm = self.CreateLMArray()
		omega = self.Sph2Grid(lm)
		lm2 = self.Grid2Sph(omega)
		

	def CreateLMArray(self):
		return zeros(self.lmsize[1], dtype=complex)

	def CreateOmegaArray(self):
		return zeros(self.omegasize[1], dtype=complex)

	def GetTransform(self, lmax):
		trans = pyprop.core.SphericalTransformTensorGrid()
		trans.Initialize(lmax)
		return trans

	def Sph2Grid(self,  lmdata):
		indata = zeros(self.lmsize, dtype=complex)
		outdata = zeros(self.omegasize, dtype=complex)
		indata[0,:] = lmdata
		self.trans.InverseTransform(indata, outdata)
		return outdata[0,:]		

	def Grid2Sph(self,  griddata):
		indata = zeros(self.omegasize, dtype=complex)
		outdata = zeros(self.lmsize, dtype=complex)
		indata[0,:] = griddata
		self.trans.ForwardTransform(indata, outdata)
		return outdata[0,:]

	def GetEigenfunction(self, l, m):
		lmdata = self.CreateLMArray()
		lmdata[:] = 0 
		lmdata[MapLmIndex(l,m)] = 1;
		omegaData = self.Sph2Grid(lmdata)

		phiCount = 2 * self.lmax + 1
		ret = fromfunction(lambda theta, phi: real(omegaData[theta*phiCount + phi]), (phiCount, phiCount))
		return ret

	def GetGrid(self):
		return self.trans.GetThetaGrid(), self.trans.GetPhiGrid()

def Ylm(l, m, theta):
	ylm = zeros(len(theta), dtype=complex)
	for i in range(len(theta)):
		curx =  cos(theta[i])
		val  = sf.legendre_sphPlm(l,abs(m),curx)
		ylm[i] = val[0] + val[1] * 1.0j
	return ylm

def TestTransform(l,m, trans=None, lmax=8):
	if trans == None:
		trans = Transform(lmax)

	a = trans.CreateLMArray()
	a[MapLmIndex(l,m)] = 1.0
	b = trans.Sph2Grid(a)
	c = trans.Grid2Sph(b)
	val = c[MapLmIndex(l,m)].real
	c[MapLmIndex(l,m)] = 0
	err = sqrt(sum(abs(c)**2))
	return val, err

def TestAllLm(lmax):
	trans = Transform(lmax)

	maxErr = -1.
	maxErrLm = -1, 0

	normErr = -1.
	normErrLm = -1, 0
	

	for l in r_[0:lmax+1]:
		for m in r_[-l:l+1]:
			norm, err = TestTransform(l, m, trans)
			if err > maxErr:
				maxErr = err
				maxErrLm = l, m

			if abs(norm - 1.0) > normErr:
				normErr = abs(norm - 1.0)
				normErrLm = l, m
			
			print "l, m, norm, err", l, m, norm, err 

	print "Max Error at ", maxErrLm, ". Error: ", maxErr
	print "Max Norm Error at ", normErrLm, ". Error: ", normErr

def TestAllLm2(lmax):
	trans = Transform(lmax)
	theta, phi = trans.GetGrid()
	
	maxErr = -1.
	maxErrLm = -1, 0

	for l in r_[0:lmax+1]:
		for m in r_[-l:l+1]:
			ylm1 = trans.GetEigenfunction(l,m)[:,0]
			ylm2 = real(Ylm(l,m,theta))

			err = sum(abs(ylm1 - ylm2))
			print "l, m, err", l, m, err
			if err > maxErr:
				maxErr = err
				maxErrLm = l,m
				maxYlm = ylm1, ylm2

	print "Max Error at ", maxErrLm, ". Error: ", maxErr
	return maxYlm


def TestOrthogonalityAssocLegendre(lmax=4):
	trans = pyprop.core.SphericalTransformTensorGrid()
	trans.Initialize(lmax)

	legendre = trans.GetAssociatedLegendrePolynomial()
	weights = trans.GetWeights()

	thetaCount = legendre.shape[0]
	lmCount = legendre.shape[1]

	assert(thetaCount == len(weights))
		
	for i in range(lmCount):
		for j in range(lmCount):
			c = dot( conj(legendre[:,i]), weights * legendre[:,j] )
			print "%i, %i, %f" % (i,j,c)


def TestOrthogonalitySphericalHarmonic(lmax=4, tolerance=1e13):
	trans = pyprop.core.SphericalTransformTensorGrid()
	trans.Initialize(lmax)

	theta = trans.GetThetaGrid()
	phi = trans.GetPhiGrid()

	sph = trans.GetSphericalHarmonic()
	w = trans.GetWeights()
	weights = array([ [w[i]]*len(phi) for i in range(len(theta)) ], dtype=double).flatten()

	omegaCount = sph.shape[0]
	lmCount = sph.shape[1]

	assert(omegaCount == len(weights))
	
	maxError = 0

	for i in range(lmCount):
		for j in range(lmCount):
			c = dot( conj(sph[:,i]), weights * sph[:,j] )
			if (i == j):
				error =  abs(c - 1) 
			else:
				error =  abs(c) 
			
			if error > maxError:
				maxError = error
			if error > tolerance:
				print "%i, %i, %g" % (i,j,error)

	print "Max Error = %g" % (maxError)
	return maxError < tolerance

def TestDifferentiationSphericalHarmonic(lmax=4, tolerance=1e-13):
	trans = pyprop.core.SphericalTransformTensorGrid()
	trans.Initialize(lmax)

	theta = trans.GetThetaGrid()
	phi = trans.GetPhiGrid()

	sph = trans.GetSphericalHarmonic()
	sphDiff = trans.GetSphericalHarmonicDerivativeTheta()
	w = trans.GetWeights()
	weights = array([ [w[i]]*len(phi) for i in range(len(theta)) ], dtype=double).flatten()
	th = array([ [theta[i]]*len(phi) for i in range(len(theta)) ], dtype=double).flatten()

	omegaCount = sph.shape[0]
	lmCount = sph.shape[1]

	assert(omegaCount == len(weights))

	#return weights, sph, sphDiff
	
	maxError = 0

	for i in range(lmCount):
		l,m = MapIndexLm(i)
		for j in range(lmCount):
			l2,m2  = MapIndexLm(j)
			c = dot( conj(sph[:,i]), weights * sphDiff[:,j] )
			if abs(c) < tolerance:
				c = 0
			if c != 0:
				print "(%i,%i), (%i,%i), %g" % (l,m,l2,m2,c)
			"""
			if (i == j):
				error =  abs(c - 1) 
			else:
				error =  abs(c) 
			
			if error > maxError:
				maxError = error
			if error > tolerance:
				print "%i, %i, %g" % (i,j,error)
			"""


def TestDifferentiationSphericalHarmonicPhi(lmax=4, tolerance=1e-13):
	trans = pyprop.core.SphericalTransformTensorGrid()
	trans.Initialize(lmax)

	theta = trans.GetThetaGrid()
	phi = trans.GetPhiGrid()

	sph = trans.GetSphericalHarmonic()
	sphDiff = trans.GetSphericalHarmonicDerivativePhi()
	w = trans.GetWeights()
	weights = array([ [w[i]]*len(phi) for i in range(len(theta)) ], dtype=double).flatten()
	th = array([ [theta[i]]*len(phi) for i in range(len(theta)) ], dtype=double).flatten()

	omegaCount = sph.shape[0]
	lmCount = sph.shape[1]

	assert(omegaCount == len(weights))

	#return weights, sph, sphDiff
	
	maxError = 0

	for i in range(lmCount):
		l,m = MapIndexLm(i)
		for j in range(lmCount):
			l2,m2  = MapIndexLm(j)
			c = dot( conj(sph[:,i]), weights * sphDiff[:,j] )
			print "(%i,%i), (%i,%i), %s" % (l,m,l2,m2,c)
			if abs(c) < tolerance:
				c = 0
			if abs(c) > tolerance:
				pass #print "(%i,%i), (%i,%i), %g" % (l,m,l2,m2,c)
			"""
			if (i == j):
				error =  abs(c - 1) 
			else:
				error =  abs(c) 
			
			if error > maxError:
				maxError = error
			if error > tolerance:
				print "%i, %i, %g" % (i,j,error)
			"""




def TestDifferentiationReducedSphericalHarmonic(lmax=4, tolerance=1e-13):
	trans = pyprop.core.ReducedSphericalTools()
	trans.Initialize(lmax)

	theta = trans.GetThetaGrid()

	sph = trans.GetAssociatedLegendrePolynomial()
	sphDiff = trans.GetAssociatedLegendrePolynomialDerivative()
	weights = trans.GetWeights()

	omegaCount = sph.shape[0]
	lmCount = sph.shape[1]

	assert(omegaCount == len(weights))

	#return weights, sph, sphDiff
	
	maxError = 0

	for i in range(lmCount):
		for j in range(lmCount):
			c = dot( conj(sph[:,i]), weights * sphDiff[:,j] )
			if abs(c) < tolerance:
				c = 0
			if c != 0:
				print "(%i), (%i), %g" % (i,j,c)

