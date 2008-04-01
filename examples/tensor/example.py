import sys

sys.path.append("./pyprop")
import pyprop

from numpy import array
from numpy import complex
from numpy import zeros

execfile("TensorGenerator.py")

#------------------------------------------------------------------------------------
#                       Test functions
#------------------------------------------------------------------------------------


from libpotential import *

def test():
	conf = pyprop.Load("test.ini")
	gen = TensorPotentialGenerator(config=conf)
	A = gen.GeneratePotential(conf, conf.TestPotential)
	B = gen.GeneratePotential(conf, conf.KineticEnergy)
	Overlap = gen.GeneratePotential(conf, conf.OverlapMatrix)
	
	#A = MultiplyInverseOverlapMatrix(A, Overlap, 0)
	A = MultiplyInverseOverlapMatrix(A, Overlap, 1)
	#B = MultiplyInverseOverlapMatrix(B, Overlap, 0)
	B = MultiplyInverseOverlapMatrix(B, Overlap, 1)

	return A, B, Overlap

def test2():
	conf = pyprop.Load("config.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	prop.AdvanceStep()
	return prop

import numpy.linalg

def MultiplyInverseOverlapMatrix(A, OverlapMatrix, bsplineRank):
	rank = len(A.shape)

	if rank == 1:
		S = OverlapMatrix
		s = UnpackMatrix(S)
		sInv = numpy.linalg.inv(S)
		a = UnpackMatrix
		b = dot(sInv, a)
		B = b.reshape((A.shape[0]))

	if rank == 2:
		B = A.copy()
		if bsplineRank == 1:
			for otherIndex in range(A.shape[0]):
				Aslice = A[otherIndex, :]
				a = UnpackMatrix(Aslice)
				S = OverlapMatrix[otherIndex, :]
				if norm(a) > 10e-5:
					s = UnpackMatrix(S)
					try:
						sInv = numpy.linalg.inv(s)
					except:
						print otherIndex
						return B
					
					b = dot(sInv, a)
					B[otherIndex, :] = b.reshape(A.shape[1])

		if bsplineRank == 0:
			for otherIndex in range(A.shape[1]):
				Aslice = A[:, otherIndex]
				a = UnpackMatrix(Aslice)
				S = OverlapMatrix[:, otherIndex]
				if norm(a) > 10e-10:
					s = UnpackMatrix(S)
					sInv = numpy.linalg.inv(s)
					
					b = dot(sInv, a)
					B[:,otherIndex] = b.reshape(A.shape[0])

	return B
	
def UnpackMatrix(A):
	rank = len(A.shape)
	if rank == 1:
		N = sqrt(len(A))
		a = A.reshape((N,N))

	if rank == 2:
		N0 = sqrt(A.shape[0])
		N1 = sqrt(A.shape[1])

		a = zeros((N0*N1, N0*N1))
		for i0 in range(N0):
			for j0 in range(N0):
				for i1 in range(N1):
					for j1 in range(N1):
						a[i0*N1 + i1, j0*N1 + j1] = A[i0 + j0*N0, i1 + j1*N1]

	return a


def MatrixMultiply(matrix, inVector, outVector, geometryList):
	rank = len(matrix.shape)

	source = inVector
	dest = outVector
	dest[:] = 0

	if rank == 1:
		pairs0 = geometryList[0].GetBasisPairs()
		for i in xrange(pairs0.shape[0]):
			r0, c0 = pairs0[i, 0], pairs0[i,1]
			dest[r0] += matrix[i] * source[c1]

	if rank == 2:
		pairs0 = geometryList[0].GetBasisPairs()
		pairs1 = geometryList[1].GetBasisPairs()
		
		for i in xrange(pairs0.shape[0]):
			r0, c0 = pairs0[i, 0], pairs0[i,1]
			for j in xrange(pairs1.shape[0]):
				r1, c1 = pairs1[j, 0], pairs1[j,1]
		
				dest[r0, r1] += matrix[i, j] * soure[c0, c1]


def TestBanded():
	conf = pyprop.Load("config.ini")
	conf.KineticEnergy0.geometry0 = "Dense"
	conf.KineticEnergy0.geometry1 = "Dense"
	conf.KineticEnergy1.geometry0 = "Dense"
	conf.KineticEnergy1.geometry1 = "Dense"
	conf.TestPotential.geometry0 = "Dense"
	conf.TestPotential.geometry1 = "Dense"
	propDense = pyprop.Problem(conf)
	propDense.SetupStep()
	tempDense = propDense.GetTempPsi()
	tempDense.GetData()[:] = 0
	propDense.MultiplyHamiltonian(tempDense)

	conf = pyprop.Load("config.ini")
	conf.KineticEnergy0.geometry0 = "Banded"
	conf.KineticEnergy0.geometry1 = "Banded"
	conf.KineticEnergy1.geometry0 = "Banded"
	conf.KineticEnergy1.geometry1 = "Banded"
	conf.TestPotential.geometry0 = "Banded"
	conf.TestPotential.geometry1 = "Banded"
	propBanded = pyprop.Problem(conf)
	propBanded.SetupStep()
	tempBanded = propBanded.GetTempPsi()
	tempBanded.GetData()[:] = 0
	propBanded.MultiplyHamiltonian(tempBanded)

	figure()
	pcolormesh(tempDense.GetData().real.copy())
	figure()
	pcolormesh(tempBanded.GetData().real.copy())


def TestStability():
	conf = pyprop.Load("config.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	for t in prop.Advance(50):
		print "t = %.4f, E(t) = %.6f" % (t, prop.GetEnergyExpectationValue())

	initPsi = prop.psi

	conf = pyprop.Load("config_radial.ini")
	conf.Propagation.timestep = abs(conf.Propagation.timestep)
	conf.Propagation.renormalization = False
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	prop.psi.GetData()[:] = initPsi.GetData()
	for t in prop.Advance(50):
		print "t = %.4f, N(t) = %.6f, P(t) = %.6f" % (t, prop.psi.GetNorm(), abs(prop.psi.InnerProduct(initPsi))**2)

def Propagate(algo=1):
	conf = pyprop.Load("config_radial.ini")
	prop = pyprop.Problem(conf)
	prop.psi.GetRepresentation().Algorithm = algo
	prop.SetupStep()
	for t in prop.Advance(10):
		print "t = %.4f, E(t) = %.6f" % (t, prop.GetEnergyExpectationValue())

	initPsi = prop.psi

	conf = pyprop.Load("config_radial.ini")
	conf.Propagation.timestep = abs(conf.Propagation.timestep)
	conf.Propagation.renormalization = False
	conf.Propagation.grid_potential_list.append("LaserPotential")
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	prop.psi.GetRepresentation().Algorithm = algo

	for i in range(prop.psi.GetData().shape[0]):
		prop.psi.GetData()[i, :] = initPsi.GetData()[i, :]
	prop.psi.Normalize()

	initPsi = prop.psi.Copy()
	
	timeList = []
	corrList = []
	normList = []

	for t in prop.Advance(20):
		n = prop.psi.GetNorm()
		if n != n:
			return prop
		c = abs(prop.psi.InnerProduct(initPsi))**2
		timeList.append(t)
		normList.append(n)
		corrList.append(c)
		print "t = %.4f, N(t) = %.6f, P(t) = %.6f" % (t, n, c)
		hold(False)
		pcolormesh(abs(prop.psi.GetData())**2)
		draw()

	prop.CorrelationList = array(corrList)
	prop.NormList = array(normList)
	prop.TimeList = array(timeList)

	return prop

def LaserFunction(conf, t):
	if t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= cos(t * conf.frequency);
	else:
		curField = 0
	return curField

import time
def TestInnerProduct():
	conf = pyprop.Load("config_radial.ini")
	conf.Propagation.grid_potential_list = []
	prop = pyprop.Problem(conf)
	
	avgCount = 100
	minCount = 10

	tempPsi = prop.GetTempPsi()
	
	prop.psi.GetData()[:] = rand(*prop.psi.GetData().shape)
	tempPsi.GetData()[:] = rand(*prop.psi.GetData().shape)

	for algo in range(1,3):
		prop.psi.GetRepresentation().Algorithm = algo
		n = prop.psi.InnerProduct(tempPsi)
		print "Norm (Algo %i) = %f" % (algo, abs(n)**2)

	for algo in range(1,3):
		prop.psi.GetRepresentation().Algorithm = algo
		minT = 1e10
		for i in range(minCount):
			t = - time.time()
			for j in range(avgCount):
				n = prop.psi.InnerProduct(tempPsi)
			t += time.time()
			if t<minT:
				minT = t / avgCount
		print "Algorithm %i: %f" % (algo, minT)


def TestMatrixMultiply():
	conf = pyprop.Load("config_radial.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	avgCount = 10
	minCount = 10

	tempPsi = prop.GetTempPsi()
	prop.psi.GetData()[:] = rand(*prop.psi.GetData().shape)
	#prop.psi.GetData()[1:,:] = 0

	pot = prop.Propagator.BasePropagator.PotentialList[0]

	algoCount = 6

	for algo in [4,5]: #range(algoCount):
		pot.MultiplyAlgorithm = algo
		tempPsi.GetData()[:] = 0
		pot.MultiplyPotential(tempPsi, 0, 0)
		#print tempPsi.GetData().real
		#figure()
		#pcolormesh(abs(tempPsi.GetData())**2)
		d = tempPsi.InnerProduct(prop.psi)
		print "Overlap (Algo %i) = %f" % (algo, abs(d)**2)

	for algo in [4,5]:
		pot.MultiplyAlgorithm = algo
		minT = 10e10
		for i in range(minCount):
			tempPsi.GetData()[:] = 0
			t = - time.time()
			for j in range(avgCount):
				pot.MultiplyPotential(tempPsi, 0, 0)
				#prop.MultiplyHamiltonian(tempPsi)
			t += time.time()
			if t<minT:
				minT = t / avgCount
		print "Algorithm %i: %f" % (algo, minT)
