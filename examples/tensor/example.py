import sys

sys.path.append("./pyprop")
import pyprop
pyprop = reload(pyprop)
pyprop.ProjectNamespace = globals()

from numpy import array
from numpy import complex
from numpy import zeros

#from pylab import *

execfile("TensorGenerator.py")

from libpotential import *

#------------------------------------------------------------------------------------
#                       Setup Functions
#------------------------------------------------------------------------------------


def SetupConfig(**args):
	configFile = args.get("configFile", "config_radial.ini")
	conf = pyprop.Load(configFile)
	
	if "silent" in args:
		silent = args["silent"]
		conf.Propagation.silent = silent

	if "imtime" in args:
		imtime = args["imtime"]
		if imtime:
			conf.Propagation.timestep = -1.0j * abs(conf.Propagation.timestep)
			conf.Propagation.renormalization = True
		else:
			conf.Propagation.timestep = abs(conf.Propagation.timestep)
			conf.Propagation.renormalization = False

	if "duration" in args:
		duration = args["duration"]
		conf.Propagation.duration = duration

	if "eigenvalueCount" in args:
		conf.Arpack.krylov_eigenvalue_count = args["eigenvalueCount"]

	additionalPotentials = args.get("additionalPotentials", [])
	conf.Propagation.grid_potential_list += additionalPotentials

	return conf


def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	return prop


def FindGroundstate(**args):
	prop = SetupProblem(imtime=True, **args)
	for t in prop.Advance(10):
		E = prop.GetEnergy()
		if pyprop.ProcId == 0:
			print "t = %02.2f, E = %2.8f" % (t, E)

	E = prop.GetEnergyExpectationValue()
	print "t = %02.2f, E = %2.8f" % (t, E)

	return prop


def FindEigenvalues(**args):
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	solver.Solve()
	print solver.GetEigenvalues()
	return solver

#------------------------------------------------------------------------------------
#                       New test-functions
#------------------------------------------------------------------------------------

def GetAutocorrelation1(**args):
	solver = FindEigenvalues(silent=True, **args)
	initPsi = solver.BaseProblem.psi
	solver.SetEigenvector(initPsi, 0)

	prop = SetupProblem(imtime=False, silent=True, **args)
	prop.psi.GetData()[:] = initPsi.GetData()

	def GetAutoCorr(prop):
		return initPsi.InnerProduct(prop.psi)

	t = []
	corr = []
	#c = array([[t, GetAutoCorr(prop)] for t in prop.Advance(True)])
	#t, corr = c[:,0].real, c[:,1]
	for curT in prop.Advance(True):
		t.append(curT)
		corr.append(GetAutoCorr(prop))

	return array(t), array(corr), prop

def GetAutocorrelation2(**args):
	solver = FindEigenvalues(silent=True, **args)
	initPsi = solver.BaseProblem.psi

	prop = SetupProblem(imtime=False, silent=True, **args)
	
	prop.psi.GetData()[:] = 0
	solver.SetEigenvector(initPsi, 0)
	prop.psi.GetData()[:] += initPsi.GetData() / sqrt(2.)
	solver.SetEigenvector(initPsi, 1)
	prop.psi.GetData()[:] += initPsi.GetData() / sqrt(2.)
	prop.psi.Normalize()

	initPsi.GetData()[:] = prop.psi.GetData()

	def GetAutoCorr(prop):
		return initPsi.InnerProduct(prop.psi)

	t = []
	corr = []
	#c = array([[t, GetAutoCorr(prop)] for t in prop.Advance(True)])
	#t, corr = c[:,0].real, c[:,1]
	for curT in prop.Advance(True):
		t.append(curT)
		corr.append(GetAutoCorr(prop))

	return array(t), array(corr), prop


def FindIonizationProbability(**args):
	solver = FindEigenvalues(silent=True, **args)
	initPsi = solver.BaseProblem.psi

	prop = SetupProblem(imtime=False, silent=True, additionalPotentials=["LaserPotentialLength"], **args)
	
	solver.SetEigenvector(initPsi, 0)
	prop.psi.GetData()[:] = initPsi.GetData() 

	def GetAutoCorr(prop):
		return initPsi.InnerProduct(prop.psi)

	for t in prop.Advance(10):
		corr = abs(GetAutoCorr(prop))**2
		print "t = %2.2f, c(t) = %f" % (t, corr)

	corr = abs(GetAutoCorr(prop))**2
	print "t = %2.2f, c(t) = %f" % (t, corr)




#------------------------------------------------------------------------------------
#                       Old test-functions
#------------------------------------------------------------------------------------

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

import pypar

def Propagate(initPsi, algo=1, **args):
	#prop = FindGroundstate(silent=True)
	#initPsi = prop.psi

	#prop = SetupProblem(imtime=False, **args)
	prop = SetupProblem(silent=False, imtime=False, additionalPotentials=["LaserPotentialVelocityDerivativeR1", "LaserPotentialVelocityDerivativeR2", "LaserPotentialVelocity"], **args)

	prop.psi.Clear()
	prop.psi.GetData()[:,:,:initPsi.GetData().shape[2]] = initPsi.GetData()
	prop.psi.Normalize()

	timeList = []
	corrList = []
	normList = []
	for t in prop.Advance(10):
		n = prop.psi.GetNorm()
		if n != n:
			return prop
		c = abs(prop.psi.InnerProduct(initPsi))**2
		timeList.append(t)
		normList.append(n)
		corrList.append(c)
		if pyprop.ProcId == 0:
			print "t = %.4f, N(t) = %.6f, P(t) = %.6f" % (t, n, c)
		#hold(False)
		#pcolormesh(abs(prop.psi.GetData())**2)
		#draw()
		#sys.stdout.flush()
		
	
	prop.Propagator.PampWrapper.PrintStatistics()
		
	c = abs(prop.psi.InnerProduct(initPsi))**2
	print "Final Correlation = %f" % c

	prop.CorrelationList = array(corrList)
	prop.NormList = array(normList)
	prop.TimeList = array(timeList)

	return prop


def SetPotential2(potentialData):
	N = sqrt(potentialData.shape[0])
	radialN = potentialData.shape[1]
	pot = potentialData.reshape(N, N, radialN)

	for radialIndex in range(radialN):
		for i in range(1,N):
			j = i - 1
			pot[i,j, radialIndex] = - 1.0j * i * pot[i,j, radialIndex]
			pot[j,i, radialIndex] =   1.0j * i * pot[j,i, radialIndex]
				

def LaserFunctionVelocity(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= - cos(t * conf.frequency);
	else:
		curField = 0
	return curField


def LaserFunctionSimpleLength(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= cos(t * conf.frequency);
	else:
		curField = 0
	return curField



def LaserFunctionLength(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		T = conf.pulse_duration
		w = conf.frequency
		curField *= sin(pi*t/T)*(-2*pi/T * cos(pi*t/T) * cos(w*t) + w*sin(pi*t/T)*sin(w*t))
	else:
		curField = 0
	return curField


import time
def TestInnerProduct():
	conf = pyprop.Load("config_radial.ini")
	conf.Propagation.grid_potential_list = []
	prop = pyprop.Problem(conf)
	
	avgCount = 10
	minCount = 10

	tempPsi = prop.GetTempPsi()
	tempPsi2 = prop.psi.Copy()
	tempPsi3 = prop.psi.Copy()
	
	prop.psi.GetData()[:] = rand(*prop.psi.GetData().shape)
	tempPsi.GetData()[:] = rand(*prop.psi.GetData().shape)

	"""

	generator = prop.Propagator.BasePropagator.TensorPotentialGenerator

	#Use TensorPotentialGenerator to construct potential in basis
	geometryList = generator.GetGeometryList(conf.InnerProductPotential1)
	potentialData = generator.GeneratePotential(conf.InnerProductPotential1)
	innerProductPotential1 = TensorPotential(prop.psi)
	conf.InnerProductPotential1.Apply(innerProductPotential1)
	innerProductPotential1.GeometryList = geometryList
	innerProductPotential1.PotentialData = potentialData
	innerProductPotential1.Name = "InnerProdPotential1"
	innerProductPotential1.SetupStep(0)

	#Use TensorPotentialGenerator to construct potential in basis
	geometryList = generator.GetGeometryList(conf.InnerProductPotential2)
	potentialData = generator.GeneratePotential(conf.InnerProductPotential2)
	innerProductPotential2 = TensorPotential(tempPsi2)
	conf.InnerProductPotential2.Apply(innerProductPotential2)
	innerProductPotential2.GeometryList = geometryList
	innerProductPotential2.PotentialData = potentialData
	innerProductPotential2.Name = "InnerProdPotential2"
	innerProductPotential2.SetupStep(0)

	#Use TensorPotentialGenerator to construct potential in basis
	geometryList = generator.GetGeometryList(conf.InnerProductPotential3)
	potentialData = generator.GeneratePotential(conf.InnerProductPotential3)
	innerProductPotential3 = TensorPotential(tempPsi3)
	conf.InnerProductPotential3.Apply(innerProductPotential3)
	innerProductPotential3.GeometryList = geometryList
	innerProductPotential3.PotentialData = potentialData
	innerProductPotential3.Name = "InnerProdPotential3"
	innerProductPotential3.SetupStep(0)

	repr = prop.psi.GetRepresentation()
	overlapMatrix0 = repr.GetGlobalOverlapMatrixBlas(0)
	overlapMatrix1 = repr.GetGlobalOverlapMatrixBlas(1)
	overlapMatrix2 = repr.GetGlobalOverlapMatrixBlas(2)

	def NewInnerProduct():
		tempPsi2.GetData()[:] = 0
		tempPsi3.GetData()[:] = 0
		innerProductPotential1.MultiplyPotential(tempPsi2, 0, 0)
		innerProductPotential2.MultiplyPotential(tempPsi3, 0, 0)
		tempPsi2.GetData()[:] = 0
		innerProductPotential3.MultiplyPotential(tempPsi2, 0, 0)
		return VectorDotProduct(tempPsi, tempPsi2)

	def NewInnerProduct2():
		tempPsi2.GetData()[:] = 0
		tempPsi3.GetData()[:] = 0
		repr.MultiplyOverlapMatrix(prop.psi, tempPsi2, 0)
		repr.MultiplyOverlapMatrix(tempPsi2, tempPsi3, 1)
		repr.MultiplyOverlapMatrix(tempPsi3, tempPsi2, 2)
		return VectorDotProduct(tempPsi, tempPsi2)

	def NewInnerProduct3():
		tempPsi2.GetData()[:] = 0
		tempPsi3.GetData()[:] = 0

		N1, N2, N3 = tempPsi.GetData().shape

		source = prop.psi.GetData().reshape(1, N1, N2*N3)
		dest = tempPsi2.GetData().reshape(1, N1, N2*N3)
		MultiplyOverlapMatrix(overlapMatrix0, source, dest)

		source = tempPsi2.GetData().reshape(N1, N2, N3)
		dest = tempPsi3.GetData().reshape(N1, N1, N3)
		MultiplyOverlapMatrix(overlapMatrix1, source, dest)

		source = tempPsi3.GetData().reshape(N1*N2, N3, 1)
		dest = tempPsi2.GetData().reshape(N1*N2, N3, 1)
		MultiplyOverlapMatrix(overlapMatrix2, source, dest)

		return VectorDotProduct(tempPsi, tempPsi2)
	"""


	for algo in range(1,4):
		prop.psi.GetRepresentation().Algorithm = algo
		if algo == 4:
			innerProduct = NewInnerProduct
		elif algo == 5:
			innerProduct = NewInnerProduct2
		elif algo == 6:
			innerProduct = NewInnerProduct3
		else:
			innerProduct = lambda: prop.psi.InnerProduct(tempPsi)

		prop.psi.GetRepresentation().Algorithm = algo
		n = innerProduct()
		print "Norm (Algo %i) = %f" % (algo, abs(n)**2)

	for algo in range(1,4):
		prop.psi.GetRepresentation().Algorithm = algo
		if algo == 4:
			innerProduct = NewInnerProduct
		elif algo == 5:
			innerProduct = NewInnerProduct2
		elif algo == 6:
			innerProduct = NewInnerProduct3
		else:
			innerProduct = lambda: prop.psi.InnerProduct(tempPsi)

		minT = 1e10
		for i in range(minCount):
			t = - time.time()
			for j in range(avgCount):
				n = innerProduct()
			t += time.time()
			if t<minT:
				minT = t / avgCount
		print "Algorithm %i: %f" % (algo, minT)



def TestMatrixMultiply():
	def SetupProblem(geometry0, geometry1):
		conf = pyprop.Load("config.ini")
		conf.Propagation.silent = True
		conf.KineticEnergy0.geometry0 = geometry0
		conf.KineticEnergy0.geometry1 = geometry1
		conf.KineticEnergy1.geometry0 = geometry0
		conf.KineticEnergy1.geometry1 = geometry1
		conf.TestPotential.geometry0 = geometry0
		conf.TestPotential.geometry1 = geometry1
		prop = pyprop.Problem(conf)
		prop.SetupStep()
		
		return prop

	algoList = [("Banded", "Banded"), ("Banded-Old", "Banded-Old"), ("Banded-Old", "Banded"), ("Hermitian", "Hermitian"), ("Dense", "Dense")]

	avgCount = 10
	minCount = 10

	baseProp = SetupProblem(*algoList[0])
	baseProp.psi.GetData()[:] = rand(*baseProp.psi.GetData().shape)

	for geom1, geom2 in algoList: 
		prop = SetupProblem(geom1, geom2)
		prop.psi.GetData()[:] = baseProp.psi.GetData()
		tempPsi = prop.GetTempPsi()
		tempPsi.GetData()[:] = 0
		pot = prop.Propagator.BasePropagator.PotentialList[0]
		#pot.MultiplyPotential(tempPsi, 0, 0)
		prop.MultiplyHamiltonian(tempPsi)
		d = tempPsi.InnerProduct(prop.psi)
		print "Overlap (Algo %i) = %f" % (algoList.index((geom1, geom2)), abs(d)**2)

	for geom1, geom2 in algoList: 
		prop = SetupProblem(geom1, geom2)
		print prop.Propagator.BasePropagator.SolveAlgorithm
		prop.psi.GetData()[:] = baseProp.psi.GetData()
		tempPsi = prop.GetTempPsi()
		tempPsi.GetData()[:] = 0
		pot = prop.Propagator.BasePropagator.PotentialList[0]
		
		minT = 10e10
		for i in range(minCount):
			tempPsi.GetData()[:] = 0
			t = - time.time()
			for j in range(avgCount):
				#pot.MultiplyPotential(tempPsi, 0, 0)
				prop.MultiplyHamiltonian(tempPsi)
			t += time.time()
			if t<minT:
				minT = t / avgCount
		print "Time (Algo %i) = %f" % (algoList.index((geom1, geom2)), minT)


def TestMatrixMultiply2():
	def SetupProblem(geometry0, geometry1):
		conf = pyprop.Load("config_radial.ini")
		conf.Propagation.silent = True
		conf.Propagation.grid_potential_list = ["LaserPotential"]
		conf.LaserPotential.geometry0 = geometry0
		conf.LaserPotential.geometry1 = geometry1
		prop = pyprop.Problem(conf)
		prop.SetupStep()
		
		return prop

	algoList = [("BandedDistributed", "Banded"), ("DipoleSelectionRule", "Banded")]

	avgCount = 10
	minCount = 10

	baseProp = SetupProblem(*algoList[0])
	baseProp.psi.GetData()[:] = rand(*baseProp.psi.GetData().shape)

	for geom1, geom2 in algoList: 
		prop = SetupProblem(geom1, geom2)
		prop.psi.GetData()[:] = baseProp.psi.GetData()
		tempPsi = prop.GetTempPsi()
		tempPsi.GetData()[:] = 0
		pot = prop.Propagator.BasePropagator.PotentialList[0]
		pot.MultiplyPotential(tempPsi, 0, 0)
		print pot.MultiplyFunction.__name__
		d = tempPsi.InnerProduct(prop.psi)
		d = tempPsi.InnerProduct(tempPsi)
		print "Overlap (Algo %i) = %f" % (algoList.index((geom1, geom2)), abs(d)**2)

	for geom1, geom2 in algoList: 
		prop = SetupProblem(geom1, geom2)
		prop.psi.GetData()[:] = baseProp.psi.GetData()
		tempPsi = prop.GetTempPsi()
		tempPsi.GetData()[:] = 0
		pot = prop.Propagator.BasePropagator.PotentialList[0]
		
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
		print "Time (Algo %i) = %f" % (algoList.index((geom1, geom2)), minT)
