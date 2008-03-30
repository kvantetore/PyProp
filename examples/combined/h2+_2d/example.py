#improt system modules
import sys
import os

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("../../.."))

#Load and reload pyprop in order to get recent changes
import pyprop
pyprop = reload(pyprop)
from libh2p import *

#numpy an pylab for good measure
from pylab import *
from numpy import *

#Choose radial grid type:
class RadialGridType:
	CARTESIAN = 1
	TRANSFORMED = 2
	ORTHOPOL = 3

def SetupProblem(**args):
	#load config file. hydrogen.ini uses ../sphericalbase.ini as a base
	#configuration file, so be sure to check out that one as well.
	conf = SetupConfig(**args)

	#Uses the same Problem class, only changes the propagator specified in the
	#config file
	prop = pyprop.Problem(conf)

	#Set up all transformations and potentials.
	prop.SetupStep()

	return prop

def FindEigenvalues(**args):
	prop = SetupProblem(**args)
	tempPsi = prop.GetTempPsi()

	solver = pyprop.PiramSolver(prop)
	solver.Solve()
	print solver.Solver.GetEigenvalues().real

	#test if the find eigenvectors are really eigenvectors
	for i in range(len(solver.Solver.GetEigenvalues())):
		solver.SetEigenvector(prop.psi, i)

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		n1 = tempPsi.GetNorm()
		n2 = prop.psi.GetNorm()
		lindep = prop.psi.InnerProduct(tempPsi) / (n1 * n2)

		E = prop.GetEnergyExpectationValue()
		
		print "E = %f, cos(theta) = %f" % (E, abs(lindep))

		#Propagate the eigenstate to see if it remains in the same state
		#args['silent'] = True
		#args['initPsi'] = prop.psi
		#Propagate(**args)

	return solver

def FindGroundstate(**args):
	#for finding ground state, we always use imaginary time
	args['imTime'] = True  #Remember: Arguments are case sensitive

	#Create the propagator
	prop = SetupProblem(**args)

	#propagate
	silent = False
	if 'silent' in args: silent = args['silent']
	for t in prop.Advance(10):
		if not silent:
			print "t = ", t, ", E = ", prop.GetEnergyExpectationValue()

	return prop

def Propagate(**args):
	args['imTime'] = False
	prop = SetupProblem(**args)

	if 'initPsi' in args:
		initPsi = args['initPsi']
	else:
		initPsi = prop.psi.Copy()

	initPsi.Normalize()
	prop.psi.GetData()[:] = initPsi.GetData()

	t = 0
	norm = prop.psi.GetNorm()
	corr = abs(prop.psi.InnerProduct(initPsi))**2
	print "t = %f, N(t) = %f, C(t) = %f" % (t, norm, corr)


	figure()
	if prop.psi.GetRank() == 2:
		pyprop.Plot2DRank(prop, 0)
	else:
		pyprop.Plot1D(prop)

	for t in prop.Advance(10):
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		print "t = %f, N(t) = %f, C(t) = %f" % (t, norm, corr)
		if prop.psi.GetRank() == 2:
			pyprop.Plot2DRank(prop, 0)
		else:
			pyprop.Plot1D(prop)

def FindEigenvaluesImtime(**args):
	eigenvalueCount = args['eigenvalueCount']
	args['imTime'] = True
	args['silent'] = True

	eigenstates = []

	for curEigenvalue in range(eigenvalueCount):
		prop = SetupProblem(**args)
		for t in prop.Advance(True):
			#Remove projection on previous states
			for state in eigenstates:
				overlap = prop.psi.InnerProduct(state)
				prop.psi.GetData()[:] -= overlap * state.GetData()
			prop.psi.Normalize()

			for state in eigenstates:
				newOverlap = prop.psi.InnerProduct(state)
				if abs(newOverlap) > 10e-15:
					print "NewOverlap = %g" % abs(newOverlap)

		E = prop.GetEnergy()

		#Remove projection on previous states
		for state in eigenstates:
			overlap = prop.psi.InnerProduct(state)
			prop.psi.GetData()[:] -= overlap * state.GetData()

		#Check whether H |psi> is linearly dependent on |psi>
		tempPsi = prop.GetTempPsi()
		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		n1 = tempPsi.GetNorm()
		n2 = prop.psi.GetNorm()
		lindep = prop.psi.InnerProduct(tempPsi) / (n1 * n2)

		print "E = %f, cos(theta) = %f" % (E, abs(lindep))
		eigenstates.append(prop.psi.Copy())
		del prop

	
	for state in eigenstates:
		args['initPsi'] = state
		Propagate(**args)


import time
def TestSphericalHarmonicTransform(lmax=32, gridSize=256):
	prop = SetupProblem(gridType=RadialGridType.CARTESIAN, gridSize=gridSize, lmax=lmax)
	prop.psi.GetData()[:] = rand(*prop.psi.GetData().shape)

	sphProp = prop.Propagator.SubPropagators[1]

	avgCount = 10
	minCount = 10

	for algo in range(-1,5):
		sphProp.Transform.transform.Algorithm = algo
		initData = prop.psi.GetData().copy()
		sphProp.Transform.InverseTransform(prop.psi)
		sphProp.Transform.ForwardTransform(prop.psi)
		print "Error (algo %i) = %f" % (algo, sum(abs(initData - prop.psi.GetData())**2))
	
	

	for algo in range(-1,5):
		sphProp.Transform.transform.Algorithm = algo
		minT = 10e10
		for i in range(minCount):
			t = - time.time()
			for j in range(avgCount):
				sphProp.Transform.InverseTransform(prop.psi)
				sphProp.Transform.ForwardTransform(prop.psi)
			t += time.time()
			if t<minT:
				minT = t / avgCount
		print "Algorithm %i: %f" % (algo, minT)

#---------------------------------------------------------------------------------

def SetupConfig(**args):
	"""
	Loads a config file, and modifies it with the arguments sendt to this function 
	i.e. conf = SetupConfig(configFile='h2p.ini', dt=0.04, imTime=True), 
	will return a conf object loading h2p.ini, setting timestep to -0.04mj, and 
	renormalization to True.
	"""
	silent = False
	if 'silent' in args:
		silent = args['silent']
	
	#config file
	if 'configFile' in args:
		configFile = args['config']
	else:
		configFile = "h2p.ini"
	conf = pyprop.Load(configFile)
	
	#Grid size
	if 'gridSize' in args:
		gridSize = args['gridSize']
		if not silent: print "Using Grid Size = ", gridSize
		#set grid size of both representations, but we will
		#only use one of them in a given run
		conf.CartesianRadialRepresentation.rank0[2] = gridSize
		conf.TransformedRadialRepresentation.n = gridSize
		conf.OrthoPolRadialRepresentation.n = gridSize

	#Grid type
	if 'gridType' in args:
		gridType = args['gridType']
	else:
		gridType = None
	SetRadialGridType(conf, gridType)
	
	#TimeStep
	if 'dt' in args:
		dt = args['dt']
		if not silent: print "Using TimeStep = ", dt
		conf.Propagation.timestep = abs(dt)

	#Imaginary Time Propagation
	if 'imTime' in args:
		imTime = args['imTime']
		if not silent: print "Using ImaginaryTime = ", imTime
		dt = conf.Propagation.timestep
		if imTime:
			conf.Propagation.renormalization = True
			conf.Propagation.timestep = -1.0j * abs(dt)
		else:
			conf.Propagation.renormalization = False
			conf.Propagation.timestep = abs(dt)

	if 'softing' in args:
		softing = args['softing']
		if not silent: print "Using softing ", softing
		conf.Potential.softing = softing

	if 'lmax' in args:
		lmax = args['lmax']
		if not silent: print "Using LMax = ", lmax
		conf.AngularRepresentation.maxl = lmax

	if 'duration' in args:
		duration = args['duration']
		conf.Propagation.duration = duration

	if 'orientation' in args:
		orientation = args['orientation']
		conf.Potential.nuclear_orientation = orientation
	
	if 'separation' in args:
		separation = args['separation']
		conf.Potential.nuclear_separation = separation

	if 'silent' in args:
		silent = args['silent']
		conf.Propagation.silent = silent
		
	print "Setup Config Complete"
	print ""
	return conf

def SetRadialGridType(conf, gridType):
	"""
	Chooses between Cartesian and Transformed grid. On a loaded config
	object.
	"""

	if gridType == None:
		gridType = RadialGridType.TRANSFORMED

	#Set RadialRepresentation section
	if gridType == RadialGridType.CARTESIAN:
		conf.RadialRepresentation = conf.CartesianRadialRepresentation
		conf.RadialPropagator.propagator = pyprop.CartesianRadialPropagator

	elif gridType == RadialGridType.TRANSFORMED:
		conf.RadialRepresentation = conf.TransformedRadialRepresentation
		conf.RadialPropagator.propagator = pyprop.TransformedRadialPropagator

	elif gridType == RadialGridType.ORTHOPOL:
		conf.RadialRepresentation = conf.OrthoPolRadialRepresentation
		conf.RadialPropagator.propagator = pyprop.OrthoPolRadialPropagator

	else:
		raise "Invalid grid type ", gridType

	#Set radialtype 
	conf.Representation.radialtype = conf.RadialRepresentation.type


#Investigate the Condition Numbers of the different discretizations
def Condition(matrix):
	U, S, V = linalg.svd(matrix)
	return S[0] / S[-1]

def ConditionScaling(args, variableName, gridSizes):
	args['silent'] = True
	conditionNumbers = []

	for gridSize in gridSizes:
		args[variableName] = gridSize
		prop = SetupProblem(**args)
		M = pyprop.GetBasisExpansionMatrix(prop)
		conditionNumber = Condition(M)
		print "%s = %i, cond = %f" % (variableName, gridSize, conditionNumber)
		conditionNumbers.append(conditionNumber)

	return array(conditionNumbers)

def ConditionScalingTransformed(**args):
	if 'gridType' not in args:
		args['gridType'] = RadialGridType.TRANSFORMED
	args['lmax'] = 10

	gridSizes = [2, 4, 8, 16, 32, 64]
	conditionNumbers = ConditionScaling(args, "gridSize", gridSizes)

	return gridSizes, conditionNumbers

def ConditionScalingAngular(**args):
	args['gridType'] = RadialGridType.TRANSFORMED
	args['gridSize'] = 20

	gridSizes = [2, 4, 8, 16, 32, 64]
	conditionNumbers = ConditionScaling(args, "lmax", gridSizes)

	return gridSizes, conditionNumbers
