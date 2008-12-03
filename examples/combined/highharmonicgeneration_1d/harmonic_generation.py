#Import system modules
import time
import sys
import os
from pylab import *
from numpy import *

#Load pyprop
sys.path.insert(1, os.path.abspath("./pyprop"))
import pyprop
pyprop = reload(pyprop)

#Load the project module
from libpotential import *

#Load constants
execfile("constants.py")

pyprop.ProjectNamespace = globals()


def SetupConfig(**args):
	#Decide which config file to use
	configFile = "config.ini"
	if "config" in args:
		configFile = args["config"]

	#Load the config file
	conf = pyprop.Load(configFile)

	#Modify the config
	if "imtime" in args:
		imtime = args["imtime"]
		propSection = conf.Propagation
		dt = abs(propSection.timestep)
		renormalize = False
		if imtime:
			print "yes"
			dt = -1.0j * dt
			renormalize = True

		propSection.timestep = dt
		propSection.renormalization = renormalize

	if "amplitude" in args:
		amplitude = args["amplitude"]
		conf.DynamicPotential.amplitude = amplitude

	
	if "xmax" in args:
		xMax = args['xmax']
		conf.BSplineRepresentation.xmax = xMax
	
	if "xsize" in args:
		xSize = args['xsize']
		conf.BSplineRepresentation.xsize = xSize
	
	if "dt" in args:
		timeStep = args['dt']
		conf.Propagation.timestep = timeStep
	
	if "duration" in args:
		duration = args['duration']
		conf.Propagation.duration = duration
		
	return conf


def SetupProblem(**args):
	"""
	Load configuration and set up problem
	"""
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop

def SetupInitialState(**args):
	
	stateIndex = args["stateIndex"]

	#Set up problem
	conf = SetupConfig(**args)

	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#Find eigenstates of desired l-subspace
	M = GetHamiltonMatrix(prop)
	print "Finding eigenvectors and eigenvalues..."
	sys.stdout.flush()
	E,V = eig(M)
	I = argsort(E)
	print "Initial state energy = ", E[I[stateIndex]].real
	sys.stdout.flush()

	#Assign initial state
	prop.psi.GetData()[:] = V[:,I[stateIndex]]
	prop.psi.Normalize()

	return prop

def FindGroundstate(**args):
	"""
	Find groundstate using imaginary time propagation.
	"""

	args['imtime'] = True
	prop = SetupProblem(**args)
	
	for t in prop.Advance(10):
		E = prop.GetEnergy()
		print "t = %f, E = %f" % (t, E)

	E = prop.GetEnergy()
	print "Ground State Energy = %f" % E

	#prop.SaveWavefunctionHDF("groundstate.h5", "/wavefunction")

	return prop


def GetHamiltonMatrix(prop):
	"""
	Set up the hamilton matrix of the specified problem, using
	MultiplyHamiltonian().
	"""

	size = prop.psi.GetData().shape[0]
	matrix = zeros((size, size), dtype=complex)
	tempPsi = prop.GetTempPsi()

	for i in range(size):
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[i] = 1

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		
		matrix[:, i] = tempPsi.GetData()[:]
		
	return matrix


def EstimateRunTime(prop, returnTime = False, repetitions = 1, sampleSize = 20):
	"""
	Estime total propagation time by doing a few initial steps
	"""

	#Print some time info
	startTime = time.time()
	for i in range(sampleSize):
		prop.AdvanceStep()
	finishTime = time.time()
	numberOfSteps = int(round(prop.Duration / abs(prop.TimeStep)))
	stepTime = (finishTime - startTime)
	totalTime = int(round(stepTime * numberOfSteps / sampleSize)) * repetitions
	totalTimefmt = time.gmtime(totalTime)
	if pyprop.ProcId == 0 and not returnTime:
		print 
		print "    Propagation time is: %.2f fs. (%f a.u.)" % \
			(prop.Duration / femtosec_to_au, prop.Duration)
		print "    Very approximate expected runtime: %i days, %i hours, %i minutes" \
			% (totalTimefmt[2]-1, totalTimefmt[3], totalTimefmt[4])
		print 
		sys.stdout.flush()

	if returnTime:
		return totalTimefmt

def GetGridLinear(conf, xmax=None, xmin=None):
	if xmin == None: xmin = conf.xmin
	if xmax == None: xmax = conf.xmax
	count = conf.count

	start = xmin
	end = xmax

	if not conf.include_left_boundary:
		count += 1
	if not conf.include_right_boundary:
		count += 1

	dx = (xmax - xmin) / float(count-1)
	if not conf.include_left_boundary:
		start += dx
	if conf.include_right_boundary:
		end += dx
	grid = r_[start:end:dx]

	return array(grid, dtype=double)


def SetupPotential(conf, rank=1):
	potential = eval(conf.classname + "_%s()" % rank)
	potential.ApplyConfigSection(conf)
	return potential


def PropagateHHG(psi, config = "highharmonicgeneration.ini", numSteps = 3000, sampleSize=20):
	#Set up propagator
	prop = SetupProblem(config = config, stateIndex = 0)
	prop.psi.GetData()[:] = psi.GetData()[:]

	#Set up gradient potential
	prop.dipoleAcc = []
	gradPotentialConf = prop.Config.DiffPotential
	gradPotential = SetupPotential(gradPotentialConf)

	#Set up dipole potential
	prop.dipoleMoment = []
	dipolePotConf = prop.Config.StarkPotential
	dipolePotConf.field_strength = 1.0
	dipolePot = SetupPotential(dipolePotConf)

	#Set up tmpPsi
	tmpPsi = prop.psi.CopyDeep()

	#Store initial wavefunction
	initPsi = prop.psi.CopyDeep()
	prop.initialCorr = []
	prop.norm = []

	#Estimate runtime
	EstimateRunTime(prop, sampleSize=sampleSize)

	#Propagate
	prop.sampleTimes = []
	infoStr = ""
	for i, t in enumerate(prop.Advance(numSteps)):
		infoStr = "\b" * len(infoStr)
		infoStr +=  "Progress: %04i/%04i" % (i, numSteps)
		sys.stdout.write(infoStr)
		sys.stdout.flush()

		prop.sampleTimes.append(t)

		#Calculate norm and initial correlation
		prop.norm.append(prop.psi.GetNorm())
		prop.initialCorr.append(abs(prop.psi.InnerProduct(initPsi))**2)

		#Calculate expectation value <grad(V)>
		gradPotExpVal = \
			prop.Propagator.CalculatePotentialExpectationValue(tmpPsi, gradPotential, 0, 0)
		prop.dipoleAcc.append(gradPotExpVal)
		
		#Calculate expectation value <z>
		dipolePotExpVal = \
			prop.Propagator.CalculatePotentialExpectationValue(tmpPsi, dipolePot, 0, 0)
		prop.dipoleMoment.append(dipolePotExpVal)
	
	return prop
