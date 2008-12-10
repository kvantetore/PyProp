#Import system modules
import time
import sys
import os
import tables
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


def FindEigenstates(useARPACK = False, **args):
	"""
	Uses pIRAM to find the the lowest eigenvectors of the
	problem specified in **args
	"""
	prop = SetupProblem(**args)

	#use custom initial residual if provided
	initialResidual = args.get("initialResidual")
	if initialResidual != None:
		prop.psi.GetData()[:] = initialResidual.GetData()

	#find eigenstates
	solver = None
	if useARPACK:
		solver = pyprop.ArpackSolver(prop)
	else:
		solver = pyprop.PiramSolver(prop)
	solver.Solve()
	return solver


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


def PropagateHHG(psi, config = "highharmonicgeneration.ini", numSteps = 3000, sampleSize=20, \
	useTensor=False, storeResult=False):

	def SetupPotentialFunction(conf, psi):
		if useTensor:
			pot = SetupTensorPotential(conf, psi)
			pot.SetupStep(timeStep)
			return pot
		else:
			return SetupPotential(conf)

	#Set up propagator
	prop = SetupProblem(config = config, stateIndex = 0)
	prop.psi.GetData()[:] = psi.GetData()[:]
	timeStep = prop.Config.Propagation.timestep

	#Set up gradient potential
	prop.dipoleAcc = []
	gradPotentialConf = prop.Config.DiffPotential
	gradPotential = SetupPotentialFunction(gradPotentialConf, prop.psi)

	#Set up dipole potential
	prop.dipoleMoment = []
	dipolePotConf = prop.Config.StarkPotential
	dipolePotConf.field_strength = 1.0
	dipolePot = SetupPotentialFunction(dipolePotConf, prop.psi)

	#Set up tmpPsi
	tmpPsi = prop.psi.CopyDeep()

	#Store initial wavefunction
	initPsi = prop.psi.CopyDeep()
	prop.initialCorr = []
	prop.norm = []

	#Store initial psi
	outputFile = ""
	if storeResult:
		outputFile = "out/data"
		for el in time.localtime():
			outputFile += "_%s" % el
		outputFile += ".h5"

		prop.SaveWavefunctionHDF(outputFile, "/initial_wavefunction")

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

	#Store result
	if storeResult:
		
		prop.SaveWavefunctionHDF(outputFile, "/wavefunction")
		h5file = tables.openFile(outputFile, "r+")
		try:
			h5file.createArray("/", "SampleTimes", prop.sampleTimes) 
			h5file.createArray("/", "Norm", prop.norm)
			h5file.createArray("/", "InitialCorrelation", prop.initialCorr)
			h5file.createArray("/", "DipoleMoment", prop.dipoleMoment)
			h5file.createArray("/", "DipoleAcceleration", prop.dipoleAcc)
		finally:
			h5file.close()
	
	return prop


def SetupProblemFromFile(file, nodeName=None):
	"""
	Set up problem object and load wavefunction from file.
	"""
	prop = None
	cfgObj = pyprop.serialization.GetConfigFromHDF5(file)
	cfgObj.set("InitialCondition", "type", "None")
	conf = pyprop.Config(cfgObj)
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	GetWavefunctionFromFile(file, prop.psi, nodeName=nodeName)
	
	return prop


def GetWavefunctionFromFile(file, psi, nodeName=None):
	h5file = tables.openFile(file, "r")
	try:
		if nodeName == None:
			for node in h5file.walkNodes():
				if node._v_name == "wavefunction":
					psi.GetData()[:] = node[:]
		else:
			psi.GetData()[:] = h5file.getNode(nodeName)[:]
	finally:
		h5file.close()

#--------------------------------------------------------------------------------------------------
# TENSOR POTENTIAL
#--------------------------------------------------------------------------------------------------
execfile("tensor/TensorGenerator.py")

from tensor.libpotential import *


def Propagate(algo=1, config="config/tensorpotential.ini", **args):
	prop = SetupProblem(silent=True, additionalPotentials=["LaserPotentialVelocity1", "LaserPotentialVelocity2", "LaserPotentialVelocity3"], config=config, **args)
	#prop = SetupProblem(additionalPotentials=["LaserPotentialLength"], config="config/tensorpotential.ini")
	#prop = SetupProblem(silent=True, config="config/tensorpotential.ini")

	initPsi = prop.psi.Copy()

	#pcolormesh(abs(prop.psi.GetData())**2)

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
		if pyprop.ProcId == 0:
			print "t = %.4f, N(t) = %.6f, P(t) = %.6f" % (t, n, c)
		#hold(False)
		#pcolormesh(abs(prop.psi.GetData())**2)
		#draw()
		#sys.stdout.flush()
		#pcolormesh(abs(prop.psi.GetData())**2)
	
	prop.Propagator.PampWrapper.PrintStatistics()
		
	c = abs(prop.psi.InnerProduct(initPsi))**2
	print "Final Correlation = %f" % c

	prop.CorrelationList = array(corrList)
	prop.NormList = array(normList)
	prop.TimeList = array(timeList)

	return prop

def SetupTensorPotential(configSection, psi):
	"""
	Generate a TensorPotential.
	"""

	#Potentials we should create on the fly
	generator = TensorPotentialGenerator(representation = psi.GetRepresentation())

	#Use TensorPotentialGenerator to construct potential in basis
	geometryList = generator.GetGeometryList(configSection)
	potentialData = generator.GeneratePotential(configSection)

	#Create PotentialWrapper for TensorPotential
	potential = TensorPotential(psi)
	configSection.Apply(potential)
	potential.GeometryList = geometryList
	potential.PotentialData = potentialData
	potential.Name = configSection.classname

	return potential

#--------------------------------------------------------------------------------------------------
# Time functions for time-dependent potentials
#--------------------------------------------------------------------------------------------------
def LaserFunctionVelocity(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= - cos(t * conf.frequency) / conf.frequency;
	else:
		curField = 0
	return curField

def LaserFunctionVelocityFromSinSqrLength(conf, t):
	w = conf.frequency
	T = conf.pulse_duration
	if 0 <= t < T:
		curField = 1/4. * conf.amplitude
		curField *= -(2*sin(w*t)/w + T*sin((w-2*pi/T)*t)/(2*pi-T*w) - \
			T*sin((w+2*pi/T)*t)/(2*pi + T*w))
	else:
		curField = 0
	return curField

def LaserFunctionSimpleLength(conf, t):
	pulseStart = getattr(conf, "pulse_start", 0)
	if pulseStart <= t < pulseStart + conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= cos(t * conf.frequency);
	else:
		curField = 0
	return curField


def LaserFunctionLength(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude / conf.frequency;
		T = conf.pulse_duration
		w = conf.frequency
		curField *= sin(pi*t/T)*(-2*pi/T * cos(pi*t/T) * cos(w*t) + w*sin(pi*t/T)*sin(w*t))
	else:
		curField = 0
	return curField


def LaserFunctionLengthFlatTop(conf, t):
	pulseStart = 0
	if conf.Exists("pulse_start"):
		pulseStart = conf.pulse_start
	curField = conf.amplitude * cos(t * conf.frequency);

	if (t > conf.pulse_duration) or (t < pulseStart):
		curField = 0
	elif 0 <= t < conf.ramp_on_time:
		curField *= sin(t * pi / (2*conf.ramp_on_time))**2;
	elif t > conf.pulse_duration - conf.ramp_off_time:
		curField *= sin((conf.pulse_duration - t) * pi / (2*conf.ramp_off_time))**2;
	else:
		curField *= 1

	return curField
