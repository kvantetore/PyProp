#import system modules
import sys
import os
from numpy import conj
import pylab
from numpy import fft

#pytables
import tables

#Pyprop itself
sys.path.append("pyprop")
import pyprop; 
pyprop = reload(pyprop)

def GetMatrixQdot4(psi, conf):
	"""
	Get couplings between the four states considered for
	the optimal control PRL (Saelen et al 2008).
	"""
	scaling = conf.scaling
	V12 = -1.0067555306 * scaling
	V13 = 0.0
	V14 = -0.069677839373 * scaling
	V23 = -1.0501555811  * scaling 
	V24 = 0.0
	V34 = -3.0758641281 * scaling
	potential = array([ [0, V12, V13, V14], \
	                    [V12, 0, V23, V24], \
	                    [V13, V23, 0, V34], \
	                    [V14, V24, V34, 0]],\
	                    dtype=complex)
	
	return potential

def GetDiagonalElementsQdot4(psi, config, potential):
	"""
	"""
	#potential[:] = [0.241995, 0.324479, 0.380237, 0.382176]
	potential[:] = [-0.71337255649561E-01, 0.11146318930044E-01, 0.66904900928925E-01, 0.68844014109733E-01]

def GetDiagonalElements(psi, config, potential):
	"""
	A funtion to provide diagonal (energy) matrix elements
	"""
	h5file = tables.openFile(config.filename, "r")
	try:
		potential[:] = h5file.getNode(config.dataset)[:]
		potential[:] *= config.scaling	
	finally:
		h5file.close()
	#potential[:] = pylab.load(config.file_name)[:config.size] * config.scaling


def Setup(**args):
	"""
	Setup optimal control problem
	"""
	configFile = 'config.ini'
	if "config" in args:
		configFile = args["config"]

	conf = pyprop.Load(configFile)

	if "dt" in args:
		conf.Propagation.timestep = args["dt"]

	controlAlgorithm = args.get("controlAlgorithm", "Krotov")
	confSection = eval("conf.%s" % controlAlgorithm)

	if args.get("perturbControl", False):
		confSection.perturb_control = 1e-15
	
	if "bwdUpdate" in args:
		confSection.update_backwards = args["bwdUpdate"]

	confSection.max_iterations = args.get("maxIter", confSection.max_iterations)

	prop = pyprop.Problem(conf)
	prop.SetupStep()

	controlSolver = eval("pyprop.%s(prop)" % controlAlgorithm)
	controlSolver.ApplyConfigSection(confSection)
	controlSolver.Setup()

	return controlSolver


def SetupZhuRabitz(**args):
	"""
	Setup ZhuRabitz problem
	"""
	conf = pyprop.Load('config.ini')

	if "timestep" in args:
		config.Propagation.timestep = args["timestep"]

	prop = pyprop.Problem(conf)
	prop.SetupStep()
	zhurabitz = pyprop.ZhuRabitz(prop)
	return zhurabitz

def SetupDegani(**args):
	"""
	Setup Degani problem
	"""
	conf = pyprop.Load('config.ini')

	if "timestep" in args:
		config.Propagation.timestep = args["timestep"]


	prop = pyprop.Problem(conf)
	prop.SetupStep()
	degani = pyprop.Degani(prop)
	return degani

def Run():
	krotov = Setup()
	krotov.Run()

	return krotov

def GetPulseSpectrum(krotov, whichControl):
	"""
	"""
	
	dw = 2 * pi * fft.fftshift(fft.fftfreq(len(krotov.ControlVectors[whichControl]), krotov.TimeGridResolution))
	pulseSpectrum = fft.fftshift(fft.fft(krotov.ControlVectors[whichControl]))

	return dw, pulseSpectrum


def Commutator(A,B):
	"""
	Computes the commutator [A,B] = AB-BA
	"""
	return dot(A,B) - dot(B, A)


def NestedCommutatorTypeI(A, B, depth):
	"""
	Computes the nested commutator [A,[A,...,[A,B]...]],
	that is, an inner commutator of C = [A,B] and then
	then recursively the commutator [A,C] until desired 
	depth is reached.
	"""
	C = Commutator(A,B)
	for d in range(depth):
		C = Commutator(A,C)
	
	return C
	

def TextToHDFDense(fileName, vectorSize, scaling):

	groupName = 'doubledot'
	datasetPath = '/' + groupName + '/matrixelements_50'
	data = pylab.load(fileName)

	fileh5 = tables.openFile(fileName + '.h5', 'w')
	
	try:
		group = fileh5.createGroup(fileh5.root, groupName)
		h5array = pyprop.serialization.CreateDataset(fileh5, datasetPath, (vectorSize,vectorSize))

		#Fill h5array with matrix element data
		for i in range(shape(data)[0]):
			row = int(data[i,0]) - 1
			col = int(data[i,1]) - 1
			matel = data[i,2] * scaling
			h5array[row, col] = matel 
			h5array[col,row] = matel

	finally:
		fileh5.close()

def RunStabilityExperiment(**args):
	args["perturbControl"] = True
		
	controlProblem = Setup(**args)

	controlProblem.Run()
	v1 = controlProblem.ControlVectors[0,:].copy()

	controlProblem = Setup(**args)
	controlProblem.Run()
	v2 = controlProblem.ControlVectors[0,:].copy()

	semilogy(linspace(0, controlProblem.PropagationTime, controlProblem.TimeGridSize), abs(v1 - v2))

	return v1, v2

def RunExperiment1():
	"""

	"""
	
	workDir = "experiments/number_1"

	#Load config and setup problem
	conf = pyprop.Load("%s/config.ini" % workDir)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	
	h5file = tables.openFile("%s/results.h5" % workDir, "w")
	try:
		#Run Krotov
		krotov = pyprop.Krotov(prop)
		krotov.Run()
		h5file.createGroup("/","krotov")
		h5file.createArray("/krotov", "J", krotov.J)
		h5file.createArray("/krotov", "Yield", krotov.Yield)
		h5file.createArray("/krotov", "Control", krotov.ControlVectors)
	
		#Run Krotov w/backward update
		prop = pyprop.Problem(conf)
		prop.SetupStep()
		prop.Config.Krotov.update_backwards = True
		krotovb = pyprop.Krotov(prop)
		krotovb.Run()
		h5file.createArray("/krotov_b", "J", krotovb.J)
		h5file.createArray("/krotov_b", "Yield", krotovb.Yield)
		h5file.createArray("/krotov_b", "Control", krotovb.ControlVectors)

		#Run Zhu-Rabitz
		prop = pyprop.Problem(conf)
		prop.SetupStep()
		zhurabitz = pyprop.ZhuRabitz(prop)
		zhurabitz.Run()
		h5file.createGroup("/","zhurabitz")
		h5file.createArray("/zhurabitz", "J", zhurabitz.J)
		h5file.createArray("/zhurabitz", "Yield", zhurabitz.Yield)
		h5file.createArray("/zhurabitz", "Control", zhurabitz.ControlVectors)

		#Run Degani
		prop = pyprop.Problem(conf)
		prop.SetupStep()
		degani = pyprop.Degani(prop)
		degani.Run()
		h5file.createGroup("/","degani")
		h5file.createArray("/degani", "J", degani.J)
		h5file.createArray("/degani", "Yield", degani.Yield)
		h5file.createArray("/degani", "Control", degani.ControlVectors)

		#Store config file
		h5file.createGround("/", "config")
		h5file.setNodeAttr("/config", "configObject", conf.cfgObj)

	finally:
		h5file.close()
