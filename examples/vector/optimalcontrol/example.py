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


def GetDiagonalElements(psi, config, potential):
	"""
	A funtion to provide diagonal (energy) matrix elements
	"""
	h5file = tables.openFile(config.filename, "r")
	try:
		potential[:] = h5file.getNode(config.dataset)[:config.size]
		potential[:] *= config.scaling	
	finally:
		h5file.close()


def Setup(**args):
	"""
	Setup optimal control problem
	"""
	print "Setting up control problem"

	configFile = 'config.ini'
	if "config" in args:
		configFile = args["config"]

	conf = pyprop.Load(configFile)

	if "dt" in args:
		conf.Propagation.timestep = args["dt"]

	controlAlgorithm = args.get("controlAlgorithm", "Krotov")
	print "    Using control algorithm: %s" % controlAlgorithm
	confSection = eval("conf.%s" % controlAlgorithm)

	if "bwdUpdate" in args:
		confSection.update_backwards = args["bwdUpdate"]
		print "Backward update: %s" % args["bwdUpdate"]

	confSection.max_iterations = args.get("maxIter", confSection.max_iterations)

	prop = pyprop.Problem(conf)
	prop.SetupStep()

	controlSolver = eval("pyprop.%s(prop)" % controlAlgorithm)
	controlSolver.ApplyConfigSection(confSection)
	controlSolver.Setup()

	return controlSolver


def Run():
	krotov = Setup()
	krotov.Run()

	return krotov

def GetPulseSpectrum(krotov, whichControl):
	dw = 2 * pi * fft.fftshift(fft.fftfreq(len(krotov.ControlVectors[whichControl]), \
		krotov.TimeGridResolution))
	pulseSpectrum = fft.fftshift(fft.fft(krotov.ControlVectors[whichControl]))
	return dw, pulseSpectrum


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
