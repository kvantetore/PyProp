#Import system modules
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
		
	return conf


def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop


def FindGroundstate(**args):
	args['imtime'] = True
	prop = SetupProblem(**args)
	
	for t in prop.Advance(10):
		E = prop.GetEnergy()
		print "t = %f, E = %f" % (t, E)

	E = prop.GetEnergy()
	print "Ground State Energy = %f" % E

	SaveWavefunction("groundstate.h5", "/wavefunction", prop.psi)

	return prop


def SaveWavefunction(filename, dataset, psi):
	if pyprop.ProcId == 0:
		if os.path.exists(filename):
			os.unlink(filename)
	pyprop.serialization.SaveWavefunctionHDF(filename, dataset, psi)


def FindIonizationProbability(**args):
	args['imtime'] = False
	args['config'] = "propagation.ini"
	prop = SetupProblem(**args)

	initialPsi = prop.psi.Copy()

	for t in prop.Advance(10):
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initialPsi))**2
		print "t = %f, Norm = %f, Corr = %f" % (t, norm, corr)

	norm = prop.psi.GetNorm()
	corr = abs(prop.psi.InnerProduct(initialPsi))**2
	print "Ionization Probability = %f" % norm
	print "Initial state correlation = %f" % corr

	return prop


def GetHamiltonMatrix(prop):
	size = prop.psi.GetData().shape[0]
	matrix = zeros((size, size), dtype=complex)
	tempPsi = prop.GetTempPsi()

	for i in range(size):
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[i,1] = 1

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		
		matrix[:, i] = tempPsi.GetData()[:,1]
		
	return matrix

