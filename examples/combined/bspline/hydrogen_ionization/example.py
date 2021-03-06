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

	if "duration" in args:
		duration = args["duration"]
		conf.Propagation.duration = duration
		
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

	prop.SaveWavefunctionHDF("groundstate.h5", "/wavefunction")

	return prop


def FindEigenvalues(**args):
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	solver.Solve()
	print solver.GetEigenvalues()

	solver.SetEigenvector(prop.psi, 0)
	prop.SaveWavefunctionHDF("groundstate.h5", "/wavefunction")

	return solver



#def SaveWavefunction(filename, dataset, psi):
#	if pyprop.ProcId == 0:
#		if os.path.exists(filename):
#			os.unlink(filename)
#	pyprop.serialization.SaveWavefunctionHDF(filename, dataset, psi)

def FindIonizationProbability(**args):
	args['imtime'] = False
	args['config'] = "propagation.ini"
	prop = SetupProblem(**args)

	initialPsi = prop.psi.Copy()

	timeList = []
	corrList = []
	hold(False)
	for t in prop.Advance(30):
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initialPsi))**2
		timeList.append(t)
		corrList.append(corr)
		print "t = %f, Norm = %f, Corr = %f" % (t, norm, corr)
		pcolormesh(abs(prop.psi.GetData()))

	norm = prop.psi.GetNorm()
	corr = abs(prop.psi.InnerProduct(initialPsi))**2
	print "Final norm = %f" % norm
	print "Initial state correlation = %f" % corr

	prop.TimeList = array(timeList)
	prop.CorrelationList = array(corrList)

	return prop

def FindIonizationProbabilityAmplitude():
	amplitudeList = r_[0:2:0.2]
	ionizationList = zeros(len(amplitudeList), dtype=double)

	for i in range(len(amplitudeList)):
		prop = FindIonizationProbability(amplitude=amplitudeList[i])
		ionizationList[i] = 1.0 - prop.psi.GetNorm()

	plot(amplitudeList, ionizationList, label="Ionization Probability")
	xlabel("Electric Field Strength")
	ylabel("Ionization Probability")
	legend()

	return amplitudeList, ionizationList


def GetHamiltonMatrixSubspace(prop, l):
	size = prop.psi.GetData().shape[0]
	matrix = zeros((size, size), dtype=complex)
	tempPsi = prop.GetTempPsi()

	for i in range(size):
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[i,l] = 1

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		
		matrix[:, i] = tempPsi.GetData()[:,l]
		
	return matrix
