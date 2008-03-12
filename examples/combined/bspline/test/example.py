#Import system modules
import sys
import os
from pylab import *
from numpy import *

#Load pyprop
sys.path.insert(1, os.path.abspath("./pyprop"))
import pyprop
pyprop = reload(pyprop)

#Load common potentials
sys.path.append(os.path.abspath('../potentials'))
from libcommonpotentials import *

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

	SaveWavefunction("groundstate.h5", "/wavefunction", prop.psi)

	return prop


def SaveWavefunction(filename, dataset, psi):
	"""
	Save wavefunction using HDF5 
	"""

	if pyprop.ProcId == 0:
		if os.path.exists(filename):
			os.unlink(filename)
	pyprop.serialization.SaveWavefunctionHDF(filename, dataset, psi)


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


def TestInnerProduct(**args):
	"""
	Compute inner product between two given states, using the 
	psi.InnerProduct() function. States are given by the 'states' 
	keyword (a two-element list). States are obtained by constructing
	the hamilton matrix and diagonalizing.
	"""

	state1 = 0
	state2 = 1
	if 'states' in args:
		state1 = args['states'][0]
		state2 = args['states'][1]
	
	prop = SetupProblem(**args)

	M = GetHamiltonMatrix(prop)
	E, V = eig(M)
	I = argsort(E)

	tempPsi = prop.psi.Copy()
	tempPsi.GetData()[:] = V[:,I[state1]]
	tempPsi.Normalize()
	prop.psi.GetData()[:] = V[:,I[state2]]
	prop.psi.Normalize()

	figure()
	plot(abs(tempPsi.GetData()))
	plot(abs(prop.psi.GetData()))

	print
	print "Energies = ", abs(E[I[state1]]), abs(E[I[state2]])
	print "Innerproduct <%i|%i> = %1.15e" % (state1, state2, abs(prop.psi.InnerProduct(tempPsi)))

	return prop


def TestPhaseAccuracy(**args):
	"""
	Integrate an initial linear combination of two eigenstates in
	time, projecting on the initial state during propagation. This
	allows us to check the phase accuracy of the integrator, which
	should give some indication of its performance (and correctness).

	Projection on initial state of a linear combination of two eigen-
	states are at time t given by:

	|< initial | psi(t) >|**2 = 0.5 * (1 + cos(dE * t)),  dE = E1 - E2
	"""

	#Select states to integrate
	state1 = 0
	state2 = 1
	if 'states' in args:
		state1 = args['states'][0]
		state2 = args['states'][1]
	
	
	#Set a sufficient integration time
	args['duration'] = 100
	args['dt'] = 0.01

	#Set up problem
	prop = SetupProblem(**args)

	#Get eigenstates
	M = GetHamiltonMatrix(prop)
	E, V = eig(M)
	I = argsort(E)
	deltaE = E[I[state1]] - E[I[state2]]

	#Create intial (linear combination) state
	prop.psi.GetData()[:] = V[:,I[state1]] + V[:,I[state2]]
	prop.psi.Normalize()
	
	#Store initial state
	initialPsi = prop.psi.Copy()

	#Integrate
	prop.initialCorr = []
	prop.outputTimes = []
	prop.analyticCorr = []
	for t in prop.Advance(300): 
		prop.initialCorr.append( prop.psi.InnerProduct(initialPsi) )
		prop.outputTimes.append( t )
		prop.analyticCorr.append( 0.5 * (1 + cos(deltaE * t)) )

	#Plot result
	figure()
	subplot(211)
	plot(prop.outputTimes, abs(array(prop.initialCorr))**2)
	plot(prop.outputTimes, array(prop.analyticCorr))
	subplot(212)
	semilogy(prop.outputTimes, abs(abs(array(prop.initialCorr))**2 - array(prop.analyticCorr)))

	return prop

