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
	E1 = E[I[state1]] 
	E2 = E[I[state2]] 
	deltaE = E1 - E2

	#Strip phase (rotate to real)
	for i in [state1, state2]:
		largestCompIndex = argmax(abs(V[:,I[i]]))
		x = V[largestCompIndex, I[i]].real
		y = V[largestCompIndex, I[i]].imag
		eigenvectorPhase = arctan(y/x) 
		V[:,I[i]] *= exp(-1j * eigenvectorPhase)

	prop.V1 = V[:,I[state1]]
	prop.V2 = V[:,I[state2]]


	#Create wavefunction from state 1 and 2
	psi1 = prop.psi.Copy()
	psi1.GetData()[:] = V[:,I[state1]].real
	psi1.Normalize()	
	psi2 = prop.psi.Copy()
	psi2.GetData()[:] = V[:,I[state2]].real
	psi2.Normalize()

	#Create initial (linear combination) state
	prop.psi.GetData()[:] = 0.0
	prop.psi.GetData()[:] = psi1.GetData()[:] + psi2.GetData()[:]
	prop.psi.Normalize()
	
	#Store initial state
	initialPsi = prop.psi.Copy()

	#Integrate
	prop.initialCorr = []
	prop.state1Corr = []
	prop.state2Corr = []
	prop.outputTimes = []
	prop.analyticCorr = []
	prop.norm = []
	for t in prop.Advance(300): 
		prop.norm.append( prop.psi.GetNorm() )
		prop.initialCorr.append( prop.psi.InnerProduct(initialPsi) )
		prop.state1Corr.append( prop.psi.InnerProduct(psi1) )
		prop.state2Corr.append( prop.psi.InnerProduct(psi2) )
		prop.outputTimes.append( t )
		prop.analyticCorr.append( 0.5 * (exp(-1j * E1 * t) + exp(-1j * E2 * t)) )

	#
	#Plot results
	#
	figure()

	# Plot projection on states 1 and 2
	subplot(211)
	#plot(prop.outputTimes, abs(array(prop.state1Corr))**2)
	#plot(prop.outputTimes, abs(array(prop.state2Corr))**2)
	plot(prop.outputTimes, abs(array(prop.initialCorr))**2)
	ylabel('Autocorrelation')
	axis('tight')

	#Plot relative phase error
	phase = arctan2(array(prop.initialCorr).imag, array(prop.initialCorr).real)
	phaseAnalytic = arctan2(array(prop.analyticCorr).imag, array(prop.analyticCorr).real)
	phaseError = abs(phase - phaseAnalytic)
	maxPhaseError = max(abs(phaseError))
	subplot(212)
	plot(prop.outputTimes, phaseError / maxPhaseError)
	xlabel('Time (a.u.)')
	ylabel('Rel. phase err. (%1.1e)' % maxPhaseError)
	axis('tight')

	return prop


def TestPhaseAccuracySingleState(**args):
	"""
	Integrate an eigenstates in	time, projecting on the initial state 
	during propagation. This allows us to check the phase accuracy of 
	the integrator, which should give some indication of its performance 
	(and correctness).
	"""

	#Select state to integrate
	state = 0
	if 'state' in args:
		state = args['state']
	
	#Set a sufficient integration time
	args['duration'] = 100
	args['dt'] = 0.01

	#Set up problem
	prop = SetupProblem(**args)

	#Get eigenstates
	M = GetHamiltonMatrix(prop)
	E, V = eig(M)
	I = argsort(E)
	energy = E[I[0]]

	#Strip phase (rotate to real)
	largestCompIndex = argmax(abs(V[:,I[state]]))
	x = V[largestCompIndex, I[state]].real
	y = V[largestCompIndex, I[state]].imag
	eigenvectorPhase = arctan(y/x) 
	V[:,I[state]] *= exp(-1j * eigenvectorPhase)
	prop.V = V[:,I[state]]

	#Create intial (linear combination) state
	prop.psi.GetData()[:] = 0.0
	prop.psi.GetData()[:] = V[:,I[state]].real
	prop.psi.Normalize()
	
	#Store initial state
	initialPsi = prop.psi.Copy()

	#Integrate
	prop.initialCorr = []
	prop.norm = []
	prop.outputTimes = []
	prop.analyticCorr = []

	for t in prop.Advance(400): 
		prop.norm.append( prop.psi.GetNorm() )
		prop.initialCorr.append( prop.psi.InnerProduct(initialPsi) )
		prop.outputTimes.append( t )
		prop.analyticCorr.append( exp(-1j * energy * t) )

	#
	#Plot results
	#
	figure()

	#Plot norm and projection on initial state
	subplot(211)
	plot(prop.outputTimes, abs(array(prop.norm)))
	plot(prop.outputTimes, abs(array(prop.initialCorr))**2)
	legend(('Norm', 'Autocorrelation'))
	axis('tight')

	#Plot relative phase error
	phase = arctan2(array(prop.initialCorr).imag, array(prop.initialCorr).real)
	phaseAnalytic = arctan2(array(prop.analyticCorr).imag, array(prop.analyticCorr).real)
	phaseError = abs(phase - phaseAnalytic)
	maxPhaseError = max(abs(phaseError))
	subplot(212)
	plot(prop.outputTimes, phaseError / maxPhaseError)
	xlabel('Time (a.u.)')
	ylabel('Rel. phase err. (%1.1e)' % maxPhaseError)
	axis('tight')

	return prop

