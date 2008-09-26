#improt system modules
import sys
import os

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("./pyprop"))

#Load and reload pyprop in order to get recent changes
import pyprop
pyprop = reload(pyprop)
from libh2p import *

#numpy an pylab for good measure
from pylab import *
from numpy import *

def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	return prop

def FindEigenvalues(**args):
	prop = SetupProblem(**args)
	tempPsi = prop.GetTempPsi()

	solver = pyprop.ArpackSolver(prop)
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
	silent = pyprop.ProcId != 0
	args["silent"] = silent

	#for finding ground state, we always use imaginary time
	args['imTime'] = True  #Remember: Arguments are case sensitive

	#Create the propagator
	prop = SetupProblem(**args)
	prop.AdvanceStep()

	#Estimate runtime by running a few timesteps
	stepCount = 50
	timestepDuration = - time.time()
	for i in range(stepCount):
		prop.AdvanceStep()
	timestepDuration += time.time()
	timestepDuration /= stepCount

	prop.psi.GetRepresentation().GetDistributedModel().GlobalBarrier()
	
	estimatedRuntime = abs(prop.Duration / prop.TimeStep) * timestepDuration
	fmt = time.gmtime(estimatedRuntime)
	PrintOut("Estimated Runtime = %i days, %i hours, %i minutes (%f secs), timestep = %f" % (fmt[2]-1, fmt[3], fmt[4], estimatedRuntime, timestepDuration))

	return 

	#propagate
	for t in prop.Advance(10):
		PrintOut( "t = %4.3f, E = %.18f" % (t, prop.GetEnergyExpectationValue()) )

	return prop

def PrintOut(output):
	if pyprop.ProcId == 0:
		print output

import time

def Propagate(**args):
	#Setup Problem
	args['imTime'] = False
	args["silent"] = pyprop.ProcId != 0
	args["additionalPotentials"] = ["LaserPotential", "AbsorbingBoundary"]
	prop = SetupProblem(**args)

	electronicRadialRank = 0
	electronicAngularRank = 1

	radialBoxSize = 80
	angularBoxSize = 16

	#Setup box potentials
	radialBox = SetupStepPotential(prop.psi, -radialBoxSize, radialBoxSize, electronicRadialRank)
	angularBox = SetupStepPotential(prop.psi, -1, angularBoxSize, electronicAngularRank)

	#Initial Wavefunction
	initPsi = args.get('initPsi', None)
	if initPsi == None:
		initPsi = prop.psi.Copy()
	initPsi.Normalize()
	prop.psi.GetData()[:] = initPsi.GetData()

	#Estimate runtime by running a few timesteps
	stepCount = 10
	timestepDuration = - time.time()
	for i in range(stepCount):
		prop.AdvanceStep()
	timestepDuration += time.time()
	timestepDuration /= stepCount

	estimatedRuntime = abs(prop.Duration / prop.TimeStep) * timestepDuration
	fmt = time.gmtime(estimatedRuntime)
	PrintOut("Estimated Runtime = %i days, %i hours, %i minutes" % (fmt[2]-1, fmt[3], fmt[4]))

	def output():
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		angularBoxPop = angularBox.GetExpectationValue(prop.psi, 0, 0)
		radialBoxPop = radialBox.GetExpectationValue(prop.psi, 0, 0)
		PrintOut( "t = %4.3f, N(t) = %.5f, C(t) = %.10f, radBox(t) = %.10f, angBox(t) = %.10f" % (prop.PropagatedTime, norm, corr, angularBoxPop, radialBoxPop) )

	#Output once before the propagation
	output()

	#Propagated to end with 100 printouts
	for t in prop.Advance(100):
		output()

	#Output once after the propagation
	output()
	
	return prop


def SetupStepPotential(psi, zeroBefore, zeroAfter, stepRank):
	"""
	Setup a step potential which is zero outside (zeroBefore -> zeroAfter)
	in rank stepRank, and one otherwise. 

	i.e. to find the amount of psi outside r=80, (and r is rank 0)

	tempPsi = prop.GetTempPsi()
	pot = SetupStepPotential(prop.psi, -80, 80, 0)

	tempPsi.Clear()
	normInside = pot.GetExpectationValue(prop.psi, 0, 0)

	"""
	class stepSection(pyprop.Section):
		def __init__(self):
			self.type = pyprop.PotentialType.Dynamic
			self.classname = "StepPotential"
			self.zero_before = zeroBefore
			self.zero_after = zeroAfter
			self.step_rank = stepRank

	section = stepSection()		
	potential = pyprop.CreatePotentialFromSection(section, "H2MaskPotential", psi)
	potential.SetupStep(0)

	return potential




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
		configFile = args['configFile']
	else:
		configFile = "config.ini"
	conf = pyprop.Load(configFile)
	
	#Grid size
	if 'electronicGridSize' in args:
		gridSize = args['electronicGridSize']
		if not silent: print "Using Electronic Grid Size = ", gridSize
		rank0 = conf.ElectronicRadialRepresentation.rank0
		rank0[2] = gridSize
		conf.SetValue("ElectronicRadialRepresentation", "rank0", rank0)

	if 'nuclearGridSize' in args:
		gridSize = args['nuclearGridSize']
		if not silent: print "Using Nuclear Grid Size = ", gridSize
		rank0 = conf.NuclearRadialRepresentation.rank0
		rank0[2] = gridSize
		conf.SetValue("NuclearRadialRepresentation", "rank0", rank0)

	if 'lmax' in args:
		lmax = args['lmax']
		if not silent: print "Using LMax = ", lmax
		conf.SetValue("ElectronicAngularRepresentation", "maxl", lmax)

	#TimeStep
	if 'dt' in args:
		dt = args['dt']
		if not silent: print "Using TimeStep = ", dt
		conf.SetValue("Propagation", "timestep", dt)

	#Imaginary Time Propagation
	if 'imTime' in args:
		imTime = args['imTime']
		if not silent: print "Using ImaginaryTime = ", imTime
		dt = conf.Propagation.timestep
		if imTime:
			conf.SetValue("Propagation", "renormalization", True)
			conf.SetValue("Propagation", "timestep",  -1.0j * abs(dt))
		else:
			conf.SetValue("Propagation", "renormalization", False)
			conf.SetValue("Propagation", "timestep",  abs(dt))

	if "additionalPotentials" in args:
		additionalPotentials = args["additionalPotentials"]
		pot = conf.Propagation.potential_evaluation + list(additionalPotentials)
		conf.SetValue("Propagation", "potential_evaluation", pot)

	if 'silent' in args:
		silent = args['silent']
		conf.SetValue("Propagation", "silent", silent)

		#H2+ Potential
	if 'softing' in args:
		softing = args['softing']
		if not silent: print "Using softing ", softing
		conf.SetValue("Potential", "softing", softing)

	#Laser Potential
	if 'duration' in args:
		duration = args['duration']
		conf.SetValue("Propagation", "duration", duration)

	
	print "Setup Config Complete"
	print ""
	return conf


