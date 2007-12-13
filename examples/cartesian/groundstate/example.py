"""
The simplest possible pyprop project. It has no compiled dependencies, only this file
and the configuration files.

Calling FindGroundstate() will calculate the ground state of 2d hydrogen in cartesian
coordinates, and save the final state to disk.

Calling Propagate() will load this state from disk, and propagate it further to 
demonstrate that it is in fact an eigenstate of the system

"""

#import system modules
import sys
import os

#numerical modules
from numpy import *
import pylab

#home grown modules
pyprop_path = "../../../"
sys.path.insert(1, os.path.abspath(pyprop_path))
import pyprop
pyprop = reload(pyprop)

def FindEnergy(prop):
	prop.psi.Normalize()
	
	psi = prop.psi.Copy()

	#potential energy:
	curE = 0
	for pot in prop.Propagator.PotentialList:
		potEval = pot.PotentialEvaluator
		curE += potEval.CalculateExpectationValue(prop.psi, prop.TimeStep, prop.PropagatedTime)
		print curE

	
	prop.Propagator.FFTTransform.ForwardTransform(prop.psi)
	prop.Propagator.ChangeRepresentation()

	pot = prop.Propagator.KineticPotential.GetPotential(prop.TimeStep)
	prop.psi.GetData()[:] *= pot

	prop.Propagator.FFTTransform.InverseTransform(prop.psi)
	prop.Propagator.ChangeRepresentation()

	curE += real(prop.psi.InnerProduct(psi))
	print curE

	prop.psi.GetData()[:] = psi.GetData()


def FindGroundstate(**args):
	"""
	Loads the configuration file "find_groundstate.ini", which contains information
	on how to find the ground state of 2d hydrogen, by imaginary time propagation.

	When the problem is fully advanced, the wavefunction is saved to the file 
	"groundstate.dat", and the ground state energy is written to screen

	finally it returns the Problem object prop back to the caller for further processing
	"""
	
	#load config
	conf = pyprop.Load("find_groundstate.ini")
	silent = False
	if 'silent' in args:
		silent = args['silent']
		conf.Propagation.silent = silent
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#propagate to find ground state
	for t in prop.Advance(10):
		E = prop.GetEnergy()
		if not silent:
			print "t =", t, "E =", E
	
	#save groundstate to disk
	prop.SaveWavefunctionHDF("groundstate.h5", "groundstate")

	#Find energy
	energy = prop.GetEnergy()
	if not silent:
		print "Groundstate energy:", energy, "a.u."

	return prop, energy

def Test():
	conf = pyprop.Load("find_groundstate.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	
	print prop.GetEnergyExpectationValue()
	print prop.GetEnergyExpectationValue()



def FindEigenvalues():
	conf = pyprop.Load("find_groundstate.ini")
	conf.Propagation.silent = pyprop.ProcId != 0
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	solver = pyprop.PiramSolver(prop)
	solver.Solve()

	print "Eigenvalues = ", solver.Solver.GetEigenvalues().real
	return solver

def GetEigenvalue(solver, num):

	prop = solver.BaseProblem
	psi = prop.psi
	psi.GetData()[:] = solver.Solver.GetEigenvectors()[num,:].reshape(psi.GetData().shape)
	psi.Normalize()

	pyprop.Plot2DFull(prop)
	print "Energy = %s" % prop.GetEnergy()



def Propagate():
	conf = pyprop.Load("propagation.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#Create a copy of the wavefunction so that we can calculate
	#the autocorrelation function during propagation
	initPsi = prop.psi.Copy()

	#Propagate the system to the time specified in propagation.ini,
	#printing the autocorrelation function, and plotting the wavefunction
	#10 evenly spaced times during the propagation
	for t in prop.Advance(10):
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		print "t = ", t, ", P(t) = ", corr
		#pyprop.Plot2DFull(prop)

	
