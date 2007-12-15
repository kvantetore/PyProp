"""
The simplest possible pyprop project. It has no compiled dependencies, only this file
and the configuration files.

Calling FindGroundstate() will calculate the ground state of 2d hydrogen in cartesian
coordinates, and save the final state to disk.

Calling Propagate() will load this state from disk to a larger grid, using a custom
initial condition, and propagate it further to demonstrate that it is in fact an 
eigenstate of the system

"""

#import system modules
import sys
import os

#numerical modules
from numpy import *
import pylab

#hdf5
import tables

#home grown modules
pyprop_path = "../../../"
sys.path.insert(1, os.path.abspath(pyprop_path))
import pyprop
pyprop = reload(pyprop)

def CustomInitialCondition(psi, conf):
	"""
	Loads a smaller wavefunction of same grid size into the current
	wavefunction
	"""
	filename = conf.filename
	dataset = conf.dataset
	rank = psi.GetRank()
	psi.GetData()[:] = 0

	localShape = psi.GetData().shape

	localSlice = []
	originalSlice = []

	repr = psi.GetRepresentation()
	distr = repr.GetDistributedModel()
	eps = 1e-10

	for i in range(psi.GetRank()):
		#size and shape of the original wavefunction
		origRank = conf.Get("orig_rank%i" % i)
		origMin = origRank[0]
		origMax = origRank[1]
		origCount = origRank[2]
		origDx =  (origMax - origMin) / float(origCount)

		#size and shape of the current wavefunction
		curRange = repr.GetRange(i)
		curDx = curRange.Dx
	
		if abs(origDx - curDx) > eps:
			raise Exception("Grid spacing in original file (%s) is different from current dx (%s) in rank %i" % (origDx, curDx, i))
		dx = curDx
	
		#origIndex is the start/stop index of the loaded wavefunction in the
		#global index space
		origStartIndex = int( (origMin - curRange.Min) / dx )
		origEndIndex = origStartIndex + origCount
	
		#localIndex is the start/stop index of the local processor in the
		#global index space
		localStartIndex = distr.GetLocalStartIndex(curRange.Count, i)
		localEndIndex = localStartIndex + localShape[i]
	
		#startIndex, endIndex are indices in the global index space
		startIndex = maximum(localStartIndex, origStartIndex)
		endIndex = minimum(localEndIndex, origEndIndex)
	
		#Check if any piece of the wavefunction is in
		if startIndex < localEndIndex and endIndex >= localStartIndex : 
			#indices in the procesor local index space
			siLocalNew  = startIndex - localStartIndex
			eiLocalNew  = endIndex - localStartIndex
			localSlice.append(slice(siLocalNew, eiLocalNew))
	
			#indices in the original index space
			siOrig = startIndex - origStartIndex
			eiOrig = endIndex - origStartIndex
			originalSlice.append(slice(siOrig, eiOrig))
	
		else:	
			#We don't have any overlap between the original grid and the
			#processor local grid, and can therefore terminate
			return
	
	#If we're here, we have at least some overlap between
	#the processor local grid and the original grid
	f = tables.openFile(filename, "r")
	try:
		origData = f.getNode(dataset)
		psi.GetData()[localSlice] = origData[tuple(originalSlice)]
	finally:
		f.close()


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
	if pyprop.ProcId == 0:
		if os.path.exists("groundstate.h5"):
			os.unlink("groundstate.h5")
	prop.SaveWavefunctionHDF("groundstate.h5", "wavefunction")

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
	for t in prop.Advance(100):
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		if pyprop.ProcId == 0:
			print "t = ", t, ", P(t) = ", corr
		#pyprop.Plot2DFull(prop)

	
