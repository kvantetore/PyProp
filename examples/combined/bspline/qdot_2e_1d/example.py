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
	if "dt" in args:
		timeStep = args['dt']
		conf.Propagation.timestep = timeStep

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

	if "omega_left" in args:
		omegaLeft = args["omega_left"]
		conf.QuantumDotPotential_0.omega_left = omegaLeft
		conf.QuantumDotPotential_1.omega_left = omegaLeft

	if "xmax" in args:
		xMax = args['xmax']
		conf.BSplineRepresentation.xmax = xMax
	
	if "xsize" in args:
		xSize = args['xsize']
		conf.BSplineRepresentation.xsize = xSize
	
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


def FindEigenvalues(**args):
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	solver.Solve()
	
	for E in solver.GetEigenvalues():
		print "%.17f" % E

	return solver


def SaveWavefunction(filename, dataset, psi):
	if pyprop.ProcId == 0:
		if os.path.exists(filename):
			os.unlink(filename)
	pyprop.serialization.SaveWavefunctionHDF(filename, dataset, psi)


def GetHamiltonMatrixSubspace(prop, I, J):
	size = prop.psi.GetData().shape[0]
	matrix = zeros((size, size), dtype=complex)
	tempPsi = prop.GetTempPsi()

	for i in range(size):
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[i,I] = 1

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		
		matrix[:, i] = tempPsi.GetData()[:,J]
		
	return matrix


def GetHamiltonMatrix(prop):
	BSplineObject = prop.psi.GetRepresentation().GetRepresentation(0).GetBSplineObject()
	k = BSplineObject.MaxSplineOrder

	size = prop.psi.GetData().size
	wfSize = prop.psi.GetData().shape[0]
	tempPsi = prop.GetTempPsi()
	matrix = zeros((size, size))
	matrix[:,:] = 0

	for i in range(wfSize):
		print "%i of %i" % (i, wfSize)
		for j in range(wfSize):
			prop.psi.GetData()[:] = 0
			prop.psi.GetData()[i,j] = 1

			tempPsi.GetData()[:] = 0
			prop.MultiplyHamiltonian(tempPsi)
			
			startI = j * wfSize
			stopI = startI + wfSize
			for k in range(wfSize):
				J = k * wfSize + i
				matrix[startI:stopI, J] = tempPsi.GetData()[:,k]

	return matrix


def FindEigenvalues(**args):
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	solver.Solve()
	
	for E in solver.GetEigenvalues():
		print "%.17f" % E

	return solver


def TestStability(**args):
	args["imtime"] = True
	args["omega_left"] = 1
	initProp = FindGroundstate(**args)

	args["imtime"] = False
	args["omega_left"] = 1.5
	conf = SetupConfig(**args)
	#conf.Propagation.potential_evaluation = ["StarkPotential"]
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	
	prop.psi.GetData()[:] = initProp.psi.GetData()
	initPsi = prop.psi.Copy()

	for t in prop.Advance(10):
		print "t = %.2f, N(t) = %.8f, P(t) = %.8f" % (t, prop.psi.GetNorm(), abs(prop.psi.InnerProduct(initPsi)**2))
	return prop

	
import time

def TestInnerProductSpeed(**args):
	prop = SetupProblem(**args)

	avgCount = 5
	minCount = 10

	for algo in range(4):
		prop.psi.GetRepresentation().Algorithm = algo
		minT = 1e10
		for i in range(minCount):
			t = - time.time()
			for j in range(avgCount):
				n = prop.psi.GetNorm()
			t += time.time()
			if t<minT:
				minT = t / avgCount
		print "Algorithm %i: %f" % (algo, minT)
	

