#!/bin/env python

import sys
sys.path.append("./pyprop")
sys.path.append("../helium_stabilization")
import pyprop
import numpy
from scipy import *
from pylab import *

import pypar
import time
import tables

#from libepetratest import *
from libpotential import *
from pyprop.core import EpetraPotential_3

execfile("../helium_stabilization/preconditioner.py")
execfile("preconditioner.py")

pyprop.ProjectNamespace = globals()

def Print(str="", whichProc = [0]):
	def printMsg():
		if pyprop.ProcId in whichProc:
			print str

	if pyprop.Redirect.redirect_stdout:
		pyprop.Redirect.Disable()
		printMsg()
		pyprop.Redirect.Enable(silent=True)
	else:
		printMsg()
		
pyprop.PrintOut = Print

def timeIt(func):
	t1 = time.time()
	func()
	t2 = time.time()
	Print("  Function '%s' took %4.1f s." % (func.func_name, (t2-t1)))


def FormatDuration(duration):
	duration = int(duration)
	h = (duration / 3600)
	m = (duration / 60) % 60
	s = (duration % 60)

	str = []
	if h>0: str += ["%ih" % h]
	if m>0: str += ["%im" % m]
	if s>0: str += ["%is" % s]

	return " ".join(str)


def SetupConfig(**args):
	configFile = args.get("config", "config.ini")
	#if configfile is a string, load it, otherwise, treat it as
	#a config parser object
	if isinstance(configFile, str):
		conf = pyprop.Load(configFile)
	elif isinstance(configFile, pyprop.Config):
		#Let's make a deep copy here to avoid changing input
		conf = objectcopier.deepcopy(configFile)
	else:
		conf = pyprop.Config(configFile)
	
	if "silent" in args:
		silent = args["silent"]
		conf.SetValue("Propagation", "silent", silent)


	#Update config object from possible changed ConfigParser object
	newConf = pyprop.Config(conf.cfgObj)

	return newConf

	
def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	return prop


def TestMultiplyOverlap():
	pyprop.PrintOut("Now testing S * psi...")
	pyprop.Redirect.Enable(silent=True)

	fileName  = "test_multiplyoverlap.h5"
	seed(0)

	conf = pyprop.Load("config-test.ini")
	psi = pyprop.CreateWavefunction(conf)
	initPsi = pyprop.CreateWavefunction(conf)


	if pyprop.ProcCount == 1:
		psi.GetData()[:] = random(psi.GetData().shape)
		print "Normalizing wavefunction..."
		psi.Normalize()
		initPsi.GetData()[:] = psi.GetData()[:]
		pyprop.serialization.SaveWavefunctionHDF(fileName, "/wavefunction", psi, conf)
		psi.GetRepresentation().MultiplyOverlap(psi)
		pyprop.serialization.SaveWavefunctionHDF(fileName, "/wavefunctionoverlap", psi, conf)
	else:
		pyprop.serialization.LoadWavefunctionHDF(fileName, "/wavefunction", initPsi)
		pyprop.serialization.LoadWavefunctionHDF(fileName, "/wavefunctionoverlap", psi)

	destPsi = initPsi.Copy()
	destPsi.Clear()
	tmpPsi = initPsi.Copy()
	tmpPsi.Clear()

	destPsi.GetData()[:] = initPsi.GetData()[:]
	destPsi.GetRepresentation().MultiplyOverlap(destPsi)

	Print()
	Print("  Proc %s: ||S * psi - S'*psi||_max = %s" % (pyprop.ProcId, numpy.max(numpy.max(psi.GetData() - destPsi.GetData()))))
	Print("  Proc %s: ||S * psi - S'*psi|| = %s" % (pyprop.ProcId, linalg.norm(psi.GetData() - destPsi.GetData())))

	#finish and cleanup
	pypar.barrier()
	pyprop.Redirect.Disable()
	pyprop.PrintOut("\n...done!")


def TestSolveOverlap():
	pyprop.PrintOut("")
	pyprop.PrintOut("Now testing S^-1 * psi...")
	pyprop.Redirect.Enable(silent=True)

	seed(0)

	fileName  = "test_solveoverlap.h5"

	conf = pyprop.Load("config-test.ini")
	psi = pyprop.CreateWavefunction(conf)
	initPsi = pyprop.CreateWavefunction(conf)

	if pyprop.ProcCount == 1:
		psi.GetData()[:] = random(psi.GetData().shape)
		psi.Normalize()
		initPsi.GetData()[:] = psi.GetData()[:]

		#Store initial (random) psi
		pyprop.serialization.SaveWavefunctionHDF(fileName, "/wavefunction", psi, conf)

		#Store S^-1 * psi
		psi.GetRepresentation().SolveOverlap(psi)
		pyprop.serialization.SaveWavefunctionHDF(fileName, "/wavefunctionoverlap", psi, conf)

		#determine overlap matrix condition number
		overlap = psi.GetRepresentation().GetGlobalOverlapMatrix(0)
		A = overlap.GetOverlapBlasBanded()
		B = pyprop.core.ConvertMatrixBlasBandedToFull(A)
		Print("  Overlap matrix condition number = %e" % cond(B))
		
	else:
		pyprop.serialization.LoadWavefunctionHDF(fileName, "/wavefunction", initPsi)
		pyprop.serialization.LoadWavefunctionHDF(fileName, "/wavefunctionoverlap", psi)

	destPsi = initPsi.Copy()
	destPsi.Clear()
	tmpPsi = initPsi.Copy()
	tmpPsi.Clear()

	#Calculate S^-1 * psi
	destPsi.GetData()[:] = initPsi.GetData()[:]
	destPsi.GetRepresentation().SolveOverlap(destPsi)
	tmpPsi.GetData()[:] = destPsi.GetData()[:]

	#Calculate S * S^-1 * psi
	destPsi.GetRepresentation().MultiplyOverlap(destPsi)
	
	Print()
	a = numpy.max(numpy.max(psi.GetData() - tmpPsi.GetData()))
	Print("  Proc %s: ||S^-1 * psi - S'^-1 * psi||_max = %s" % (pyprop.ProcId, a), range(pyprop.ProcCount))
	Print()
	b = numpy.max(numpy.max(initPsi.GetData() - destPsi.GetData()))
	Print("  Proc %s: ||S * S^-1 * psi - I * psi||_max = %s" % (pyprop.ProcId, b), range(pyprop.ProcCount))
	
	c = linalg.norm(initPsi.GetData() - destPsi.GetData())
	Print("  Proc %s: ||S * S^-1 * psi - I * psi|| = %s" % (pyprop.ProcId, c), range(pyprop.ProcCount))

	#finish and cleanup
	pypar.barrier()
	pyprop.Redirect.Disable()
	pyprop.PrintOut("\n...done!")


def TestInnerProduct():
	pyprop.PrintOut("")
	pyprop.PrintOut("Now testing innerproduct...")
	pyprop.Redirect.Enable(silent=True)

	seed(0)

	fileName  = "test_innerproduct.h5"

	conf = pyprop.Load("config-test.ini")
	psi = pyprop.CreateWavefunction(conf)
	tmpPsi = pyprop.CreateWavefunction(conf)

	if pyprop.ProcCount == 1:
		psi.GetData()[:] = random(psi.GetData().shape)
		psi.Normalize()

		tmpPsi.GetData()[:] = random(psi.GetData().shape)
		tmpPsi.Normalize()

		pyprop.serialization.SaveWavefunctionHDF(fileName, "/wavefunction1", psi, conf)
		pyprop.serialization.SaveWavefunctionHDF(fileName, "/wavefunction2", tmpPsi, conf)
	else:
		pyprop.serialization.LoadWavefunctionHDF(fileName, "/wavefunction1", psi)
		pyprop.serialization.LoadWavefunctionHDF(fileName, "/wavefunction2", tmpPsi)

	
	inner1 = psi.InnerProduct(tmpPsi)
	inner2 = psi.InnerProduct(psi)

	Print()
	Print("<psi|tmpPsi> = %s" % inner1, range(pyprop.ProcCount))
	Print("<psi|psi> = %s" % inner2, range(pyprop.ProcCount))

	#finish and cleanup
	pypar.barrier()
	pyprop.Redirect.Disable()
	pyprop.PrintOut("\n...done!")


def TestFindEigenvalues():
	#Custom matvec product
	def matvec(psi, tmpPsi, t, dt):
		tmpPsi.Clear()
		potMatrix.Multiply(psi, tmpPsi)
		tmpPsi.GetRepresentation().SolveOverlap(tmpPsi)

	pyprop.PrintOut("")
	pyprop.PrintOut("Now testing eigenvalue computation...")

	#Setup problem
	Print("  Setting up problem...", [0])
	prop = SetupProblem(config='config_eigenvalues.ini', silent = True)
	psi = prop.psi
	pyprop.Redirect.Enable(silent=True)

	#Setup Epetra potential and copy Tensor potential data into it
	Print("  Converting tensor potential to Epetra matrix...", [0])
	potMatrix = EpetraPotential_3()
	potMatrix.Setup(psi)
	for pot in prop.Propagator.BasePropagator.PotentialList:
		localBasisPairs = pot.BasisPairs
		potMatrix.AddTensorPotentialData(pot.PotentialData, localBasisPairs, 0)
		del pot.PotentialData
	potMatrix.GlobalAssemble()

	#Setup PIRAM
	Print("  Setting up Piram...", [0])
	prop.Config.Arpack.matrix_vector_func = matvec
	solver = pyprop.PiramSolver(prop)

	#Find eigenvalues
	Print("  Calculating eigenvalues...", [0])
	solver.Solve()
	eigs = solver.GetEigenvalues()
	Print("    Eigenvalues = %s" % str(eigs), [0])

	#Test first eigenvalue
	prop.psi.Clear()
	tmpPsi = prop.psi.Copy()
	solver.SetEigenvector(prop.psi, 0)
	prop.psi.Normalize()
	matvec(prop.psi, tmpPsi, 0, 0)
	eigRes = abs(prop.psi.InnerProduct(tmpPsi) - solver.GetEigenvalues()[0])
	Print("    ||H * v - lambda * v|| = %s" % eigRes,[0])

	#finish and cleanup
	pypar.barrier()
	pyprop.Redirect.Disable()
	pyprop.PrintOut("\n...done!")

def TestEpetraMatrix():
	#Setup problem
	prop = SetupProblem(config='config_eigenvalues.ini')
	psi = prop.psi
	tmpPsi = psi.Copy()
	tmpPsi.Clear()

	#Setup Epetra potential
	pot = prop.Propagator.BasePropagator.PotentialList[0]
	localBasisPairs = pot.BasisPairs
	potMatrix = EpetraPotential_3()
	potMatrix.Setup(psi, pot.PotentialData, localBasisPairs, 0)

	#Tensor potential properties
	potNorm = sqrt(sum(abs(pot.PotentialData)**2))
	numElements = pot.PotentialData.size
	
	print
	print "||pot||_F = %s" % potNorm
	print "Number of elements in potential = %s" % numElements
	print

	#Custom matvec product
	def matvec(psi, tmpPsi, t, dt):
		tmpPsi.Clear()
		potMatrix.Multiply(psi, tmpPsi)
		tmpPsi.GetRepresentation().SolveOverlap(tmpPsi)

	#pyprop.serialization.LoadWavefunctionHDF(fileName, "/wavefunction1", psi)
#	h5file = tables.openFile("eig_test.h5")
#	try:
#		for i in range(h5file.root.wavefunction.shape[0]):
#			psi.GetData()[:,:,i] = h5file.root.wavefunction[i,:,:]
#	finally:
#		h5file.close()

	matvec(psi, tmpPsi, 0, 0)
	eigv = psi.InnerProduct(tmpPsi)
	psiNorm = sqrt(abs(psi.InnerProduct(psi)))

	print "Energy expectation value = %s" % eigv
	print "||psi|| = %s" % psiNorm


def TestSolveOverlapSpeed():
	def timeIt(func):
		t1 = time.time()
		func()
		t2 = time.time()
		Print("  Function '%s' took %4.1f s." % (func.func_name, (t2-t1)))

	numSolves = 100

	Print("")
	Print("Now testing multiple S^-1 * psi...")
	pyprop.Redirect.Enable(silent=True)

	seed(0)

	conf = pyprop.Load("config_eigenvalues.ini")
	psi = pyprop.CreateWavefunction(conf)
	tmpPsi = psi.Copy()

	Print("  Size of wavefunction is: %s" % repr(psi.GetData().shape)) 

	#Calculate S^-1 * psi
	Print("  Performing %i solves..." % numSolves)
	def solve():
		for i in range(numSolves):
			psi.GetRepresentation().MultiplyOverlap(tmpPsi)
	timeIt(solve)
	
	#finish and cleanup
	pypar.barrier()
	pyprop.Redirect.Disable()
	Print("\n...done!")


def TestEpetraMatvecSpeed():
	numMatVecs = 500

	Print("")
	Print("Now testing Epetra matvec speed...")
	pyprop.Redirect.Enable(silent=True)

	#Test
	conf = pyprop.Load("config_propagation.ini")
	psi = pyprop.CreateWavefunction(conf)
	Print("  Size of wavefunction is: %s" % repr(psi.GetData().shape)) 

	#Setup problem
	Print("  Setting up propagator w/potentials...")
	prop = SetupProblem(config='config_propagation.ini')
	psi = prop.psi
	tmpPsi = psi.Copy()
	tmpPsi.Clear()

	Print("  Local size of wavefunction is: %s" % str(prop.psi.GetData().shape)) 
	Print("  Global size of wavefunction is: %s" % str(prop.psi.GetRepresentation().GetFullShape())) 

	#Get Epetra potential
	#pot = prop.Propagator.BasePropagator.PotentialList[1]
	Print("  Number of potentials: %s" % len(prop.Propagator.BasePropagator.PotentialList))

	#Calculate S^-1 * psi
	Print("  Performing %i matvecs..." % numMatVecs)
	def matvecs():
		for i in range(numMatVecs):
			#pot.MultiplyPotential(psi, tmpPsi, 0, 0)
			prop.Propagator.BasePropagator.MultiplyHamiltonianNoOverlap(psi, tmpPsi, 0, 0)
			#tmpPsi.GetRepresentation().SolveOverlap(tmpPsi)

	timeIt(matvecs)
	
	#finish and cleanup
	pypar.barrier()
	pyprop.Redirect.Disable()
	pyprop.PrintOut("\n...done!")


def TestPropagation(inputFile = "groundstate_propagation.h5", **args):
	pyprop.PrintMemoryUsage("Before TestPropagation")
	
	Print("")
	Print("Now testing propagation...")
	pyprop.Redirect.Enable(silent=True)

	#Set up propagation problem
	potList = []
	#if not args.get("laserOff", False):
	#	Print("Setting up new problem with laser potentials...")
	#	potList += ["LaserPotentialVelocityDerivativeR1", "LaserPotentialVelocityDerivativeR2", "LaserPotentialVelocity"]
	#else:
	#	Print("Setting up new problem WITHOUT laser potentials...")

	#if not args.get("absorberOff", False):
	#	potList += ["Absorber"]
	#	Print("Setting up new problem with absorber...")
	#else:
	#	Print("Setting up new problem WITHOUT absorber...")

	pyprop.PrintMemoryUsage("Before SetupProblem")
	args["config"] = args.get("config", "config_propagation_nonorthdistr.ini")
	#args["config"] = args.get("config", "config_propagation.ini")
	prop = SetupProblem(**args)
	Print("Proc %i has wavefunction shape = %s" % (pyprop.ProcId, list(prop.psi.GetData().shape)), [pyprop.ProcId])
	pyprop.PrintMemoryUsage("After SetupProblem")

	#Load initial state
	Print("Loading intial state...")
	pyprop.PrintMemoryUsage("Before Loading InitialState")
	prop.psi.Clear()
	pyprop.serialization.LoadWavefunctionHDF(inputFile, "/wavefunction", prop.psi)
	prop.psi.Normalize()
	initPsi = prop.psi.Copy()
	initialEnergyCalculated = prop.GetEnergyExpectationValue()
	Print("Initial State Energy = %s" % initialEnergyCalculated)
	pyprop.PrintMemoryUsage("After Loading InitialState")

	Print("Done setting up problem!")
	
	#Propagate
	Print("Starting propagation")
	outputCount = args.get("outputCount", 10)
	startTime = time.time()
	for step, t in enumerate(prop.Advance(outputCount)):
		#calculate values
		norm = prop.psi.GetNorm()
		corr = abs(initPsi.InnerProduct(prop.psi))**2

		#estimate remaining time
		curTime = time.time() - startTime
		totalTime = (curTime / t) * prop.Duration
		eta = totalTime - curTime

		#Print stats
		Print("t = %.2f; N = %.15f; Corr = %.10f, ETA = %s" % (t, norm, corr, FormatDuration(eta)))
		Print("GMRES errors = %s" % prop.Propagator.Solver.GetErrorEstimateList())
		pyprop.PrintMemoryUsage("At t = %s" % t)
	
	endTime = time.time()
	Print("Propagation time = %s" % FormatDuration(endTime-startTime))

	#Final output
	pyprop.PrintMemoryUsage("After Propagation")
	norm = prop.psi.GetNorm()
	corr = abs(initPsi.InnerProduct(prop.psi))**2

	Print("Final status: N = %.10f; Corr = %.10f" % (norm, corr))

	Print("")
	prop.Propagator.Solver.PrintStatistics()

	#finish and cleanup
	pypar.barrier()
	pyprop.Redirect.Disable()
	pyprop.PrintMemoryUsage("After TestPropagation")
	Print("\n...done!")


def TestPreconditioner():
    prop = SetupProblem(config = "config_propagation_nonorthdistr.ini")
    prec = TwoElectronPreconditioner(prop.psi)
    prec.ApplyConfigSection(prop.Config.Preconditioner)
    dt = prop.Config.Propagation.timestep
    prec.SetHamiltonianScaling(1.0j * dt / 2.0)
    prec.SetOverlapScaling(1.0)
    prec.Setup(prop.Propagator)

    return prec


if __name__ == "__main__":
	#TestMultiplyOverlap()
	#TestSolveOverlap()
	#TestInnerProduct()
	#TestFindEigenvalues()
	#TestEpetraMatrix()
	#TestSolveOverlapSpeed()
	#TestEpetraMatvecSpeed()
	TestPropagation()
