#!/bin/env python

import sys
sys.path.append("./pyprop")
sys.path.append("../helium_stabilization")
import pyprop
import numpy
from scipy import *
from pylab import *

import pypar

from libepetratest import *
from libpotential import *

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


if __name__ == "__main__":
	TestMultiplyOverlap()
	TestSolveOverlap()
	TestInnerProduct()
	TestFindEigenvalues()
	#TestEpetraMatrix()
