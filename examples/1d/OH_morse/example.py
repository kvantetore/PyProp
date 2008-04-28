#Import system modules
import sys
import os
from pylab import *
import numpy
import tables

#Load pyprop
sys.path.insert(1, os.path.abspath("./pyprop"))
import pyprop
pyprop = reload(pyprop)

#Load the project module
from libpotential import *

#Load unit conversion tools
execfile('constants.py')

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
		if pyprop.ProcId == 0:
			print "t = %f, E = %f" % (t, E)

	E = prop.GetEnergy()
	if pyprop.ProcId == 0:
		print "Ground State Energy = %f" % E

	SaveWavefunction("groundstate.h5", "/wavefunction", prop.psi)

	return prop


def FindEigenstates(**args):
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	numEigs = prop.Config.Arpack.krylov_eigenvalue_count
	outFile = prop.Config.Arpack.hdf5_filename

	#Find eigenvectors
	solver.Solve()

	if pyprop.ProcId == 0:
		h5file = tables.openFile(outFile, "w")
		try:
			#Save eigenvectors
			for i in range(numEigs):
				prop.psi.GetData()[:] = solver.GetEigenvector(i)
				prop.psi.Normalize()
				h5file.createArray("/", "eigenvector%03i" % i, prop.psi.GetData())

			#Save eigenvalues
			h5file.createArray("/", "eigenvalues", solver.GetEigenvalues()[:])

			#Store config object
			h5file.setNodeAttr("/", "configObject", prop.Config.cfgObj)
		finally:
			h5file.close()
	
	return solver


def CalculateEigenBasisMatrixElements(prop, potential, eigFileName, **args):
	"""
	From eigenvectors of the 1D morse problem and a given 1D potential, 
	calculate matrix elements in the eigenfunction basis. Result is
	stored as a dense matrix in a HDF5 file.
	"""

	outFile = args.has_key("outFile") and args["outFile"] or "out/morse_matrix.h5"

	#Set up buffer
	tmpPsi = prop.psi.CopyDeep()

	eigFile = tables.openFile(eigFileName)
	
	numberOfEigs = eigFile.root.eigenvalues[:].size

	M = zeros((numberOfEigs, numberOfEigs), dtype=complex)
	for i in range(numberOfEigs):
		#Load eigenvector i
		prop.psi.GetData()[:] = eigFile.getNode("/eigenvector%03i" % i)
		for j in range(i, numberOfEigs):
			#Load eigenvector j
			tmpPsi.GetData()[:] = eigFile.getNode("/eigenvector%03i" % j)

			#Multiply potential:  X|j>
			potential.MultiplyPotential(prop.psi, tmpPsi, 0, 0)
			
			#Calculate <i|X|j>
			M[i,j] = prop.psi.InnerProduct(tmpPsi)
			if i != j:
				M[j,i] = numpy.conj(M[i,j])

	eigFile.close()

	#Store matrix
	h5file = tables.openFile(outFile, "w")
	try:
		h5file.createArray("/", "morse_matrix_elements", M)
	finally:
		h5file.close()


def MorseOscillatorEnergies(howMany, conf):
	"""
	
	E = v_0 * (v + 1/2) - (v_0 * (v + 1/2))**2 / (4 * D_e)

	v_0 = a/(2*pi*c)*sqrt(2*D_e/m)
	a = sqrt(k_e/(2*D_e))
	"""

	D_e = conf.MorsePotential.strength
	mass = conf.Propagation.mass
	theta = conf.MorsePotential.theta
	v_0 = theta * sqrt(2 * D_e / mass)
	E = [v_0 * (v + 0.5) - (v_0 * (v + 0.5))**2 / (4 * D_e) - D_e for v in range(howMany)]
	
	return E

def SetupPotential(confSection):
	pot = eval(confSection.classname + "_1()")
	pot.ApplyConfigSection(confSection)

	return pot

def SaveWavefunction(filename, dataset, psi):
	if pyprop.ProcId == 0:

		if type(filename) == str:
			if os.path.exists(filename):
				os.unlink(filename)

	pyprop.serialization.SaveWavefunctionHDF(filename, dataset, psi)


def GetDiagonalElements(psi, config, potential):
	"""
	A funtion to provide diagonal (energy) matrix elements
	"""
	h5file = tables.openFile(config.filename, "r")
	try:
		potential[:] = h5file.getNode(config.dataset)[:]
		potential[:] *= config.scaling	
	finally:
		h5file.close()


def ComputeTargetStateZhuRabitzExperiment(prop, numEigs=0, eigFile=None, \
	eigDataSet = "/eigenvector", outFileName = "zhu_rabitz_final_state.h5"):
	"""
	Zhu and Rabitz use a gaussian projection operator to characterize their target space,

	    P = gamma / sqrt(pi) * exp[-gamma**2 * (x - x')**2 ]

	This function computes |ZR> = sum(<i|P|i>|i>, i), where |i> is the i'th eigenfunction 
	of the Morse oscillator.
	"""

	#Setup Zhu-Rabtiz operator
	ZhuRabitzOperator = eval(prop.Config.ZhuRabitzOperator.classname + "_1()")
	ZhuRabitzOperator.ApplyConfigSection(prop.Config.ZhuRabitzOperator)

	#Get tmpPsi
	tmpPsi = prop.psi.Copy()

	#Create vector to hold Zhu-Rabitz state
	ZhuRabitzState = numpy.zeros(numEigs, dtype=complex)

	h5file = tables.openFile(eigFile, "r")

	for i in range(numEigs):
		tmpPsi.Clear()
		prop.psi.GetData()[:] = h5file.getNode("%s%03i" % (eigDataSet, i))
		ZhuRabitzOperator.MultiplyPotential(prop.psi, tmpPsi, 0, 0)
		ZhuRabitzState[i] = prop.psi.InnerProduct(tmpPsi)
		
	outFile = tables.openFile(outFileName, "w")
	try:
		outFile.createArray("/", "ZhuRabitzFinalState", ZhuRabitzState)
	finally:
		outFile.close()


def GetZhuRabitzOperator(conf, grid):
	h5file = tables.openFile(conf.filename, "r")
	potential = 0
	try:
		potential = h5file.getNode(conf.dataset)[:]
	finally:
		h5file.close()

	return potential

def SetupKrotov(config, **args):
	"""
	Setup Krotov problem
	"""

	conf = pyprop.Load(config)

	if "timestep" in args:
		config.Propagation.timestep = args["timestep"]

	prop = pyprop.Problem(conf)
	prop.SetupStep()
	krotov = pyprop.Krotov(prop)
	return krotov


def SetupZhuRabitz(config, **args):
	"""
	Setup ZhuRabitz problem
	"""

	conf = pyprop.Load(config)

	if "timestep" in args:
		config.Propagation.timestep = args["timestep"]

	prop = pyprop.Problem(conf)
	prop.SetupStep()
	zhurabitz = pyprop.ZhuRabitz(prop)
	return zhurabitz


def SetupDegani(config, **args):
	"""
	Setup Degani problem
	"""

	conf = pyprop.Load(config)

	prop = pyprop.Problem(conf)
	prop.SetupStep()
	degani = pyprop.Degani(prop)
	return degani
