#Import system modules
import sys
import os
from pylab import *
from numpy import *
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
	solver.Solve()

	if pyprop.ProcId == 0:
		for i in range(size(solver.GetEigenvalues())):
			prop.psi.GetData()[:] = solver.GetEigenvector(i)[:]
			SaveWavefunction("eigenvector%02i.h5" % i, "/eigenvector%02i" % i, prop.psi)

		h5file = tables.openFile("eigenvalues.h5", "w")
		try:
			h5file.createArray("/", "eigenvalues", solver.GetEigenvalues()[:])
		finally:
			h5file.close()
	
	return solver

def SaveWavefunction(filename, dataset, psi):
	if pyprop.ProcId == 0:

		if type(filename) == str:
			if os.path.exists(filename):
				os.unlink(filename)

	pyprop.serialization.SaveWavefunctionHDF(filename, dataset, psi)

