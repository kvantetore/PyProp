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

#Load movie stuff
try:
	execfile('pyprop/pyprop/plot/makemovie.py')
	execfile('pyprop/pyprop/utilities/cartesian/rebuildwavefunction.py')
except:
	print "Could not load movie utilities"

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
	"""
	Use PIRAM to calculate some eigenstates. Then store these and the 
	corresponding energies to HDF5 files.
	"""

	#Setup problem and call PIRAM solver
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	solver.Solve()

	#Create output dir if it does not exist
	try:
		os.mkdir("eigenvectors")
	except:
		pass

	#Proc 0 saves the eigenvalues
	if pyprop.ProcId == 0:
		print solver.GetEigenvalues()
		h5file = tables.openFile("eigenvectors/eigenvalues.h5", "w")
		h5file.createArray("/", "eigenvalues", solver.GetEigenvalues())
		h5file.close()

	for i in range(size(solver.GetEigenvalues())):
		solver.SetEigenvector(prop.psi, i)
		SaveWavefunction("eigenvectors/eigenvector%02i.h5" % i, "/wavefunction%02i" % i, prop.psi)
	
	return solver

def SaveWavefunction(filename, dataset, psi):
	if pyprop.ProcId == 0:

		if type(filename) == str:
			if os.path.exists(filename):
				os.unlink(filename)

	pyprop.serialization.SaveWavefunctionHDF(filename, dataset, psi)


def IonizationMovie(**args):
	args['imtime'] = False
	args['config'] = "propagation.ini"
	prop = SetupProblem(**args)

	initialPsi = prop.psi.Copy()
	x = prop.psi.GetRepresentation().GetLocalGrid(0)
	y = prop.psi.GetRepresentation().GetLocalGrid(1)
	
	rcParams['interactive'] = False
	fig = figure(figsize=(7,7))

	for t in enumerate(prop.Advance(200)):
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initialPsi))**2
		print "t = %f, Norm = %f, Corr = %f" % (t[1], norm, corr)

		SaveWavefunction("out/propagation%03i" % t[0], "/wavefunction", prop.psi)

		imshow(abs(prop.psi.GetData())**2, extent=(x[0],x[-1],y[0],y[-1]))
		savefig('figs/wavefunction%03i.png' % t[0])

	rcParams['interactive'] = True
	close(fig)

	return prop

def IonizationMovie2(**args):
	args['imtime'] = False
	args['config'] = "propagation.ini"
	prop = SetupProblem(**args)

	#Setup wavefunction rebuilder
	rebuilder = WavefunctionRebuilderCartesian(prop)
	prop.Config.Rebuilder.Apply(rebuilder)
	rebuilder.Setup()

	#Setup movie maker
	moviemaker = MakeMovie(rebuilder)
	prop.Config.Movie.Apply(moviemaker)
	
	#Create frames
	moviemaker.CreateFrames()

	#Create movie
	moviemaker.CreateMovie()

	return moviemaker



def FindIonizationProbabilityAmplitude():
	amplitudeList = r_[0:2:0.2]
	ionizationList = zeros(len(amplitudeList), dtype=double)

	for i in range(len(amplitudeList)):
		prop = FindIonizationProbability(amplitude=amplitudeList[i])
		ionizationList[i] = 1.0 - prop.psi.GetNorm()

	plot(amplitudeList, ionizationList, label="Ionization Probability")
	xlabel("Electric Field Strength")
	ylabel("Ionization Probability")
	legend()

	return amplitudeList, ionizationList


def GetHamiltonMatrixSubspace(prop, l):
	size = prop.psi.GetData().shape[0]
	matrix = zeros((size, size), dtype=complex)
	tempPsi = prop.GetTempPsi()

	for i in range(size):
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[i,l] = 1

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		
		matrix[:, i] = tempPsi.GetData()[:,l]
		
	return matrix


def TestSoftParameter():

	conf = SetupConfig()

	E = []
	softParams = [0.005, 0.01, 0.02, 0.04, 0.06]

	for s in softParams:
		conf.TwoElectronCorrelation.soft_param = s
		prop = pyprop.Problem(conf)
		prop.SetupStep()

		for t in prop.Advance(True):
			pass

		E.append(prop.GetEnergyExpectationValue())

	return E
