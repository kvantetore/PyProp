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

import potential_data
import spline

#Physical constants
eV_to_au = 0.036764706
au_to_eV = 27.2114  

inverse_cm_to_au = 4.5563276e-06
au_to_inverse_cm = 219475.

volt_per_metre_to_au = 1.944689e-12
femtosec_to_au = 41.36338

def ReducedMass(M1, M2):
	return M1 * M2 / (M1 + M2)

mass_proton_au = 1836.152666
mass_deuteron_au = 3670.482954
hdp_reduced_mass_au = ReducedMass(mass_proton_au, mass_deuteron_au)

def wavelength_to_frequency(nm):
	return inverse_cm_to_au/(nm*1e-7)

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

	if 'species' in args:
		species = args["species"]
		conf.VibrationalPotential.species = species

	if "molecule" in args:
		molecule = args["molecule"].lower()
		if molecule == "h2+":
			conf.RadialPropagator.mass = ReducedMass(mass_proton_au, mass_proton_au)
			conf.ElectronicCoupling.static_dipole_moment = 0
		elif molecule == "d2+":
			conf.RadialPropagator.mass = ReducedMass(mass_deuteron_au, mass_deuteron_au)
			conf.ElectronicCoupling.static_dipole_moment = 0
		elif molecule == "hd+":
			conf.RadialPropagator.mass = ReducedMass(mass_deuteron_au, mass_proton_au)
			conf.ElectronicCoupling.static_dipole_moment = 1. / 3.
		else:
			raise Exception("Unknown molecule '%s'" % molecule)

	if 'radialScaling' in args:
		radialScaling = args["radialScaling"]
		rank0 = conf.RadialRepresentation.rank0
		radialRange = rank0[1] - rank0[0]
		radialSize = rank0[2]
	
		rank0[1] = radialRange / radialScaling + rank0[0]
		rank0[2] = radialSize / radialScaling

	if 'pulseDelay' in args:
		pulseDelay = args["pulseDelay"]
		conf.ElectronicCoupling.delay = pulseDelay

	if "duration" in args:
		duration = args["duration"]
		conf.Propagation.duration = duration
		
	return conf

def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop

def FindGroundstate(**args):
	args['config'] = "groundstate.ini"
	args['imtime'] = True
	prop = SetupProblem(**args)
	
	for t in prop.Advance(10):
		E = prop.GetEnergy()
		print "t = %f, E = %.17f" % (t/femtosec_to_au, E)

	E = prop.GetEnergy()
	print "Ground State Energy = %f" % E

	return prop

def SetupEigenstates(**args):
	args["config"] = "groundstate.ini"
	args["species"] = "Ion"
	outputfile = args["outputfile"]

	solver = FindEigenvalues(**args)
	E = solver.GetEigenvalues()
		
	prop = solver.BaseProblem
	shape = (len(E), ) + prop.psi.GetData().shape

	pyprop.serialization.RemoveExistingDataset(outputfile, "/eigenstates")
	pyprop.serialization.RemoveExistingDataset(outputfile, "/eigenvalues")

	f = tables.openFile(outputfile, "a")
	try:
		f.createArray(f.root, "eigenvalues", E)
		atom = tables.ComplexAtom(itemsize=16)
		eigenstates = f.createCArray(f.root, "eigenstates", atom, shape)
		
		for i in range(len(E)):
			solver.SetEigenvector(prop.psi, i)
			eigenstates[i,:] = prop.psi.GetData()

	finally:
		f.close()

def SetupInitialState(**args):
	outputfile = args["outputfile"]

	#Setup initial state
	prop = FindGroundstate(**args)
	prop.SaveWavefunctionHDF(outputfile, "/initial_state")


def GetInputFile(**args):
	if not "molecule" in args:
		raise Exception("Please specify molecule")

	conf = SetupConfig(**args)
	molecule = args["molecule"]
	gridSize = conf.RadialRepresentation.rank0[1]
	boxSize = conf.RadialRepresentation.rank0[2]

	return "inputfiles/%s_%i_%i" % (molecule, gridSize, boxSize)


def SetupInputFile(**args):
	outputfile = GetInputFile(**args)
	args["outputfile"] = outputfile
	SetupEigenstates(**args)
	SetupInitialState(**args)


def Propagate(**args):
	args['config'] = "config.ini"
	args['imtime'] = False
	inputfile = GetInputFile(**args)
	prop = SetupProblem(**args)

	f = tables.openFile(inputfile, "r")
	try:
		#Setup initial state
		prop.psi.GetData()[:,0] = f.root.initial_state[:]
		prop.psi.GetData()[:,1] = 0
		initPsi = prop.psi.Copy()

		#Load eigenstates
		eigenstates = conj(f.root.eigenstates[:])
		f.close()
	finally:
		f.close()
	del f

	if len(eigenstates.shape) != 2:
		raise Exception("Please implement for Rank!=2")
	eigenstateSize = len(eigenstates[0,:])
	weight = prop.psi.GetRepresentation().GetRepresentation(0).GetRange(0).Dx

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	timeList = []
	corrList = []
	radialData = []

	tPrev = 0
	for t in prop.Advance(True):
		psiSlice = prop.psi.GetData()[0:eigenstateSize,0]
		corrList += [abs(dot(eigenstates, psiSlice)*weight)**2]
		radialData += [sum(abs(prop.psi.GetData())**2, axis=1)]
		
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		timeList.append(t/femtosec_to_au)
		if t-tPrev > 100*femtosec_to_au:
			print "t = %f, N = %f, Corr = %.17f" % (t/femtosec_to_au, norm, corr) 
			tPrev = t
		#plot(r, abs(prop.psi.GetData())**2)

	timeList.append(prop.PropagatedTime / femtosec_to_au)
	psiSlice = prop.psi.GetData()[0:eigenstateSize,0]
	corrList += [abs(dot(eigenstates, psiSlice)*weight)**2]

	norm = prop.psi.GetNorm()
	print "t = %f, N = %f" % (t, norm)

	return array(timeList), array(corrList), array(r), array(radialData)

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


def GetVibrationalPotential(psi, config, potential):
	#Decide which potential	to use
	species = config.species.lower()
	if species == "neutral":	
		grid = potential_data.GridNeutral
		data1 = potential_data.PotentialNeutral
		data2 = potential_data.PotentialNeutral
	elif species == "ion":
		grid = potential_data.GridIon
		data1 = potential_data.Potential_1s_sigma_G
		data2 = potential_data.Potential_2p_sigma_U
	else:
		raise Exception("Unknown species '%s'" % config.species)

	#mirror round r=0
	grid = array(list(-grid[::-1]) + list(grid))
	data1 = array(list(data1[::-1]) + list(data1))
	data2 = array(list(data2[::-1]) + list(data2))

	#Interpolate the potential
	interp1 = spline.Interpolator(grid, data1)
	interp2 = spline.Interpolator(grid, data2)

	r = psi.GetRepresentation().GetLocalGrid(0)
	dr = r[1] - r[0]
	if psi.GetRank() == 1:
		potential[:] = [interp1.Evaluate(x) +  minimum(2./dr, 1./ abs(x)) for x in r]
	
	if psi.GetRank() == 2:
		potential[:,0] = [interp1.Evaluate(x) +  minimum(2./dr, 1./ abs(x)) for x in r]
		potential[:,1] = [interp2.Evaluate(x) +  minimum(2./dr, 1./ abs(x)) for x in r]

	#plot(r, potential)

def GetElectronicCouplingBates(psi, config):
	r = psi.GetRepresentation().GetLocalGrid(0)
	rho = (1 + abs(r) + r**2/3.) * exp(-abs(r)) 

	coupling = -1/(2 + 1.4*abs(r)) + r/(2*sqrt(1 - rho**2))
	coupling[where(r==0)] = 0
	return coupling


def GetElectronicCouplingBunkinTugov(psi, config):
	r = psi.GetRepresentation().GetLocalGrid(0)
	r_zero = 2.
	alpha = 0.72
	d = 1.07
	d_prime = 0.396
	xp = -0.055
	
	return d + d_prime*(abs(r)-r_zero) - alpha*xp*d_prime*(abs(r)-r_zero)**2

def GetElectronicCoupling(psi, config):
	if config.coupling_model == "Bates":
		coupling = GetElectronicCouplingBates(psi, config)
	elif config.coupling_model == "BunkinTugov":
		coupling = GetElectronicCouplingBunkinTugov(psi, config)
	else:
		raise "Invalid coupling model %s" % config.coupling_model

	shape = psi.GetData().shape
	if psi.GetRank() == 1:
		raise Exception("Must be a 2D wavefunction in order to use the coupling")
	if shape[1] != 2:
		raise Exception("Electronic coupling rank (2) must be of size 2")

	r = psi.GetRepresentation().GetLocalGrid(0)
	matrix = zeros((shape[0], 2, 2), dtype=complex)
	matrix[:,0,0] = config.static_dipole_moment * abs(r) / 2.0
	matrix[:,0,1] = coupling
	matrix[:,1,0] = conj(coupling)
	matrix[:,1,1] = config.static_dipole_moment * abs(r) / 2.0

	return matrix
		

def GetLaserField(config, t):
	E0 = 2744 * volt_per_metre_to_au * sqrt(config.intensity) 
	t0 = config.delay
	scaling = 4 * log(2) 

	envelope = exp( - scaling * (t - t0)**2 / config.duration**2 )
	field = E0 * cos(config.frequency*(t - t0)) * envelope

	return field
	
	
def PlotEigenstate(solver, index):
	prop = solver.BaseProblem
	solver.SetEigenvector(prop.psi, index)

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	data = prop.psi.GetData()
	plot(r, abs(data)**2 -0.6)
