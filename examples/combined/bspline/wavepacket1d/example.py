#Import system modules
import sys
import os
from numpy import double, r_, fft
from numpy import max as nmax
from numpy import max as nmin
from numpy import where as nwhere
from time import sleep

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
	if "imtime" in args:
		imtime = args["imtime"]
		propSection = conf.Propagation
		dt = abs(propSection.timestep)
		renormalize = False
		if imtime:
			dt = -1.0j * dt
			renormalize = True
		else:
			conf.DynamicPotential.frequency = eval(conf.DynamicPotential.frequency)

		propSection.timestep = dt
		propSection.renormalization = renormalize

	if "amplitude" in args:
		amplitude = args["amplitude"]
		conf.DynamicPotential.amplitude = amplitude

	if "xmax" in args:
		xMax = args['xmax']
		conf.BSplineRepresentation.xmax = xMax
	
	if "xsize" in args:
		xSize = args['xsize']
		conf.BSplineRepresentation.xsize = xSize
	
	if "dt" in args:
		timeStep = args['dt']
		conf.Propagation.timestep = timeStep

	if "pulseDuration" in args:
		T = args['pulseDuration']
		conf.Propagation.duration = T
		DynamicPotential.pulse_duration = T

	if "lmax" in args:
		lmax = args['lmax']
		conf.AngularRepresentation.lmax = lmax
		conf.BSplinePropagator.lmax = lmax
		
	return conf


def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop


def FindGroundstate(**args):
	args['imtime'] = True
	prop = SetupProblem(**args)
	
	for t in prop.Advance(5):
		pass
		E = prop.GetEnergy()
		if pyprop.ProcId == 0:
			print "t = %f, E = %f" % (t, E)

	E = prop.GetEnergy()
	if pyprop.ProcId == 0:
		print "Ground State Energy = %f" % E

	if "savewavefunction" in args:
		if args["savewavefunction"]:
			SaveWavefunction("groundstate.h5", "/wavefunction", prop.psi)

	return prop


def SetupInitialState(**args):

	if 'stateIndex' in args:
		stateIndex = args['stateIndex']
	else:
		stateIndex = 0
		args['stateIndex'] = stateIndex

	#Set up problem
	conf = SetupConfig(**args)

	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#Find eigenstates of desired l-subspace
	M = GetHamiltonMatrix(prop)
	E,V = eig(M)
	I = argsort(E)
	print "Initial state energy = ", E[I[stateIndex]].real

	#Assign initial state
	prop.psi.GetData()[:] = 0
	prop.psi.GetData()[:] = V[:,I[stateIndex]]
	prop.psi.Normalize()

	return prop


def SaveWavefunction(filename, dataset, psi):
	if pyprop.ProcId == 0:
		if os.path.exists(filename):
			os.unlink(filename)
	pyprop.serialization.SaveWavefunctionHDF(filename, dataset, psi)


def FindIonizationProbability(**args):
	prop = SetupInitialState(**args)

	initialPsi = prop.psi.Copy()

	prop.norm = []
	prop.corr = []
	prop.times =  []
	for t in prop.Advance(200):
		prop.times += [t]
		prop.norm += [prop.psi.GetNorm()]
		prop.corr += [abs(prop.psi.InnerProduct(initialPsi))**2]
		print "t = %f, Norm = %f, Corr = %f" % (t, prop.norm[-1], prop.corr[-1])

	norm = prop.psi.GetNorm()
	corr = abs(prop.psi.InnerProduct(initialPsi))**2
	print "Ionization Probability = %f" % norm
	print "Initial state correlation = %f" % corr

	return prop


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


def GetHamiltonMatrix(prop):
	size = prop.psi.GetData().size
	matrix = zeros((size, size), dtype=complex)
	tempPsi = prop.GetTempPsi()

	for i in range(size):
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[i] = 1

		tempPsi.GetData()[:] = 0
		prop.MultiplyHamiltonian(tempPsi)
		
		matrix[:, i] = tempPsi.GetData()[:]
		
	return matrix


def PropagateWavePacket(**args):
	#Set up problem
	prop = SetupProblem(**args)
	conf = prop.Config

	#Setup traveling wavepacket initial state
	f = lambda x: conf.Wavepacket.function(conf.Wavepacket, x)
	bspl = prop.psi.GetRepresentation().GetRepresentation(0).GetBSplineObject()
	c = bspl.ExpandFunctionInBSplines(f)
	prop.psi.GetData()[:] = c
	prop.psi.Normalize()
	initialPsi = prop.psi.Copy()

	#Get x-grid
	subProp = prop.Propagator.SubPropagators[0]
	subProp.InverseTransform()
	grid = prop.psi.GetRepresentation().GetLocalGrid(0)
	subProp.ForwardTransform()

	#Setup equispaced x grid
	x_min = conf.BSplineRepresentation.xmin
	x_max = conf.BSplineRepresentation.xmax
	grid_eq = linspace(x_min, x_max, grid.size)
	x_spacing = grid_eq[1] - grid_eq[0]
	
	#Set up fft grid
	k_spacing = 1.0 / grid_eq.size
	k_max = pi / x_spacing
	k_min = -k_max
	k_spacing = (k_max - k_min) / grid_eq.size
	grid_fft = zeros(grid_eq.size, dtype=double)
 	grid_fft[:grid_eq.size/2+1] = r_[0.0:k_max:k_spacing]
   	grid_fft[grid_eq.size/2:] = r_[k_min:0.0:k_spacing]
	print "Momentum space resolution = %f a.u." % k_spacing
	
	k0 = conf.Wavepacket.k0
	k0_trunk = 5 * k0
	trunkIdx = list(nwhere(abs(grid_fft) <= k0_trunk)[0])

	rcParams['interactive'] = True
	figure()
	p1 = subplot(211)
	p2 = subplot(212)
	p1.hold(False)

	psi_eq = zeros((grid.size), dtype=complex)

	for t in prop.Advance(40):
		print "t = %f, norm = %.15f, P = %.15f " % \
			( t, prop.psi.GetNorm(), abs(prop.psi.InnerProduct(initialPsi))**2 )
		sys.stdout.flush()

		subProp.InverseTransform()
		p1.plot(grid, abs(prop.psi.GetData())**2)
		subProp.ForwardTransform()
		bspl.ConstructFunctionFromBSplineExpansion(prop.psi.GetData(), grid_eq, psi_eq)
		psi_fft = (abs(fft.fft(psi_eq))**2)
		psi_fft_max = nmax(psi_fft)
		psi_fft /= psi_fft_max

		#Plot momentum space |psi|**2
		p2.hold(False)
		p2.semilogy(grid_fft[trunkIdx], psi_fft[trunkIdx] + 1e-21)
		p2.hold(True)
		p2.semilogy([0,0], [1e-20, psi_fft_max], 'r-')
		p2.semilogy([-k0,-k0], [1e-20, psi_fft_max], 'g--')
		p2.semilogy([k0,k0], [1e-20, psi_fft_max], 'g--')

		#Set subplot axis
		p1.axis([grid[0],grid[-1],0,0.1])
		p2.axis([-k0_trunk, k0_trunk, 1e-20, psi_fft_max])
		show()

	hold(True)
		

	return prop
