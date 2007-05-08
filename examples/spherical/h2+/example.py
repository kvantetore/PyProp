#improt system modules
import sys
import os

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("../../.."))

#Load and reload pyprop in order to get recent changes
import pyprop
#pyprop = reload(pyprop)
from libh2p import *

#numpy an pylab for good measure
from pylab import *
from numpy import *

#Choose radial grid type:
class RadialGridType:
	CARTESIAN = 1
	TRANSFORMED = 2

def SetupProblem(**args):
	#load config file. hydrogen.ini uses ../sphericalbase.ini as a base
	#configuration file, so be sure to check out that one as well.
	conf = SetupConfig(**args)

	#Uses the same Problem class, only changes the propagator specified in the
	#config file
	prop = pyprop.Problem(conf)

	#Set up all transformations and potentials.
	prop.SetupStep()

	return prop

def FindGroundstate(**args):
	#for finding ground state, we always use imaginary time
	args['imTime'] = True  #Remember: Arguments are case sensitive

	#Create the propagator
	prop = SetupProblem(**args)

	#propagate
	silent = False
	if 'silent' in args: silent = args['silent']
	for t in prop.Advance(20):
		if not silent:
			print "t = ", t, ", E = ", prop.GetEnergy()

	return prop

def FindGroundstateEnergy(**args):
	prop = FindGroundstate(**args)
	return prop.GetEnergy()

def PropagateState(**args):
	#Never use imtime for propagation
	args['imTime'] = False

	#Create propagator
	prop = SetupProblem(**args)

	if not 'initPsi' in args:
		raise "Argument initPsi not supplied"
	initPsi = args['initPsi']
	if initPsi.__class__ == str:
		prop.LoadWavefunction(initPsi)
		initPsi = prop.psi.Copy()

	else:
		#if initpsi and prop.psi has different number of 
		#spherical harmonics, it is no problem: we only load
		#the first n, where n is the smallest number of spherical
		#harmonics in initpsi and prop.psi
		lmlength1 = len(initPsi.GetData()[0,:])
		lmlength2 = len(prop.psi.GetData()[0,:])
		lmlength = min(lmlength1, lmlength2)
		prop.psi.GetData()[:,:lmlength] = initPsi.GetData()[:,:lmlength]

	for t in  prop.Advance(20):
		Corr = abs(prop.psi.InnerProduct(initPsi))**2
		Norm = prop.psi.GetNorm()
		print "t = ", t, ", Corr = ", Corr, ", Norm = ", Norm

	return prop


def FindConvergenceGridSizeCartesian(**args):
	args['gridType'] = RadialGridType.CARTESIAN	
	gridSize = r_[64:512:50]
	E = zeros(len(gridSize), dtype=double)
	for i in range(len(gridSize)):
		args['gridSize'] = gridSize[i]
		E[i] = FindGroundstateEnergy(**args)

	return gridSize, E
		
def FindNormLoss(**args):
	if 'dt' in args:
		dt = args['dt']
	else:
		dt = 0.01
	args['imTime'] = False

	if 'dtSteps' in args:
		dtSteps = args['dtSteps']
	else:
		dtSteps = 10

	dtList = r_[dt:dt/dtSteps:-dt/dtSteps]
	normLoss = zeros(len(dtList), dtype=double)
	for i in range(len(dtList)):
		dt = dtList[i]
		args['dt'] = dt
		prop = SetupProblem(**args)
		prop.psi.Normalize()
		n1 = prop.psi.GetNorm()
		for j in range(10): prop.AdvanceStep()
		n2 = prop.psi.GetNorm()
		normLoss[i] = (n2 / n1) ** ( 1 / (10 * dt))
		
	return dtList, normLoss	
#---------------------------------------------------------------------------------

def SetupConfig(**args):
	"""
	Loads a config file, and modifies it with the arguments sendt to this function 
	i.e. conf = SetupConfig(configFile='h2p.ini', dt=0.04, imTime=True), 
	will return a conf object loading h2p.ini, setting timestep to -0.04mj, and 
	renormalization to True.
	"""

	#config file
	if 'configFile' in args:
		configFile = args['config']
	else:
		configFile = "h2p.ini"
	conf = pyprop.Load(configFile)
	
	#Grid size
	if 'gridSize' in args:
		gridSize = args['gridSize']
		print "Using Grid Size = ", gridSize
		#set grid size of both representations, but we will
		#only use one of them in a given run
		conf.CartesianRadialRepresentation.rank0[2] = gridSize
		conf.TransformedRadialRepresentation.n = gridSize

	#Grid type
	if 'gridType' in args:
		gridType = args['gridType']
	else:
		gridType = None
	SetRadialGridType(conf, gridType)
	
	#TimeStep
	if 'dt' in args:
		dt = args['dt']
		print "Using TimeStep = ", dt
		conf.Propagation.timestep = abs(dt)

	#Imaginary Time Propagation
	if 'imTime' in args:
		imTime = args['imTime']
		print "Using ImaginaryTime = ", imTime
		dt = conf.Propagation.timestep
		if imTime:
			conf.Propagation.renormalization = True
			conf.Propagation.timestep = -1.0j * abs(dt)
		else:
			conf.Propagation.renormalization = False
			conf.Propagation.timestep = abs(dt)

	if 'softing' in args:
		softing = args['softing']
		print "Using softing ", softing
		conf.Potential.softing = softing

	if 'lmax' in args:
		lmax = args['lmax']
		print "Using LMax = ", lmax
		conf.AngularRepresentation.maxl = lmax

	if 'duration' in args:
		duration = args['duration']
		conf.Propagation.duration = duration

	if 'orientation' in args:
		orientation = args['orientation']
		conf.Potential.nuclear_orientation = orientation
		

	print "Setup Config Complete"
	print ""
	return conf

def SetRadialGridType(conf, gridType):
	"""
	Chooses between Cartesian and Transformed grid. On a loaded config
	object.
	"""

	if gridType == None:
		gridType = RadialGridType.TRANSFORMED

	#Set RadialRepresentation section
	if gridType == RadialGridType.CARTESIAN:
		conf.RadialRepresentation = conf.CartesianRadialRepresentation
		conf.RadialPropagator.propagator = pyprop.CartesianRadialPropagator

	elif gridType == RadialGridType.TRANSFORMED:
		conf.RadialRepresentation = conf.TransformedRadialRepresentation
		conf.RadialPropagator.propagator = pyprop.TransformedRadialPropagator

	else:
		raise "Invalid grid type ", gridType

	#Set radialtype 
	conf.Representation.radialtype = conf.RadialRepresentation.type

def MapLmIndex(l, m):
	return (l + 1) * l + m

def CalculateAngularMomentumDistribution(prop):
	radialRepr = prop.psi.GetRepresentation().GetRadialRepresentation()
	if radialRepr.__class__ == pyprop.core.CartesianRepresentation_1:
		weights = radialRepr.GetRange(0).Dx
	elif radialRepr.__class__ == pyprop.core.TransformedRadialRepresentation:
		weights = radialRepr.Range.GetWeights()	
	else:
		raise "Invaliud representation: " + str(radialRepr)

	lmax = prop.Propagator.LmRepresentation.Range.MaxL

	data = prop.psi.GetData()
	lDistrib = zeros(lmax+1, dtype=double)
	for l in r_[0:lmax+1]:
		for m in r_[-l:l+1]:
			idx = MapLmIndex(l,m)
			distrib = sum(weights * abs(data[:,idx])**2)
			lDistrib[l] += distrib

	return lDistrib

def CalculateMDistribution(prop):
	radialRepr = prop.psi.GetRepresentation().GetRadialRepresentation()
	if radialRepr.__class__ == pyprop.core.CartesianRepresentation_1:
		weights = radialRepr.GetRange(0).Dx
	elif radialRepr.__class__ == pyprop.core.TransformedRadialRepresentation:
		weights = radialRepr.Range.GetWeights()	
	else:
		raise "Invaliud representation: " + str(radialRepr)

	lmax = prop.Propagator.LmRepresentation.Range.MaxL

	data = prop.psi.GetData()
	mDistrib = zeros(lmax+1, dtype=double)
	for l in r_[0:lmax+1]:
		for m in r_[-l:l+1]:
			idx = MapLmIndex(l,m)
			distrib = sum(weights * abs(data[:,idx])**2)
			mDistrib[abs(m)] += distrib

	return mDistrib
	
def CalculateLMDistribution(prop):
	radialRepr = prop.psi.GetRepresentation().GetRadialRepresentation()
	if radialRepr.__class__ == pyprop.core.CartesianRepresentation_1:
		weights = radialRepr.GetRange(0).Dx
	elif radialRepr.__class__ == pyprop.core.TransformedRadialRepresentation:
		weights = radialRepr.Range.GetWeights()	
	else:
		raise "Invaliud representation: " + str(radialRepr)

	lmax = prop.Propagator.LmRepresentation.Range.MaxL

	data = prop.psi.GetData()
	mDistrib = zeros((lmax+1)**2, dtype=double)
	for l in r_[0:lmax+1]:
		for m in r_[-l:l+1]:
			idx = MapLmIndex(l,m)
			distrib = sum(weights * abs(data[:,idx])**2)
			mDistrib[idx] = distrib

	return mDistrib
	
def PlotAngularMomentumDistribution(prop):
	lDistrib = CalculateAngularMomentumDistribution(prop)
	bar(r_[0:len(lDistrib)], lDistrib)

	xlabel("Angular Momentum (l)")
	ylabel("Probability")
	title("Angular momentum distribution")
	return lDistrib


