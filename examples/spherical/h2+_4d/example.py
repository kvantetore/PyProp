#improt system modules
import sys
import os

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("../../.."))

#Load and reload pyprop in order to get recent changes
import pyprop
pyprop = reload(pyprop)
from libh2p import *

#numpy an pylab for good measure
from pylab import *
from numpy import *

#Choose radial grid type:
class GridType:
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
	
	#Radial Grid size
	if 'radialGridSize' in args:
		radialGridSize = args['radialGridSize']
		print "Using Grid Size = ", radialGridSize
		#set grid size of both representations, but we will
		#only use one of them in a given run
		conf.CartesianRadialRepresentation.rank0[2] = radialGridSize
		conf.TransformedRadialRepresentation.n = radialGridSize

	#Radial Grid type
	if 'radialGridType' in args:
		radialGridType = args['radialGridType']
	else:
		radialGridType = None
	SetRadialGridType(conf, radialGridType)

	#Nuclear Grid size
	if 'nuclearGridSize' in args:
		nuclearGridSize = args['nuclearGridSize']
		print "Using Nuclear Grid Size = ", nuclearGridSize
		#set grid size of both representations, but we will
		#only use one of them in a given run
		conf.CartesianNuclearRepresentation.rank0[2] = nuclearGridSize
		conf.TransformedNuclearRepresentation.n = nuclearGridSize

	#Nuclear Grid type
	if 'nuclearGridType' in args:
		nuclearGridType = args['nuclearGridType']
	else:
		nuclearGridType = None
	SetNuclearGridType(conf, nuclearGridType)
	
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

def SetNuclearGridType(conf, gridType):
	"""
	Chooses between Cartesian and Transformed grid. On a loaded config
	object.
	"""
	if gridType == None:
		gridType = GridType.TRANSFORMED

	#Set RadialRepresentation section
	if gridType == GridType.CARTESIAN:
		conf.NuclearRepresentation = conf.CartesianNuclearRepresentation
		conf.NuclearPropagator.propagator = pyprop.CartesianRadialPropagator

	elif gridType == GridType.TRANSFORMED:
		conf.NuclearRepresentation = conf.TransformedNuclearRepresentation
		conf.NuclearPropagator.propagator = pyprop.TransformedRadialPropagator

	else:
		raise "Invalid nuclear grid type ", gridType

def SetRadialGridType(conf, gridType):
	"""
	Chooses between Cartesian and Transformed grid. On a loaded config
	object.
	"""
	if gridType == None:
		gridType = GridType.TRANSFORMED

	#Set RadialRepresentation section
	if gridType == GridType.CARTESIAN:
		conf.RadialRepresentation = conf.CartesianRadialRepresentation
		conf.RadialPropagator.propagator = pyprop.CartesianRadialPropagator

	elif gridType == GridType.TRANSFORMED:
		conf.RadialRepresentation = conf.TransformedRadialRepresentation
		conf.RadialPropagator.propagator = pyprop.TransformedRadialPropagator

	else:
		raise "Invalid radial grid type ", gridType

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


