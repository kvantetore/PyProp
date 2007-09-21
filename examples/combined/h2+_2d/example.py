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
	for t in prop.Advance(10):
		if not silent:
			print "t = ", t, ", E = ", prop.GetEnergyExpectationValue()

	return prop

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



