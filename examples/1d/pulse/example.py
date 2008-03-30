#improt system modules
import sys
import os

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("../../.."))

#Load and reload pyprop in order to get recent changes
import pyprop
pyprop = reload(pyprop)
from libpotential import *

#numpy an pylab for good measure
try:
	from pylab import *
except: pass
from numpy import *

#Choose radial grid type:
class GridType:
	CARTESIAN = 1
	TRANSFORMED = 2

def SetGridType(conf, gridType):
	if gridType == None:
		print "Grid type not specified, using cartesian grid"
		gridType = GridType.CARTESIAN

	#Set Representation section
	if gridType == GridType.CARTESIAN:
		conf.Representation = conf.CartesianRepresentation
		conf.Propagation.propagator = pyprop.CartesianPropagator

	elif gridType == GridType.TRANSFORMED:
		conf.Representation = conf.TransformedRepresentation
		conf.Propagation.propagator = pyprop.TransformedGridPropagator
	else:
		raise "Invalid grid type ", gridType

	#Set radialtype 
	conf.Representation.radialtype = conf.Representation.type

def SetupConfig(conf, **args):
	#Grid size
	if 'gridSize' in args:
		gridSize = args['gridSize']
		print "Using Grid Size = ", gridSize
		#set grid size of both representations, but we will
		#only use one of them in a given run
		conf.CartesianRepresentation.rank0[2] = gridSize
		conf.TransformedRepresentation.n = gridSize

	if 'silent' in args:
		conf.Propagation.silent = args['silent']

	#Grid type
	if 'gridType' in args:
		gridType = args['gridType']
	else:
		gridType = None
	SetGridType(conf, gridType)
	
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
	print "Setup Config Complete"
	print ""

def test(**args):
	conf = pyprop.Load("config.ini")
	SetupConfig(conf, **args)

	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop

def MapLmIndex(l, m):
	return (l + 1) * l + m

def FindGroundstate(**args):
	conf = pyprop.Load("groundstate.ini")
	SetupConfig(conf, **args)

	prop = pyprop.Problem(conf)
	prop.SetupStep()

	for t in prop.Advance(5):
		print "t = ", t, " E = ", prop.GetEnergyImTime()

	return prop

def RunPulseExperiment(gridType=GridType.CARTESIAN):
	#Find groundstate
	initPsi = FindGroundstate(gridType=gridType)

	#load config file. hydrogen.ini uses ../sphericalbase.ini as a base
	#configuration file, so be sure to check out that one as well.
	conf = pyprop.Load("config.ini")
	SetGridType(conf, gridType)

	#Uses the same Problem class, only changes the propagator specified in the
	#config file
	prop = pyprop.Problem(conf)

	#Set up all transformations and potentials.
	prop.SetupStep()

	#load wavefunction from groundstate
	prop.psi.GetData()[:] = initPsi.GetData()

	#propagate through the problem
	index = 0
	for t in prop.Advance(10):
		pyprop.Plot1D(prop)
		ax = axis()
		axis([ax[0], 10, ax[2], ax[3]])
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		print "t = ", t, "; Norm = ", norm, "; Corr = ", corr
		
		index += 1

	norm = prop.psi.GetNorm()
	corr = abs(prop.psi.InnerProduct(initPsi))**2
	print "Final: Norm = ", norm, "; Corr = ", corr
	

	return prop


def CalculateAngularMomentumDistribution(prop):
	dr = prop.psi.GetRepresentation().GetRadialRepresentation().GetRange(0).Dx
	lmax = prop.Propagator.LmRepresentation.Range.MaxL

	lmDistrib = sum(abs(prop.psi.GetData())**2, 0) * dr
	lDistrib = zeros(lmax+1, dtype=double)
	for l in r_[0:lmax+1]:
		for m in r_[-l:l+1]:
			lDistrib[l] += lmDistrib[MapLmIndex(l,m)]

	return lDistrib
	
	
def PlotAngularMomentumDistribution(prop):
	lDistrib = CalculateAngularMomentumDistribution(prop)
	bar(r_[0:len(lDistrib)], lDistrib)

	xlabel("Angular Momentum (l)")
	ylabel("Probability")
	title("Angular momentum distribution")
	return lDistrib


