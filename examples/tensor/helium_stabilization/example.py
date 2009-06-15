#System Modules
from __future__ import with_statement
import sys
import os
import time
import commands
from datetime import timedelta

#Pyprop
sys.path.append("./pyprop")
import pyprop
pyprop = reload(pyprop)
pyprop.ProjectNamespace = globals()
from pyprop import PrintOut
from pyprop.serialization import RemoveExistingDataset

#Import numpy and such
import numpy
import pylab
from numpy import *
from pylab import *
from libpotential import *
import tables

#Try to import scipy and pysparse
try:
	import pysparse
	import scipy.linalg
except:
	pyprop.PrintOut("Could not load pysparse/scipy.linalg.")

try:
	import scipy
	import scipy.interpolate
	import scipy.special
	import scipy.sparse
	import scipy.linsolve
	scipy.linsolve.use_solver(useUmfpack=False)
except:
	pyprop.PrintOut("Could not load scipy")

#Setup INSTALLATION variable to allow us to identify the
#machine we're currently running on
INSTALLATION = os.environ.get("INSTALLATION", "local")
if INSTALLATION == "hexagon":
	import pyprop.utilities.submitpbs_hexagon as submitpbs
elif INSTALLATION == "stallo":
	import pyprop.utilities.submitpbs_stallo as submitpbs

execfile("stabilization.py")
execfile("twoelectron_test.py")
execfile("benchmark.py")
execfile("eigenvalues.py")
execfile("analysis.py")
execfile("bigmatrix.py")
execfile("preconditioner.py")
execfile("spectrum_finder.py")
execfile("analysis_eigenstate.py")
execfile("submit.py")
execfile("serialization.py")
execfile("analysis_single.py")
execfile("analysis_double.py")
execfile("tpdi.py")
execfile("constants.py")

#------------------------------------------------------------------------------------
#                       Setup Functions
#------------------------------------------------------------------------------------

def SetupConfig(**args):
	configFile = args.get("config", "config.ini")
	#if configfile is a string, load it, otherwise, treat it as
	#a config parser object
	if isinstance(configFile, str):
		conf = pyprop.Load(configFile)
	elif isinstance(configFile, pyprop.Config):
		conf = configFile
	else:
		conf = pyprop.Config(configFile)
	
	if "silent" in args:
		silent = args["silent"]
		conf.Propagation.silent = silent

	if "lmax" in args or "L" in args:
		lmax = args["lmax"]
		L = args["L"]
		indexIterator = pyprop.DefaultCoupledIndexIterator(lmax=lmax, L=L)
		conf.SetValue("AngularRepresentation", "index_iterator", indexIterator)

	if "imtime" in args:
		imtime = args["imtime"]
		if imtime:
			conf.SetValue("Propagation", "timestep", 1.0j * abs(conf.Propagation.timestep))
			conf.SetValue("Propagation", "renormalization", True)
		else:
			conf.SetValue("Propagation", "timestep", abs(conf.Propagation.timestep))
			conf.SetValue("Propagation", "renormalization", False)

	if "duration" in args:
		duration = args["duration"]
		conf.SetValue("Propagation", "duration", duration)
	
	if "pulse_duration" in args:
		pulse_duration = args["pulse_duration"]
		conf.SetValue("PulseParameters", "pulse_duration", pulse_duration)
		conf.SetValue("PulseParameters", "cycles", None)
	
	if "timestep" in args:
		dt = args["timestep"]
		conf.SetValue("Propagation", "timestep", dt)

	if "eigenvalueCount" in args:
		conf.SetValue("Arpack", "krylov_eigenvalue_count", args["eigenvalueCount"])

	if "eigenvalueBasisSize" in args:
		conf.SetValue("Arpack", "krylov_basis_size", args["eigenvalueBasisSize"])

	if "eigenvalueShift" in args:
		conf.SetValue("GMRES", "shift", args["eigenvalueShift"])
		print "Using shift ", args["eigenvalueShift"]

	if "shift" in args:
		conf.SetValue("GMRES", "shift", args["shift"])
		print "Using shift ", args["shift"]

	if "index_iterator" in args:
		conf.SetValue("AngularRepresentation", "index_iterator", args["index_iterator"])

	if "amplitude" in args:
		conf.SetValue("PulseParameters", "amplitude", args["amplitude"])

	if "frequency" in args:
		conf.SetValue("PulseParameters", "frequency", args["frequency"])
		#PrintOut("    Setting new field frequency: %s" % args["frequency"])

	if "multipoleCutoff" in args:
		conf.SetValue("ElectronicCouplingPotential", "geometry0", "SelectionRule_R12_%i" % args["multipoleCutoff"])

	if "radialGrid" in args:
		radialGrid = args["radialGrid"]
		def setvalue(variable):
			conf.SetValue("RadialRepresentation", variable, radialGrid[variable])
		
		gridType = radialGrid["bpstype"]
		setvalue("bpstype")
		setvalue("order")
		setvalue("xmax")
		setvalue("xsize")

		if gridType == "linear":
			pass
		elif gridType == "exponentiallinear":
			setvalue("xpartition")
			setvalue("gamma")
		elif gridType == "exponential":
			setvalue("gamma")
		else:
			raise Exception("Invalid Grid Type '%s'" % gridType)

		conf.SetValue("Absorber", "absorber_start", float(radialGrid["xmax"]) - conf.Absorber.absorber_length)

	if args.get("useDefaultPotentials", True):
		potentials = conf.Propagation.grid_potential_list 
	else:
		potentials = []
	potentials += args.get("additionalPotentials", [])
	if "xmax" in args:
		conf.SetValue("RadialRepresentation", "xmax", args["xmax"])
	
	if "xsize" in args:
		conf.SetValue("RadialRepresentation", "xsize", args["xsize"])
	
	if "order" in args:
		conf.SetValue("RadialRepresentation", "order", args["order"])

	if "phase" in args:
		phase = args['phase']
		if phase == "zero":
			conf.SetValue("LaserPotentialVelocityBase", "phase", 0.0)
		elif phase == "pihalf":
			conf.SetValue("LaserPotentialVelocityBase", "phase", pi/2.0)
		else:
			PrintOut("Unknown phase, using the one specified in the config file!")

	conf.SetValue("Propagation", "grid_potential_list", potentials)

	if args.get("useDefaultPreconditionPotentials", True):
		precondPotentials = conf.RadialPreconditioner.potential_evaluation 
	else:
		precondPotentials = []
	precondPotentials += args.get("additionalPreconditionPotentials", [])
	conf.SetValue("RadialPreconditioner", "potential_evaluation", precondPotentials)

	useStoredPotentials = args.get("useStoredPotentials", False)
	if useStoredPotentials:
		postfix = "_".join(GetRadialGridPostfix(config=conf) + GetAngularGridPostfix(config=conf))
		folder = "output/potentials/%s" % (postfix)

		potentialNames =potentialNames = [\
			"RadialKineticEnergy1", \
			"RadialKineticEnergy2", \
			"AngularKineticEnergy", \
			"CoulombPotential", \
			"ElectronicCouplingPotential", \
			"ElectronicCouplingPotentialMonopoleTerm", \
			"OverlapPotential", \
			"Absorber", \
			"LaserPotentialVelocityDerivativeR1", \
			"LaserPotentialVelocityDerivativeR2", \
			"LaserPotentialVelocity", \
			"DipolePotentialLength", \
			"SingleIonizationBox", \
			"DoubleIonizationBox", \
			]
		for potName in potentialNames:
			conf.SetValue(potName, "filename", os.path.join(folder, "%s.h5" % potName))
			conf.SetValue(potName, "dataset", "potential")
 

	#Update config object from possible changed ConfigParser object
	newConf = pyprop.Config(conf.cfgObj)

	return newConf


def GetRadialGrid(**args):
	"""
	returns the radial grid info as a dict
	"""
	conf = SetupConfig(**args)
	cfg = conf.RadialRepresentation

	#Get some of the RadialRepresentation values into a dict
	names = ["bpstype", "xmax", "xsize", "order", "xpartition", "gamma"]
	radialGrid = dict(map(lambda name: (name, getattr(cfg, name)), names))


	conf.SetValue("Propagation", "grid_potential_list", potentials)

	#Update config object from possible changed ConfigParser object
	newConf = pyprop.Config(conf.cfgObj)

	return newConf


def GetRadialGrid(**args):
	"""
	returns the radial grid info as a dict
	"""
	conf = SetupConfig(**args)
	cfg = conf.RadialRepresentation

	#Get some of the RadialRepresentation values into a dict
	names = ["bpstype", "xmax", "xsize", "order", "xpartition", "gamma"]
	radialGrid = dict(map(lambda name: (name, getattr(cfg, name)), names))

	return radialGrid


def GetRadialGridPostfix(**args):
	"""
	Returns a "unique" list of strings string identifying the radial grid
	implied by the specified args
	"""
	conf = SetupConfig(**args)
	cfg = conf.RadialRepresentation

	gridType = cfg.bpstype
	postfix = ["grid", gridType, "xmax%i" % cfg.xmax, "xsize%i" % cfg.xsize, "order%i" % cfg.order]
	if gridType == "linear":
		pass
	elif gridType == "exponentiallinear":
		postfix.append("xpartition%i" % cfg.xpartition)
		postfix.append("gamma%.1f" % cfg.gamma)
	elif gridType == "exponential":
		postfix.append("gamma%.1f" % cfg.gamma)

	return postfix


def GetAngularBasisSize(**args):
	"""
	Returns the number of spherical harmonic basis functions for a 
	given argument list
	"""
	conf = SetupConfig(**args)
	return len([1 for i in conf.AngularRepresentation.index_iterator])


def GetAngularGridPostfix(**args):
	"""
	Returns a "unique" list of strings string identifying the angular grid
	implied by the specified args
	"""
	conf = SetupConfig(**args)
	cfg = conf.AngularRepresentation

	lmax = max([l1 for l1, l2, L, M in cfg.index_iterator])
	Llist = unique([L for l1, l2, L, M in cfg.index_iterator])
	Mlist = unique([M for l1, l2, L, M in cfg.index_iterator])

	def getSortedListString(l):
		if len(l) == 1:
			string = "%i" % l[0]
		else:
			if (diff(l) == 1).all():
				string = "%i-%i" % (l[0], l[-1]+1)
			else:
				string = "%s" % ("_".join(map(str, l)))

		return string

	postfix = ["angular"]
	postfix += ["lmax%i" % lmax]
	postfix += ["L%s" % getSortedListString(Llist)]	
	postfix += ["M%s" % getSortedListString(Mlist)]	

	return postfix


def CheckCompatibleRadialGrid(conf1, conf2):
	return GetRadialGridPostfix(config=conf1) == GetRadialGridPostfix(config=conf2)


def CheckCompatibleAngularGrid(conf1, conf2, check_L=True, check_M=True):
	postfix1 = GetAngularGridPostfix(config=conf1)
	postfix2 = GetAngularGridPostfix(config=conf2)

	if postfix1[1] != postfix2[1]:
		return False
	if check_L and postfix1[2] != postfix2[2]:
		return False
	if check_M and postfix1[3] != postfix2[3]:
		return False
	return True

	
def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	return prop


#------------------------------------------------------------------------------------
#                       Basic Propagation Functions
#------------------------------------------------------------------------------------


def FindGroundstate(**args):
	prop = SetupProblem(imtime=True, **args)
	for t in prop.Advance(True):
		E = prop.GetEnergy()
		if pyprop.ProcId == 0:
			print "t = %02.2f, E = %2.8f" % (t, E)

	E = prop.GetEnergyExpectationValue()
	print "t = %02.2f, E = %2.8f" % (t, E)

	return prop


def FindEigenvalues(useArpack = False, **args):
	prop = SetupProblem(**args)
	if useArpack:
		solver = pyprop.ArpackSolver(prop)
	else:
		solver = pyprop.PiramSolver(prop)
	solver.Solve()
	print solver.GetEigenvalues()
	return solver



#------------------------------------------------------------------------------------
#                       Laser Time Functions
#------------------------------------------------------------------------------------

def LaserFunctionVelocity(conf, t):
	phase = 0.0
	if hasattr(conf, "phase"):
		phase = conf.phase
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= - cos(t * conf.frequency + phase);
	else:
		curField = 0
	return curField

def LaserFunctionFlatTop(conf, t):
	pulseStart = 0
	if conf.Exists("pulse_start"):
		pulseStart = conf.pulse_start
	curField = conf.amplitude * cos(t * conf.frequency);

	if (t > conf.pulse_duration) or (t < pulseStart):
		curField = 0
	elif 0 <= t < conf.ramp_on_time:
		curField *= sin(t * pi / (2*conf.ramp_on_time))**2;
	elif t > conf.pulse_duration - conf.ramp_off_time:
		curField *= sin((conf.pulse_duration - t) * pi / (2*conf.ramp_off_time))**2;
	else:
		curField *= 1

	return curField

def LaserFunctionSimpleLength(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= cos(t * conf.frequency);
	else:
		curField = 0
	return curField


def LaserFunctionLength(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		T = conf.pulse_duration
		w = conf.frequency
		curField *= sin(pi*t/T)*(-2*pi/T * cos(pi*t/T) * cos(w*t) + w*sin(pi*t/T)*sin(w*t))
	else:
		curField = 0
	return curField


#------------------------------------------------------------------------------------
#                       Debug Functions
#------------------------------------------------------------------------------------

def GetBasisPairs(selectionRule, indexIterator):
	class reprConfigSection(pyprop.Section):
		def __init__(self):
			self.index_iterator = indexIterator

	cfg = reprConfigSection()
	repr = pyprop.core.CoupledSphericalHarmonicRepresentation()
	cfg.Apply(repr)

	return selectionRule.GetBasisPairs(repr)


#------------------------------------------------------------------------------------
#                      Misc Functions
#------------------------------------------------------------------------------------
def CalculatePolarizabilityGroundState(fieldRange = linspace(0.001,0.01,10), **args):
	"""
	Calculate polarizability of ground state using the Jacobi-Davidson method. We use
	the formula

	    polarizability = -2 * E_shift / field**2,
	
	thus assuming that the hyperpolarizability term is negligible.
	"""
	def JacobiDavidson(TotalMatrix):
		numConv, E, V, numIter, numIterInner = pysparse.jdsym.jdsym( \
			TotalMatrix, OverlapMatrix, None, 1, -3.0, 1e-10, 100, \
			pysparse.itsolvers.qmrs) 

		return E

	def SetupTotalMatrix(fieldStrength):
		#Add H matrix
		idx = H.keys()
		for key in idx:
			I = key[0]
			J = key[1]
			TotalMatrixLL[I,J] = H[I,J]

		#Add field matrix
		idx2 = FieldMatrix.keys()
		for key in idx2:
			I = key[0]
			J = key[1]
			TotalMatrixLL[I,J] = fieldStrength * FieldMatrix[I,J]

		TotalMatrix = TotalMatrixLL.to_sss()
		return TotalMatrix

	#Setup problem
	prop = SetupProblem(**args)

	#Set up the field-free hamilton matrix
	H = SetupPotentialMatrixLL(prop, [0,1])

	#Set up the field matrix
	FieldMatrix = SetupPotentialMatrixLL(prop, [2])

	#Set up overlap matrix
	S_ll = SetupPotentialMatrixLL(prop,[3])
	OverlapMatrix = S_ll.to_sss()

	#Setup total system matrix
	TotalMatrixLL = H.copy()

	#Get reference energy
	TotalMatrix = SetupTotalMatrix(0.0)
	refEnergy = JacobiDavidson(TotalMatrix)
	print "Reference energy = %.15f" % refEnergy

	#Store shifted energies and polarizabilities
	energies = []
	polarizabilities = []
	
	for fieldStrength in fieldRange:
		print "Calculating energy shift for field strength = %s..." % fieldStrength
		#Reset potential with current intensity
		TotalMatrix = SetupTotalMatrix(fieldStrength)
		
		#Find shifted eigenvalue
		energy = JacobiDavidson(TotalMatrix)

		#Compute polarizability
		energyShift = energy - refEnergy
		polarizability = -2.0 * energyShift / fieldStrength**2

		#Store current energy and polarizability
		energies.append(energy)
		polarizabilities.append(polarizability)

	return fieldRange, polarizabilities, energies

#------------------------------------------------------------------------------------
#                      Wavefunction Analysis Function
#------------------------------------------------------------------------------------
def FindIonizationProbability(datafile, boundstateFiles, ionizationThreshhold=-2.0):
	"""
	Find total single and double ionization of Helium by projecting on states with 
	energy < 2.0 a.u.
	"""

	conf = pyprop.Config(pyprop.serialization.GetConfigFromHDF5(datafile))
	lmax = conf.AngularRepresentation.index_iterator.lmax
	Lmax = conf.AngularRepresentation.index_iterator.L[-1]

	conf.Propagation.grid_potential_list = []
	conf.Propagation.preconditioner = None

	#h5file = tables.openFile(datafile)
	#try:
	#	ionizationProbability = h5file.root.Norm[0]
	#finally:
	#	h5file.close()
	ionizationProbability = 1.0
		
	#Set up problem
	#conf.AngularRepresentation.index_iterator = pyprop.DefaultCoupledIndexIterator(lmax=lmax, L=L)
	prop = pyprop.Problem(conf)
	tmpPsi = prop.psi.Copy()
	totalIdxIterator = pyprop.DefaultCoupledIndexIterator(lmax=lmax, L=range(Lmax))

	#Load wavefunction
	h5file = tables.openFile(datafile, "r")
	try:
		prop.psi.GetData()[:] = h5file.root.wavefunction[:]
	finally:
		h5file.close()
	for L in range(Lmax + 1):
		#Project on all bound states for current L
		print "    L = %i" % L
		h5file = tables.openFile(boundstateFiles.pop(0), "r")
		numEigs = size(h5file.root.Eig.Eigenvalues)
		for i in range(numEigs):
			tmpPsi.Clear()
			for j,cur in enumerate(totalIdxIterator):
				if cur.L == L and h5file.root.Eig.Eigenvalues[i] < ionizationThreshhold:
					datasetPath = GetEigenvectorDatasetPath(i)
					tmpPsi.GetData()[j,:,:] += array(h5file.getNode(datasetPath))[cur.l1, :, :]
			ionizationProbability -= abs(prop.psi.InnerProduct(tmpPsi))**2

		h5file.close()

	return ionizationProbability


def ReconstructRadialDensityOnGrid(datafileName, radialGrid, bsplineObject, angularIndexIterator):
	"""
	Construct grid radial density for a two-electron wavefunction.
	"""

	#Total number of B-splines
	numBsplines = bsplineObject.NumberOfBSplines

	#Array for final radial density
	radialDensity = zeros((radialGrid.size, radialGrid.size), dtype=double)

	#A buffer for reconstructing 
	buffer1D = zeros(radialGrid.size, dtype=complex)
	buffer2D = zeros((numBsplines, radialGrid.size), dtype=complex)
	psiSliceGrid = zeros(radialGrid.size, dtype=complex)
	psiSlice = zeros(numBsplines, dtype=complex)

	#Calculate radial density
	with tables.openFile(datafileName) as f:
		for lIdx, st in enumerate(angularIndexIterator):
			for i in range(numBsplines):
				psiSlice[:] = f.root.wavefunction[lIdx, i, :]
				buffer1D[:] = 0.0
				bsplineObject.ConstructFunctionFromBSplineExpansion(psiSlice, radialGrid, buffer1D)
				buffer2D[i,:] = buffer1D[:]

			for j in range(radialGrid.size):
				psiSlice[:] = buffer2D[:,j]
				buffer1D[:] = 0.0
				bsplineObject.ConstructFunctionFromBSplineExpansion(psiSlice, radialGrid, buffer1D)

				#Incoherent sum over partial waves (but coherent over b-spline coefficients)
				radialDensity[:,j] += abs(buffer1D[:])**2

	#Perform smoothing by a simple running average
	#runAvgNum = args.get("runAvgNum", 1)
	#avgStartIdx = args.get("avgStartIdx", 0)
	#gridSize = radialGrid.size
	#A = array([sum(radialDensity[start:start+runAvgNum,:],axis=0)/runAvgNum for start in range(avgStartIdx, gridSize-runAvgNum)])
	#B = array([sum(A[:, start:start+runAvgNum],axis=1)/runAvgNum for start in range(avgStartIdx, gridSize-runAvgNum)])

	#return radialDensity, B
	return radialDensity

def AssertSingleProc():
	if not pyprop.IsSingleProc():
		raise Exception("Method is only supported in single processor mode")

def CreateRadialDensityFromFile(datafile, radialGrid):
	#Get config from first data file
	conf = pyprop.serialization.GetConfigFromHDF5(datafile)
	indexIt = conf.get("AngularRepresentation", "index_iterator")

	#Set up radial problem w/o potential
	conf.set("Propagation", "preconditioner", 'None')
	conf.set("InitialCondition", "type", 'None')
	conf.set("Representation", "rank", "1")
	conf.set("Representation", "representation0", "'RadialRepresentation'")
	conf.set("Representation", "type", "core.CombinedRepresentation_1")
	conf.set("GMRES", "preconditioner", "None")
	prop = SetupProblem(config = conf, useDefaultPotentials = False, silent = False)

	#Get bspline object
	bspline = prop.psi.GetRepresentation().GetRepresentation(0).GetBSplineObject()
	
	#Construct radial density
	radialDensity = ReconstructRadialDensityOnGrid(datafile, radialGrid, bspline, indexIt)

	return radialDensity


def CreateRadialDensitySeries(datafilesPath, radialGrid):
	#Get file names
	fileList = sort(["%s/%s" % (datafilesPath, file) for file in os.listdir(datafilesPath) if file.find("doubleduration") == -1 ])

	#Get config from first data file
	conf = pyprop.serialization.GetConfigFromHDF5(fileList[0])
	indexIt = conf.get("AngularRepresentation", "index_iterator")

	#Set up radial problem w/o potential
	conf.set("Propagation", "preconditioner", 'None')
	conf.set("InitialCondition", "type", 'None')
	conf.set("Representation", "rank", "1")
	conf.set("Representation", "representation0", "'RadialRepresentation'")
	conf.set("Representation", "type", "core.CombinedRepresentation_1")
	conf.set("GMRES", "preconditioner", "None")
	prop = SetupProblem(config = conf, useDefaultPotentials = False, silent = False)

	#Get bspline object
	bspline = prop.psi.GetRepresentation().GetRepresentation(0).GetBSplineObject()

	#Outfile path
	outfilePath = "out/radial_densities_rmax_%s_rsize_%s" % (radialGrid[-1], len(radialGrid))
	try:
		os.makedirs(outfilePath)
	except:
		print "Folder already exists!"
		pass

	#Create the frames
	for file in fileList:
		print "Now processing: %s" % file.split("/")[-1]

		#Get radial density
		radialDensity = ReconstructRadialDensityOnGrid(file, radialGrid, bspline, indexIt)

		#Store radial density
		h5file = tables.openFile("%s/%s" % (outfilePath, file.split("/")[-1]), "w")
		try:
			h5file.createArray("/", "RadialDensity", radialDensity)
			h5file.createArray("/", "RadialGrid", radialGrid)
		finally:
			h5file.close()


def RunningAverage2D(data2d, radialGrid, **args):
	"""
	Perform a two-dimensional running average
	"""

	#Number of points to average over
	runAvgNum = args.get("runAvgNum", 1)

	#Index at which to start averaging
	avgStartIdx = args.get("avgStartIdx", 0)

	gridSize = radialGrid.size

	#Set up window function
	windowWidth = args.get("windowWidth", 5.0)
	gridSlice = linspace(-runAvgNum / 2.0, runAvgNum / 2.0, runAvgNum)
	windowFunc = outer(ones(gridSize), exp(-gridSlice**2 / windowWidth))

	#Average over first dim
	A = array([sum(transpose(windowFunc) * data2d[start:start+runAvgNum,:],axis=0)/runAvgNum for start in range(avgStartIdx, gridSize-runAvgNum)])

	#Average over second dim
	B = array([sum(windowFunc[0:shape(A)[0],:] * A[:, start:start+runAvgNum],axis=1)/runAvgNum for start in range(avgStartIdx, gridSize-runAvgNum)])

	#Normalize
	B /= linalg.norm(B)

	return B
