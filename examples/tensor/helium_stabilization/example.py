import sys
import os
import time
try:
	import pysparse
except:
	pass

from datetime import timedelta

sys.path.append("./pyprop")
import pyprop
pyprop = reload(pyprop)
pyprop.ProjectNamespace = globals()

from numpy import *
#from pylab import *
from libpotential import *
from pyprop import PrintOut

execfile("stabilization.py")

#------------------------------------------------------------------------------------
#                       Setup Functions
#------------------------------------------------------------------------------------


def SetupConfig(**args):
	configFile = args.get("config", "config.ini")
	conf = pyprop.Load(configFile)
	
	if "silent" in args:
		silent = args["silent"]
		conf.Propagation.silent = silent

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

	if "eigenvalueCount" in args:
		conf.SetValue("Arpack", "krylov_eigenvalue_count", args["eigenvalueCount"])

	if "index_iterator" in args:
		conf.SetValue("AngularRepresentation", "index_iterator", args["index_iterator"])

	if "amplitude" in args:
		conf.SetValue("PulseParameters", "amplitude", args["amplitude"])

	potentials = conf.Propagation.grid_potential_list + args.get("additionalPotentials", [])
	conf.SetValue("Propagation", "grid_potential_list", potentials)

	#Update config object from possible changed ConfigParser object
	newConf = pyprop.Config(conf.cfgObj)

	return newConf


def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	return prop


def FindGroundstate(**args):
	prop = SetupProblem(imtime=True, **args)
	for t in prop.Advance(10):
		E = prop.GetEnergy()
		if pyprop.ProcId == 0:
			print "t = %02.2f, E = %2.8f" % (t, E)

	E = prop.GetEnergyExpectationValue()
	print "t = %02.2f, E = %2.8f" % (t, E)

	return prop


def FindEigenvalues(**args):
	prop = SetupProblem(**args)
	solver = pyprop.PiramSolver(prop)
	solver.Solve()
	print solver.GetEigenvalues()
	return solver


def SetupPotentialMatrix(prop, whichPotentials):
	print "Setting up potential matrix..."
	matrixSize = prop.psi.GetData().size
	
	#Allocate the potential matrix
	print "    Allocating potential matrix of size [%i, %i]  ~%.0f MB" \
		% (matrixSize, matrixSize, matrixSize**2 * 16 / 1024.**2)
	BigMatrix = zeros((matrixSize, matrixSize), dtype="complex")

	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)

		for idxL, idxR, i, j, k in TensorPotentialIndexMap(prop.psi.GetData().shape, potential):
			BigMatrix[idxL, idxR] += potential.PotentialData[i, j, k]

	return BigMatrix


def SetupPotentialMatrixLL(prop, whichPotentials, eps=1e-14):
	def SetElement():
		MatrixLL[idxL, idxR] += potential.PotentialData[i, j, k].real

	print "Setting up potential matrix..."
	matrixSize = prop.psi.GetData().size

	#Find potential size for largest potential
	potList = prop.Propagator.BasePropagator.PotentialList
	potentialSize = max([potList[w].PotentialData.size for w in whichPotentials])

	#Set up linked list matrix
	MatrixLL = pysparse.spmatrix.ll_mat_sym(matrixSize, potentialSize)

	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)

		count = 0
		outStr = ""
		for idxL, idxR, i, j, k in TensorPotentialIndexMap(prop.psi.GetData().shape, potential):
			#Skip upper triangle (we have a symmetric matrix)
			if idxL < idxR:
				continue

			#Skip current element if less than eps
			if abs(potential.PotentialData[i, j, k]) < eps:
				continue

			#Print progress info
			countSize = 1e4
			if mod(count, countSize) == 0:
				outStr = " " * 8
				outStr += "Progress: %i/%i" % (count/countSize, round(potentialSize/2./countSize))
				sys.stdout.write("\b"*len(outStr) + outStr)
				sys.stdout.flush()

			#Store element in linked-list matrix
			#MatrixLL[idxL, idxR] += potential.PotentialData[i, j, k].real
			SetElement()
			count += 1

		print


	return MatrixLL


def TensorPotentialIndexMap(psiShape, tensorPotential):
	"""
	Returns a generator for a map between indices in an m x m matrix and 
	"""
	basisPairs0 = tensorPotential.BasisPairs[0]
	basisPairs1 = tensorPotential.BasisPairs[1]
	basisPairs2 = tensorPotential.BasisPairs[2]
	
	Count0 = psiShape[0]
	Count1 = psiShape[1]
	Count2 = psiShape[2]

	for i, (x0,x0p) in enumerate(basisPairs0):
		for j, (x1,x1p) in enumerate(basisPairs1):
			for k, (x2,x2p) in enumerate(basisPairs2):
				indexLeft = x2 + (x1 * Count2) + (x0 * Count1 * Count2) 
				indexRight = x2p + (x1p * Count2) + (x0p * Count1 * Count2) 
				yield indexLeft, indexRight, i, j, k


#------------------------------------------------------------------------------------
#                       Laser Time Functions
#------------------------------------------------------------------------------------

def LaserFunctionVelocity(conf, t):
	if 0 <= t < conf.pulse_duration:
		curField = conf.amplitude;
		curField *= sin(t * pi / conf.pulse_duration)**2;
		curField *= - cos(t * conf.frequency);
	else:
		curField = 0
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
#                       Job Submit Functions
#------------------------------------------------------------------------------------
installation = os.environ.get("INSTALLATION", "local")
if installation == "hexagon":
	import pyprop.utilities.submitpbs_hexagon as submitpbs
if installation == "stallo":
	import pyprop.utilities.submitpbs_stallo as submitpbs


def Submit(executable=None, writeScript=False, installation="hexagon", **args):
	"""
	Set up job scripts and other necessary stuff to run ionization rate
	cycle scan experiment.
	"""

	#Create jobscript 
	jscript = submitpbs.SubmitScript()
	jscript.jobname = args.get("jobname", "pyprop")
	jscript.walltime = timedelta(hours=args.get("runHours",1), minutes=0, seconds=0)
	jscript.ppn = args.get("ppn", 4)
	jscript.proc_memory = args.get("proc_memory", "1000mb")
	jscript.nodes = args.get("nodes", "1")
	jscript.interconnect = args.get("interconnect", "ib")
	numProcs = args.get("numProcs", jscript.ppn*jscript.nodes)
	if installation == "stallo":
		jscript.workingdir = args.get("workingDir", "/home/nepstad/proj/argon/")
		jscript.executable = "mpirun -n %s " % (jscript.ppn*jscript.nodes)
		jscript.executable += "python %s" % executable
	elif installation == "hexagon":
		jscript.workingdir = args.get("workingDir", "/work/nepstad/dev/argon/")
		jscript.executable = "aprun -n %s -N %s " % (numProcs, jscript.ppn)
		jscript.executable += "./pyprop-exec %s" % executable
	jscript.parameters = commands.mkarg(repr(args))
	jscript.account = args.get("account", "fysisk")

	#Submit this job
	if writeScript:
		print "\n".join(jscript.CreateScript())
	else:
		jscript.Submit()


#------------------------------------------------------------------------------------
#                       Serialization Functions
#------------------------------------------------------------------------------------
def SetupProblemFromFile(file, nodeName=None):
	"""
	Set up problem object and load wavefunction from file.
	"""
	prop = None
	cfgObj = pyprop.serialization.GetConfigFromHDF5(file)
	cfgObj.set("InitialCondition", "type", "None")
	conf = pyprop.Config(cfgObj)
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	GetWavefunctionFromFile(file, prop.psi, nodeName=nodeName)
	
	return prop


def GetWavefunctionFromFile(file, psi, nodeName=None):
	h5file = tables.openFile(file, "r")
	try:
		if nodeName == None:
			for node in h5file.walkNodes():
				if node._v_name == "wavefunction":
					psi.GetData()[:] = node[:]
		else:
			psi.GetData()[:] = h5file.getNode(nodeName)[:]
	finally:
		h5file.close()


def GetArrayFromFile(file, nodeName):
	h5file = tables.openFile(file, "r")
	try:
		for node in h5file.walkNodes():
			if node._v_name == nodeName:
				dataArray = node[:]
	finally:
		h5file.close()
	
	return dataArray


def StoreTensorPotentialMTX(prop, whichPotentials, outFileName, eps = 1e-14):
	fh = open(outFileName, "w")
	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)

		basisPairs0 = potential.BasisPairs[0]
		basisPairs1 = potential.BasisPairs[1]
		basisPairs2 = potential.BasisPairs[2]
		
		Count0 = prop.psi.GetData().shape[0]
		Count1 = prop.psi.GetData().shape[1]
		Count2 = prop.psi.GetData().shape[2]

		for i, (x0,x0p) in enumerate(basisPairs0):
			for j, (x1,x1p) in enumerate(basisPairs1):
				for k, (x2,x2p) in enumerate(basisPairs2):
					indexLeft = x2 + (x1 * Count2) + (x0 * Count1 * Count2) 
					indexRight = x2p + (x1p * Count2) + (x0p * Count1 * Count2) 

					#Skip current element if less than eps
					if abs(potential.PotentialData[i, j, k]) < eps:
						continue

					#Write data line to file
					potReal = potential.PotentialData[i, j, k].real
					potImag = potential.PotentialData[i, j, k].imag
					outStr = "%i %i %1.16f %1.16f\n" % (indexLeft, indexRight, potReal, potImag)
					fh.write(outStr)

	fh.close()
