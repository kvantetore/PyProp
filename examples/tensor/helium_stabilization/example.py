import sys
import os
import time

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
	jscript.jobname = args.get(jobname, "pyprop")
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


