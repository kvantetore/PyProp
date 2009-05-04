import sys

sys.path.append("./pyprop")
import pyprop
pyprop = reload(pyprop)
pyprop.ProjectNamespace = globals()

from pyprop import PrintOut

import numpy
import pylab
import time
import tables
import scipy.interpolate
import scipy.special

from pylab import *
from numpy import *
from libpotential import *

execfile("stabilization.py")
execfile("eigenvalues.py")
execfile("analysis.py")

#------------------------------------------------------------------------------------
#                       Setup Functions
#------------------------------------------------------------------------------------


def SetupConfig(**args):
	configFile = args.get("config", "config.ini")
	conf = pyprop.Load(configFile)
	
	if "silent" in args:
		silent = args["silent"]
		conf.SetValue("Propagation", "silent", args["silent"])

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

	if "dt" in args:
		conf.SetValue("Propagation", "timestep", args["dt"])

	if "eigenvalueCount" in args:
		conf.SetValue("Arpack", "krylov_eigenvalue_count", args["eigenvalueCount"])

	if "eigenvalueBasisSize" in args:
		conf.SetValue("Arpack", "krylov_basis_size", args["eigenvalueBasisSize"])
	
	if "radialGrid" in args:
		radialGrid = args["radialGrid"]
		def setvalue(variable):
			conf.SetValue("RadialRepresentation", variable, radialGrid[variable])
		
		gridType = radialGrid["bpstype"]
		setvalue("bpstype")
		setvalue("order")
		setvalue("xmax")
		setvalue("xsize")

		#update absorber to be at the end of the grid
		conf.Absorber.absorber_start = radialGrid["xmax"] - conf.Absorber.absorber_length

		if gridType == "linear":
			pass
		elif gridType == "exponentiallinear":
			setvalue("xpartition")
			setvalue("gamma")
		else:
			raise Exception("Invalid Grid Type '%s'" % gridType)

	if "model" in args:
		model = args["model"].lower()
		if model == "he+":
			z = 2.0
			a1 = a2 = a3 = a4 = a5 = a6 = 0

		elif model == "h":
			z = 1.0
			a1 = a2 = a3 = a4 = a5 = a6 = 0

		elif model == "he":
			z = 1.0
			a1 = 1.231
			a2 = 0.662
			a3 = -1.325
			a4 = 1.236 
			a5 = -0.231
			a6 = 0.480

		elif model == "he2":
			z = sqrt(2 * 1.45)
			a1 = a2 = a3 = a4 = a5 = a6 = 0

		else:
			raise Exception("Unknown model '%s'" % model)

		conf.SetValue("SAEPotential", "z", z)
		conf.SetValue("SAEPotential", "a1", a1)
		conf.SetValue("SAEPotential", "a2", a2)
		conf.SetValue("SAEPotential", "a3", a3)
		conf.SetValue("SAEPotential", "a4", a4)
		conf.SetValue("SAEPotential", "a5", a5)
		conf.SetValue("SAEPotential", "a6", a6)


	if "eigenvalueShift" in args:
		conf.SetValue("GMRES", "shift", args["eigenvalueShift"])

	if "amplitude" in args:
		conf.SetValue("PulseParameters", "amplitude", args["amplitude"])
	
	if "frequency" in args:
		conf.SetValue("PulseParameters", "frequency", args["frequency"])

	potentials = conf.Propagation.grid_potential_list + args.get("additionalPotentials", [])
	conf.SetValue("Propagation", "grid_potential_list", potentials)

	#Update config object from possible changed ConfigParser object
	newConf = pyprop.Config(conf.cfgObj)

	return newConf


def GetRadialGridPostfix(**args):
	"""
	Returns "unique" list of strings string identifying the radial grid
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

	return postfix

def GetModelPostfix(**args):
	return ["model_%s" % args.get("model", "he")]


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


def LaserFunctionFlatTop(conf, t):
	curField = 0
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
#------------------------------------------------------------------------------------
#                       Preconditioner for Cayley Propagator
#------------------------------------------------------------------------------------

class BSplinePreconditioner:
	def __init__(self, psi):
		self.Rank = psi.GetRank()
		#self.Solver = pyprop.CreateInstanceRank("FiniteDifferenceSolver", self.Rank, globals=globals())
		self.Solver = pyprop.CreateInstanceRank("core.BSpline_BSplineSolver", self.Rank)
		self.psi = psi

	def ApplyConfigSection(self, conf):
		self.PreconditionRank = conf.rank
		self.PotentialSections = [conf.Config.GetSection(s) for s in conf.potential_evaluation]

	def SetHamiltonianScaling(self, scalingH):
		self.HamiltonianScaling = scalingH

	def SetOverlapScaling(self, scalingS):
		self.OverlapScaling = scalingS

	def GetHamiltonianScaling(self):
		return self.HamiltonianScaling

	def GetOverlapScaling(self):
		return self.OverlapScaling

	def Setup(self, prop):
		#Add all potentials to solver
		for conf in self.PotentialSections:
			#Setup potential in basis
			potential = prop.BasePropagator.GeneratePotential(conf)
			for i, geom in enumerate(potential.GeometryList):
				if i==self.PreconditionRank and geom.GetStorageId() != "BandNH":
					raise "Potentials must be banded-nonhermitian in precondition rank"
				elif i!=self.PreconditionRank and geom.GetStorageId() != "Diag":
					raise "Potentials must be diagonal in non-precondition ranks"
		
			#Add potential to solver
			self.Solver.AddTensorPotential(potential.PotentialData)

		#Setup solver
		scalingS = self.GetOverlapScaling()
		scalingH = self.GetHamiltonianScaling() 
		self.Solver.Setup(prop.psi, self.PreconditionRank, scalingS, scalingH)

	def Solve(self, psi):
		self.Solver.Solve(psi)

