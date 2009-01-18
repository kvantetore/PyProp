import sys

sys.path.append("./pyprop")
import pyprop
pyprop = reload(pyprop)
pyprop.ProjectNamespace = globals()

from pyprop import PrintOut

import numpy
import pylab
import time

from numpy import *
from libpotential import *

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
			conf.Propagation.timestep = -1.0j * abs(conf.Propagation.timestep)
			conf.Propagation.renormalization = True
		else:
			conf.Propagation.timestep = abs(conf.Propagation.timestep)
			conf.Propagation.renormalization = False

	if "duration" in args:
		duration = args["duration"]
		conf.Propagation.duration = duration

	if "eigenvalueCount" in args:
		conf.Arpack.krylov_eigenvalue_count = args["eigenvalueCount"]

	if "amplitude" in args:
		amp = args["amplitude"]
		freq = conf.LaserPotentialVelocity1.frequency
		conf.LaserPotentialVelocity1.amplitude = amp / freq
		conf.LaserPotentialVelocity2.amplitude = amp / freq
		conf.LaserPotentialVelocity3.amplitude = amp / freq

	additionalPotentials = args.get("additionalPotentials", [])
	conf.Propagation.grid_potential_list += additionalPotentials

	return conf


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

	def Setup(self, prop, dt):
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
		scalingS = 1.0
		scalingH = (1.0j*dt/2.)
		self.Solver.Setup(prop.psi, self.PreconditionRank, scalingS, scalingH)
		print "Setting up!"

	def Solve(self, psi):
		self.Solver.Solve(psi)

