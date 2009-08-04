from __future__ import with_statement
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

	if "pulseDuration" in args:
		T = args["pulseDuration"]
		conf.SetValue("PulseParameters", "pulse_duration", T)

	if "dt" in args:
		conf.SetValue("Propagation", "timestep", args["dt"])

	if "eigenvalueCount" in args:
		conf.SetValue("Arpack", "krylov_eigenvalue_count", args["eigenvalueCount"])

	if "eigenvalueBasisSize" in args:
		conf.SetValue("Arpack", "krylov_basis_size", args["eigenvalueBasisSize"])
	
	if "amplitude" in args:
		conf.SetValue("PulseParameters", "amplitude", args["amplitude"])
	
	if "frequency" in args:
		conf.SetValue("PulseParameters", "frequency", args["frequency"])
  
	if "phase" in args:
		phase = args['phase']
		if phase == "zero":
			conf.SetValue("LaserPotentialVelocityBase", "phase", 0.0)
		elif phase == "pihalf":
			conf.SetValue("LaserPotentialVelocityBase", "phase", pi/2.0)
		elif phase == "pi":
			conf.SetValue("LaserPotentialVelocityBase", "phase", pi)
		else:
			PrintOut("Unknown phase, using the one specified in the config file!")

	potentials = conf.Propagation.grid_potential_list + args.get("additionalPotentials", [])
	conf.SetValue("Propagation", "grid_potential_list", potentials)

	#Update config object from possible changed ConfigParser object
	newConf = pyprop.Config(conf.cfgObj)

	return newConf


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
#                       Propagation functions
#------------------------------------------------------------------------------------

def FormatDuration(duration):
	duration = int(duration)
	h = (duration / 3600)
	m = (duration / 60) % 60
	s = (duration % 60)

	str = []
	if h>0: str += ["%ih" % h]
	if m>0: str += ["%im" % m]
	if s>0: str += ["%is" % s]

	return " ".join(str)


def Propagate(**args):
	#Set up propagation problem
	potList = ["LaserPotentialLength"]
	prop = SetupProblem(additionalPotentials=potList, **args)

	#Setup initial state
	gsFilename = args.get("groundstateFilename", "helium_groundstate.h5")
	gsDatasetPath = args.get("groundstateDatasetPath", "/wavefunction")
	prop.psi.Clear()
	pyprop.serialization.LoadWavefunctionHDF(gsFilename, gsDatasetPath, prop.psi)
	initPsi = prop.psi.Copy()

	timeList = []
	corrList = []
	normList = []
	outsideAbsorberList = []

	tempPsi = prop.psi.Copy()

	#Setup absorber tests
	#absorberStart = prop.Config.Absorber.absorber_start
	#absorberEnd = absorberStart + prop.Config.Absorber.absorber_length
	#outsideAbsorberBox = SetupRadialMaskPotential(prop, absorberEnd, 100)
	#insideAbsorberBox = SetupRadialMaskPotential(prop, 0, absorberStart)

	#Propagate
	PrintOut("Starting propagation")
	outputCount = args.get("outputCount", 100)
	startTime = time.time()

	#Propagate until end of pulse
	duration = prop.Duration
	prop.Duration = duration
	for t in prop.Advance(outputCount, yieldEnd=True):
		#calculate values
		norm = prop.psi.GetNorm()
		corr = abs(initPsi.InnerProduct(prop.psi))**2

		#Calculate amount of stuff outside absorber
		#outsideAbsorber = abs(outsideAbsorberBox.GetExpectationValue(prop.psi, prop.GetTempPsi(), 0, 0))

		timeList.append(t)
		corrList.append(corr)
		normList.append(norm)

		#estimate remaining time
		curTime = time.time() - startTime
		totalTime = (curTime / t) * prop.Duration
		eta = totalTime - curTime
	
		#Print stats
		#PrintOut("t = %.2f; N = %.10f; Corr = %.10f, bound = %.10f, outside = %.10f, ETA = %s" % (t, norm, corr, boundTotal, outsideAbsorber, FormatDuration(eta)))
		PrintOut("t = %.2f; N = %.10f; Corr = %.10f, ETA = %s" % (t, norm, corr, FormatDuration(eta)))
		#PrintOut(prop.Propagator.Solver.GetErrorEstimateList())

	#Save the time-valued variables
	prop.TimeList = timeList
	prop.CorrList = corrList
	prop.NormList = normList
	#prop.OutsideAbsorberList = outsideAbsorberList

	PrintOut("")
	prop.Propagator.PampWrapper.PrintStatistics()

	#Saving final wavefunction
	#outputFilename = args.get("outputFilename", "final.h5")
	#outputDatasetPath = args.get("outputDatasetPath", "/wavefunction")
	#prop.SaveWavefunctionHDF(outputFilename, outputDatasetPath)

	return prop

