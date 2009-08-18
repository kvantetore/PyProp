from __future__ import with_statement
import sys
import os

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

from pyprop.serialization import RemoveExistingDataset

execfile("analysis.py")
execfile("analysis_single.py")
execfile("analysis_double.py")
execfile("analysis_eigenvalues.py")
execfile("eigenvalues.py")
execfile("../helium_stabilization/constants.py")

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

	if "grid" in args:
		radialGrid = args["grid"]
		def setvalue(variable):
			conf.SetValue("BsplineRepresentation", variable, radialGrid[variable])
		
		gridType = radialGrid["bpstype"]
		setvalue("bpstype")
		setvalue("order")
		setvalue("xmin")
		setvalue("xmax")
		setvalue("xsize")

		#update absorber to be at the end of the grid
		conf.Absorber.absorber_start = radialGrid["xmax"] - conf.Absorber.absorber_length

		if gridType == "linear":
			pass
		elif gridType == "exponentiallinear":
			setvalue("xpartition")
			setvalue("gamma")
		elif gridType == "centerexponentiallinear":
			setvalue("xpartition")
			setvalue("gamma")
		else:
			raise Exception("Invalid Grid Type '%s'" % gridType)	

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
	useArpack = args.get("useArpack", False)
	prop = SetupProblem(**args)
	if useArpack:
		solver = pyprop.ArpackSolver(prop)
	else:
		solver = pyprop.PiramSolver(prop)
	solver.Solve()
	print solver.GetEigenvalues()
	return solver


def GetGridPostfix(**args):
	"""
	Returns "unique" list of strings string identifying the grid
	implied by the specified args
	"""
	conf = SetupConfig(**args)
	cfg = conf.BsplineRepresentation

	gridType = cfg.bpstype
	postfix = ["grid", gridType, "xmax%i" % cfg.xmax, "xsize%i" % cfg.xsize, "order%i" % cfg.order]
	if gridType == "linear":
		pass
	elif gridType == "exponentiallinear" or gridType == "centerexponentiallinear":
		postfix.append("xpartition%i" % cfg.xpartition)
		postfix.append("gamma%.1f" % cfg.gamma)

	return postfix


def GetPulsePostfix(**args):
	"""
	Returns "unique" list of strings string identifying the pulse parameters
	implied by the specified args
	"""
	conf = SetupConfig(**args)
	cfg = conf.PulseParameters
	
	freq = round(cfg.frequency, 2)
	amplitude = round(cfg.amplitude, 2)
	cycles = cfg.cycles
	postfix = ["frequency%s" % freq, "amplitude%s" % amplitude, "cycles%i" % cycles]

	return postfix


def GroundstateFilename(**args):
	gridPostfix = GetGridPostfix(**args)
	filenamebase = "eigenstates/groundstate" + "_".join(gridPostfix) + ".h5"
	return filenamebase


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


def Propagate(saveWavefunction = False, **args):
	#Set up propagation problem
	potList = ["LaserPotentialLength"]
	prop = SetupProblem(additionalPotentials=potList, **args)

	#Setup initial state
	gsFilename = args.get("groundstateFilename", GroundstateFilename(**args))
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

	#Generate name of output file and path
	krylovSize = prop.Config.Propagation.krylov_basis_size
	T = round(prop.Duration, 1)
	dt = round(prop.TimeStep.real, 3)
	gridpostfix = GetGridPostfix(config=prop.Config)
	pulsepostfix = GetPulsePostfix(config=prop.Config)
	outputFilename = "output/convergence/propagation_krylov_%s_%s_%s_T_%s_dt_%1.2f" % (krylovSize, "_".join(gridpostfix), "_".join(pulsepostfix), T, dt)
	#outputFilename = args.get("outputFilename", outputFilename)
	outputDatasetPath = args.get("outputDatasetPath", "/wavefunction")

	#Propagate
	PrintOut("Starting propagation")
	outputCount = args.get("outputCount", 100)
	startTime = time.time()

	#Propagate until end of pulse
	duration = prop.Duration
	prop.Duration = duration
	for stepNum, t in enumerate(prop.Advance(outputCount, yieldEnd=True)):
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
		PrintOut("t = %.2f; N = %.15f; Corr = %.15f, ETA = %s" % (t, norm, corr, FormatDuration(eta)))
		#PrintOut(prop.Propagator.Solver.GetErrorEstimateList())

		if saveWavefunction:
			fname = outputFilename + "_%03i.h5" % stepNum
			prop.SaveWavefunctionHDF(fname, outputDatasetPath)

	#Save the time-valued variables
	prop.TimeList = timeList
	prop.CorrList = corrList
	prop.NormList = normList
	#prop.OutsideAbsorberList = outsideAbsorberList

	PrintOut("")
	prop.Propagator.PampWrapper.PrintStatistics()

	#Save wavefunction
	finalOutputFilename = outputFilename + ".h5"
	prop.SaveWavefunctionHDF(finalOutputFilename, outputDatasetPath)

	#Save propagation data
	h5file = tables.openFile(finalOutputFilename, "r+")
	try:
		h5file.createArray("/", "norm", prop.NormList)
		h5file.createArray("/", "corr", prop.CorrList)
		h5file.createArray("/", "times", prop.TimeList)
	finally:
		h5file.close()

	return prop


#------------------------------------------------------------------------------------
#                                     Misc I/O
#------------------------------------------------------------------------------------
def GetGroundstate(**args):
	groundstateFilename = GroundstateFilename(**args)
	h5file = tables.openFile(groundstateFilename, "r")
	try:
		groundstate = h5file.root.wavefunction[:]
	finally:
		h5file.close()
	
	return groundstate


def SaveGroundstate(**args):
	"""
	Find ground state using Piram, and save it to a HDF5 file with
	name generated by 'GroundstateFilename'
	"""
	#Find some bound states
	solver = FindEigenvalues(**args)

	#Put the groundstate into a wavefunction
	solver.SetEigenvector(solver.BaseProblem.psi, 0)

	#Store groundstate
	fname = GroundstateFilename(**args)
	solver.BaseProblem.SaveWavefunctionHDF(fname, "/wavefunction")

