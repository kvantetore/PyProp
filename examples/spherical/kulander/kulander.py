#improt system modules
import sys
import os

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("../../.."))

#Load and reload pyprop in order to get recent changes
import pyprop
pyprop = reload(pyprop)
from libkulander import *

import submitpbs
from datetime import timedelta

#numpy an pylab for good measure
from pylab import *
from numpy import *

#exponential fit:
def ExponentialFit(samplePoints, sampleValues):
	x = samplePoints
	y = log(sampleValues)
	[c1, c0] = polyfit(x, y, 1)

	decayRate = - c1
	normFactor = exp(c0)
	return decayRate, normFactor

def PlotExpFit(t, n):
	decayRate, const = ExponentialFit(t, n)
	plot(t, const * exp(- decayRate * t))

#Choose radial grid type:
class RadialGridType:
	CARTESIAN = 1
	TRANSFORMED = 2

def SetRadialGridType(conf, gridType):
	#Set RadialRepresentation section
	if gridType == RadialGridType.CARTESIAN:
		conf.RadialRepresentation = conf.CartesianRadialRepresentation

	elif gridType == RadialGridType.TRANSFORMED:
		conf.RadialRepresentation = conf.TransformedRadialRepresentation

	else:
		raise "Invalid grid type ", gridType

	#Set radialtype 
	conf.Representation.radialtype = conf.RadialRepresentation.type

def test(gridType=RadialGridType.CARTESIAN):
	conf = pyprop.Load("kulander.ini")
	SetRadialGridType(conf, gridType)

	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop

def MapLmIndex(l, m):
	return (l + 1) * l + m

def RunKulanderExperiment(gridType=RadialGridType.CARTESIAN):
	#load config file. hydrogen.ini uses ../sphericalbase.ini as a base
	#configuration file, so be sure to check out that one as well.
	conf = pyprop.Load("kulander.ini")
	SetRadialGridType(conf, gridType)

	#Uses the same Problem class, only changes the propagator specified in the
	#config file
	prop = pyprop.Problem(conf)

	#Set up all transformations and potentials.
	prop.SetupStep()

	#Make sure the wavefunction is normalized
	prop.psi.Normalize()

	#save the initial wavefunction
	initPsi = prop.psi.Copy()

	#propagate through the problem, and do something 
	#every timestep
	
	index = 0
	for t in prop.Advance(500):
		#pyprop.Plot2DRank(prop, 0)
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		print "t = ", t, "; Norm = ", norm, "; Corr = ", corr
		
		ldist = CalculateAngularMomentumDistribution(prop)
		pyprop.AppendMatlabArray("output/ldist", ldist)

		index += 1

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

def SubmitJob():
	script = submitpbs.SubmitScript()
	script.jobname = "kulander"
	
	#resources
	script.ppn = 1
	script.nodes = 1
	script.walltime = timedelta(hours=10)
	script.proc_memory = "500mb"
		
	#stdio
	script.stdout = "output/kulander_stdout"
	script.stderr = "output/kulander_stderr"
		
	#mpi
	script.mpi = "openmpi"
	script.compiler = "intel"
	script.mpirun_enabled = False
		
	#exec
	script.executable = "./run.py"
	script.parameters = ""

	#go
	script.Submit()
#	for line in script.CreateScript():
#		print line 

def LoadOutput(outputFile):
	f = open(outputFile,'r')
	
	lines = f.readlines()
	del f
	
	tVec = zeros(len(lines), dtype=double)
	normVec = zeros(len(lines), dtype=double)
	corrVec = zeros(len(lines), dtype=double)
	
	i = 0
	for line in lines:
	    if line.startswith("t ="):
	        try:
	            for segment in line.split(","):
	                segment = segment.strip()
	                exec(segment, globals(), globals())
	            tVec[i] = t
	            normVec[i] = Norm
	            corrVec[i] = Corr
	            i = i+1
	        except:
				print "error"

	lastIndex = where(corrVec == 0)[0].min()-1
	tVec = tVec[:lastIndex]
	normVec = normVec[:lastIndex]
	corrVec = corrVec[:lastIndex]

	return tVec, normVec, corrVec


