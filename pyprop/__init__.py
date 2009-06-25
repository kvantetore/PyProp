import sys
import os
import time

#import mpi
try:
	import pylab
	from pylab import *
except:
	print "Warning: Unable to load matplotlib. Plotting will not be available"

try:
	import tables
except:
	print "Warning: Could not load module tables (pytables). Serialization to HDF files will not be available"

import numpy
from numpy import *

#The main namespace used for accessing functions defined in the "main"
#project (i.e. example.py) file from configuration files and the like. 
#Must be overridden if the main project is not in the __main__ namespace 
#(i.e. when using IPython1)
ProjectNamespace = sys.modules["__main__"].__dict__

__DisableMPI = False
if 'PYPROP_SINGLEPROC' in os.environ:
	__DisableMPI = bool(os.environ['PYPROP_SINGLEPROC'])

ProcId = 0
ProcCount = 1
if __DisableMPI:
	print "MPI Disabled"
else:
	try:
		import pypar
		ProcId = pypar.rank()
		ProcCount = pypar.size()
	except:
		print "Warning: unable to load mpi."
	
import core
reload(core)
if __DisableMPI:
	for name, obj in core.__dict__.iteritems():
		if name.startswith("DistributedModel_"):
			obj.ForceSingleProc()
	
StaticStorageModel = core.StaticPotential_1.StorageModel


DEBUG_PRINT_MEMORY_USAGE = False

import utilities
reload(utilities)

import serialization
reload(serialization)

import plotting
reload(plotting)

execfile(__path__[0] + "/Distribution.py")
execfile(__path__[0] + "/Enum.py")
execfile(__path__[0] + "/Potential.py")
execfile(__path__[0] + "/Absorber.py")
execfile(__path__[0] + "/Problem.py")

execfile(__path__[0] + "/CreateInstance.py")
execfile(__path__[0] + "/Config.py")
execfile(__path__[0] + "/Plot.py")
execfile(__path__[0] + "/Utility.py")
execfile(__path__[0] + "/Redirect.py")
execfile(__path__[0] + "/Interrupt.py")
execfile(__path__[0] + "/Timer.py")

execfile(__path__[0] + "/BasisExpansion.py")
execfile(__path__[0] + "/CoupledSphericalHarmonics.py")

execfile(__path__[0] + "/bspline/BSpline.py")

execfile(__path__[0] + "/GridGeneration.py")

#Load propagators
execfile(__path__[0] + "/propagator/init.py")
execfile(__path__[0] + "/solver/init.py")
execfile(__path__[0] + "/tensorpotential/init.py")


if hasattr(core, "PapiSetup"):
	core.PapiSetup()

def GetMemoryUsage():
	if hasattr(core, "PapiGetDynamicMemoryInfo"):
		mem = core.PapiGetDynamicMemoryInfo()
	return mem.size, mem.peak


def SerialPrint(str, proc=-1):
	if ProcCount == 1:
		print str
	else:
		if proc==-1: procList = range(ProcCount)
		else: procList = [proc]
		for i in procList:
			if i == ProcId:
				print "Proc %4i: %s" % (ProcId, str,)
			sys.stdout.flush()	
			pypar.barrier()

def PrintMemoryUsage(header, proc=-1):
	if DEBUG_PRINT_MEMORY_USAGE:
		mem, peak = GetMemoryUsage()
		SerialPrint(header, proc=0)
		SerialPrint("Memory Usage = %10iKB, Peak = %10iKB" % (mem, peak), proc=proc)



