import sys
import os

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


import utilities
reload(utilities)

import serialization
reload(serialization)

execfile(__path__[0] + "/Distribution.py")
execfile(__path__[0] + "/Enum.py")
execfile(__path__[0] + "/Potential.py")
execfile(__path__[0] + "/Problem.py")

execfile(__path__[0] + "/CreateInstance.py")
execfile(__path__[0] + "/Config.py")
execfile(__path__[0] + "/Plot.py")
execfile(__path__[0] + "/Utility.py")
execfile(__path__[0] + "/Redirect.py")
execfile(__path__[0] + "/Interrupt.py")

#Load propagators
execfile(__path__[0] + "/propagator/init.py")


