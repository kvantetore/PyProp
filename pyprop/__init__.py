#import mpi
import pylab
import numpy

from pylab import *
from numpy import *

#import mpi.pympi as pympi
import cPickle as pickle

import core
reload(core)

import utilities
reload(utilities)
from utilities import *
del utilities

execfile(__path__[0] + "/Distribution.py")
execfile(__path__[0] + "/Enum.py")
execfile(__path__[0] + "/Potential.py")
execfile(__path__[0] + "/Problem.py")

execfile(__path__[0] + "/CreateInstance.py")
execfile(__path__[0] + "/Config.py")
execfile(__path__[0] + "/Serialize.py")
execfile(__path__[0] + "/Plot.py")
execfile(__path__[0] + "/Utility.py")

execfile(__path__[0] + "/Propagator.py")
execfile(__path__[0] + "/propagator/CartesianPropagator.py")
execfile(__path__[0] + "/propagator/CartesianMixedPropagator.py")
execfile(__path__[0] + "/propagator/ExponentialFiniteDifferencePropagator.py")
execfile(__path__[0] + "/propagator/SphericalPropagator.py")
execfile(__path__[0] + "/propagator/TransformedGridPropagator.py")


#set up ProcId and ProcCoun. if pympi is not imported, 
#we are on a single process"
ProcId = 0
ProcCount = 1
try:
	ProcId = pympi.rank
	ProcCount = pympi.size
except: pass
