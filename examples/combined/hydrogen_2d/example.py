representation0 = "LinearGrid"
import sys
import os
import time

import pylab
import numpy

sys.path.insert(0, "./pyprop")
import pyprop
pyprop = reload(pyprop)

from libpotential import *

execfile("benchmark.py")
execfile("initialization.py")
execfile("basic-propagation.py")
execfile("grid-generation.py")

def SetupConfig(**args):
	conf = CommonSetupConfig(**args)

	if "differenceOrder" in args:
		differenceOrder = args["differenceOrder"]
		conf.SetValue("ElectronPropagator", "difference_order", differenceOrder)

	if "polarCount" in args:
		polarCount = args["polarCount"]
		conf.SetValue("ElectronRepresentation", "rank0", [0, 2*pi, polarCount])

	return conf

def Propagate(initPsi = None, **args):
	args["imtime"] = False
	#args["potentialList"] = ["LaserPotential"]
	args["duration"] = 40

	prop = SetupProblem(**args)

	if initPsi:
		prop.psi.GetData()[:] = initPsi.GetData()
	else:
		prop.psi.Normalize()
		initPsi = prop.psi.Copy()

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	#phi = prop.psi.GetRepresentation().GetLocalGrid(1)

	for t in prop.Advance(20):
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		norm = prop.psi.GetNorm()
		
		#pcolormesh(phi/pi, r, abs(prop.psi.GetData())**2)
		#imshow(abs(prop.psi.GetData())**2, vmax=0.05)
		
		if pyprop.ProcId == 0:
			print "t = %03.2f, corr(t) = %1.8f, N = %01.15f" % (t, corr, norm)

	return prop

