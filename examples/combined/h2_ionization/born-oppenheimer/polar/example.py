import sys
import os
import time

import pylab
import numpy

sys.path.insert(0, "./pyprop")
import pyprop
pyprop = reload(pyprop)

from libpotential import *

execfile("../common/benchmark.py")
execfile("../common/initialization.py")
execfile("../common/basic-propagation.py")
execfile("../common/grid-generation.py")

def SetupConfig(**args):
	conf = CommonSetupConfig(**args)

	if "bsplineOrder" in args:
		bsplineOrder = args["bsplineOrder"]
		conf.SetValue("ElectronRadialRepresentation", "order", bsplineOrder)

	if "bsplineCount" in args:
		bsplineCount = args["bsplineCount"]
		conf.SetValue("ElectronRadialRepresentation", "xsize", bsplineCount)

	if "differenceOrder" in args:
		differenceOrder = args["differenceOrder"]
		conf.SetValue("ElectronRadialPropagator", "difference_order", differenceOrder)

	if "polarCount" in args:
		polarCount = args["polarCount"]
		conf.SetValue("ElectronPolarRepresentation", "rank0", [0, 2*pi, polarCount])

	return conf

def Propagate(**args):
	initPsi = args["initPsi"]
	
	args["imtime"] = False
	#args["potentialList"] = ["LaserPotential", "AbsorbingPotential"]
	args["duration"] = 40

	prop = SetupProblem(**args)
	prop.psi.GetData()[:] = initPsi.GetData()

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	#phi = prop.psi.GetRepresentation().GetLocalGrid(1)

	for t in prop.Advance(20):
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		norm = prop.psi.GetNorm()
		
		#pcolormesh(phi/pi, r, abs(prop.psi.GetData())**2, vmax=0.05)
		#imshow(abs(prop.psi.GetData())**2, vmax=0.05)

		print "t = %03.2f, corr(t) = %1.8f, N = %01.8f" % (t, corr, norm)

	return prop

def test():
	prop = SetupProblem(configFile="test.ini")
	#return prop.Propagator.SubPropagators[0].Propagator.GetBackwardPropagationLapackBanded().copy()
	#return prop.Propagator.SubPropagators[0].Propagator.GetDifferenceCoefficients().copy()

	x = prop.psi.GetRepresentation().GetLocalGrid(0)
	prop.psi.GetData()[:] = sin(x)

	destPsi = prop.GetTempPsi()
	destPsi.GetData()[:] = 0

	cnPropagator = prop.Propagator.SubPropagators[0].Propagator
	cnPropagator.MultiplyKineticEnergyOperator(prop.psi, destPsi)

	plot(x, prop.psi.GetData().copy())
	plot(x, destPsi.GetData().copy())


