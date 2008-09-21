import sys
import os
import time

import pylab
import numpy
from numpy import array, r_, double, complex, exp

sys.path.insert(0, "./pyprop")
import pyprop
pyprop = reload(pyprop)

from libpotential import *

import time

execfile("../common/benchmark.py")
execfile("../common/initialization.py")
execfile("../common/basic-propagation.py")
execfile("../common/grid-generation.py")
execfile("h2p.py")

def SetupConfig(**args):
	conf = CommonSetupConfig(**args)

	if "differenceOrder" in args:
		differenceOrder = args["differenceOrder"]
		conf.SetValue("ElectronPropagator", "difference_order", differenceOrder)

	if "polarCount" in args:
		polarCount = args["polarCount"]
		conf.SetValue("ElectronRepresentation", "rank0", [0, 2*pi, polarCount])

	return conf


def Propagate(**args):
	initPsi = args.get("initPsi", None)
	
	args["imtime"] = False
	args["potentialList"] = ["LaserPotential"]
	args["duration"] = 40

	prop = SetupProblem(**args)
	if initPsi != None: 
		prop.psi.GetData()[:] = initPsi.GetData()
	else:
		prop.psi.Normalize()
		initPsi = prop.psi.Copy()

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	#phi = prop.psi.GetRepresentation().GetLocalGrid(1)

	curtime = - time.time()
	for t in prop.Advance(20):
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		norm = prop.psi.GetNorm()
		
		#pcolormesh(phi/pi, r, abs(prop.psi.GetData())**2, vmax=0.05)
		#clf()
		#imshow(log(abs(prop.psi.GetData())**2))
		#draw()


		#errorEstimate = prop.Propagator.PampWrapper.GetPropagationErrorEstimate()
		errorEstimate = 0
		
		if pyprop.ProcId == 0:
			print "t = %03.2f, corr(t) = %1.8f, N = %01.8f, Err = %s" % (t, corr, norm, errorEstimate)

	curtime += time.time()	
	if pyprop.ProcId == 0: print "Duration = %s" % (curtime)

	return prop

def GetH2Groundstate(**args):
	assert("configFile" not in args)
	args["configFile"] = "config.ini"
	args["eigenvalueCount"] = 1
	
	prop = SetupProblem(**args)
	solver = pyprop.ArpackSolver(prop)
	solver.Solve()

	E = solver.GetEigenvalues()[0]
	solver.SetEigenvector(prop.psi, 0)

	return E, prop.psi



def SetupH2MaskPotential(psi, cutoffDistance):
	class cutoffSection(object):
		def __init__(self, cutoffDistance):
			self.type = PotentialType.Static
			self.classname = "H2MaskPotential"
			self.cutoff_distance = cutoffDistance

	section = cutoffSection(cutoffDistance)		
	potential = pyprop.CreatePotentialFromSection(section, "H2MaskPotential", psi)

	return potential



