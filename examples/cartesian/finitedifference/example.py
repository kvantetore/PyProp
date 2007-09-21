import sys
import os

from pylab import *
from numpy import *

pyprop_path = "../../../"
sys.path.insert(1, os.path.abspath(pyprop_path))
import pyprop
pyprop = reload(pyprop)

def FindGroundstate():
	#load config
	conf = pyprop.Load("find_groundstate.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#propagate to find ground state
	for t in prop.Advance(10):
		print "t = ", t, ", E =", prop.GetEnergy()
	
	#save groundstate to disk
	prop.SaveWavefunctionHDF("groundstate.h5", "/wavefunction")

	#Find energy
	E1 = prop.GetEnergyImTime()
	E2 = prop.GetEnergyExpectationValue()
	print "Groundstate energy:\n\t %s a.u.\n\t %s" % (E1, E2)
	pyprop.Plot1D(prop)

	return prop


	
