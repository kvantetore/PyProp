import sys
import os
import profile

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
	prop.SaveWavefunction("groundstate.dat")

	#Find energy
	print "Groundstate energy:", prop.GetEnergy(), "a.u."
	pyprop.Plot1D(prop)

	return prop


	
