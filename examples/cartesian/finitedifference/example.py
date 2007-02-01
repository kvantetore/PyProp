import sys
import pyprop
import profile

from numpy import *
from pylab import *

reload(sys.modules['pyprop'])

def FindGroundstate():
	#load config
	conf = pyprop.Load("find_groundstate.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	#propagate to find ground state
	for t in prop.Advance(10):
		print "t = ", t
	
	#save groundstate to disk
	prop.SaveWavefunction("groundstate.dat")

	#Find energy
	print "Groundstate energy:", prop.GetEnergy(), "a.u."
	pyprop.Plot1D(prop)

	return prop


	
