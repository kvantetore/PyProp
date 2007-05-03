#improt system modules
import sys
import os

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("../../.."))

#Load and reload pyprop in order to get recent changes
import pyprop
pyprop = reload(pyprop)
from libqdot import *

#numpy an pylab for good measure
import pylab
from numpy import *

def FindGroundstate():
	prop = SetupProblem()
	for t in prop.Advance(10):
		print "t = ", t, ", E = ", prop.GetEnergy()

def Propagate():
	prop = SetupProblem()
	prop.psi.Normalize()
	x = prop.psi.GetRepresentation().GetLocalGrid(0)
	for t in prop.Advance(50):
		print "t = ", t, ", Norm = ", prop.psi.GetNorm()
		data = abs(prop.psi.GetData())**2
		pylab.plot(x, data)

def SetupProblem(**args):
	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()	
	return prop

def SetupConfig(**args):
	conf = pyprop.Load("config.ini")
	return conf

