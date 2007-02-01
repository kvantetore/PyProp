#system modules
import sys
import os

#add pyprop to our search path
PYPROP_PATH = os.path.abspath('../../../')
sys.path.insert(1, PYPROP_PATH)

#import and reload pyprop
import pyprop
pyprop = reload(pyprop)

#pylab and numpy should be imported in this order.
from pylab import *
from numpy import *

def test():
	conf = pyprop.Load("config.ini")
	prop = pyprop.Problem(conf)
	prop.SetupStep()
	return prop

def FindGroundstateEnergy(transformType, N, dt):
	conf = pyprop.Load("config.ini")
	conf.Representation.n = N
	conf.Representation.transform_type = transformType
	conf.Propagation.timestep = abs(dt) * -1.0j
	
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	for t in prop.Advance(1): pass
	return prop.GetEnergy()

def CalculateEnergyConvergenceGridSize(transformType, dt=0.01):
	N = r_[4:40:4]
	E = zeros(len(N), dtype=double)
	for i in range(0,len(N)):
		E[i] = FindGroundstateEnergy(transformType, N[i], dt)
		
	return N, E

def CalculateEnergyConvergenceTimeStep(transformType, N=16):
	dt = r_[0.001:1:0.01]
	E = zeros(len(dt), dtype=double)
	for i in range(0,len(dt)):
		E[i] = FindGroundstateEnergy(transformType, N, dt[i])
		
	return dt, E
	

def run():
	prop = test()
	prop.psi.Normalize()
	initPsi = prop.psi.Copy()
	
	for t in prop.Advance(20):
		plot(abs(prop.psi.GetData())**2)
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		print "t = ", t,", corr = ", corr, "norm = ", prop.psi.GetNorm()
		print "E =", prop.GetEnergy()

	return prop
