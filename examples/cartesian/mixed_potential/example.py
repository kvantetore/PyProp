import os
import sys

sys.path.insert(1, "../../../")
import pyprop
pyprop = reload(pyprop)

from libMixedPotential import *

conf = pyprop.Load("config.ini")
prop = pyprop.Problem(conf)
prop.SetupStep()
prop.psi.Normalize()
initPsi = prop.psi.Copy()

def run():
	for t in prop.Advance(10):
		corr = abs(initPsi.InnerProduct(prop.psi))**2
		print "t = ", t, ", corr = ", corr

