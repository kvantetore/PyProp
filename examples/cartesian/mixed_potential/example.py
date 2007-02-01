import pyprop
pyprop = reload(pyprop)
from libMixedPotential import *

conf = pyprop.Load("config.ini")
prop = pyprop.Problem(conf)
prop.SetupStep()

def run():
	for t in prop.Advance(10):
		print "t = ", t

