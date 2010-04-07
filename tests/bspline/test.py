import unittest
import sys
import os.path

#import pyprop
sys.path.append(os.path.join("..", ".."))
import pyprop.config
import pyprop.problem

#referenced from config file
import pyprop.core as core
import pyprop.modules.discretizations.bspline as bspline
import pyprop.propagator.combined as combined
import pyprop.modules.potentials.oscillator as oscillator

class PropagatorTestCase(unittest.TestCase):
	def runTest(self):
		conf = pyprop.config.Load("config.ini")
		prop = pyprop.problem.Problem(conf)
		prop.SetupStep()
		prop.AdvanceStep()
		for t in prop.Advance(10):
			E = prop.GetEnergy()
			print t, prop.GetEnergy()
			
		assert(abs(E-0.5) < 0.01)

	
if __name__ == "__main__":
	unittest.main()

