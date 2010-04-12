import unittest
import sys
import os.path

#import pyprop
sys.path.append(os.path.join("..", ".."))
import pyprop.config
import pyprop.problem

#referenced from config file
import pyprop.core as core
import pyprop.propagator.basispropagator as basispropagator
import pyprop.modules.discretizations.bspline as bspline
import pyprop.modules.potentials.oscillator as oscillator
import pyprop.modules.potentials.tensorpotentialbase as tensorpotentialbase

class MultiplyHamiltonianTestCase(unittest.TestCase):
	def runTest(self):
		conf = pyprop.config.Load("config.ini")
		prop = pyprop.problem.Problem(conf)
		prop.SetupStep()
		
		srcPsi = prop.psi.Copy()
		dstPsi = prop.psi.Copy()
		dstPsi.Clear()
		
		prop.MultiplyHamiltonian(srcPsi, dstPsi)

if __name__ == "__main__":
	unittest.main()

