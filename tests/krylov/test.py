import unittest
import sys
import os.path

from numpy import sort

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
import pyprop.modules.solvers.krylov.piram as piram
import pyprop.modules.solvers.krylov.pamp as pamp
import pyprop.modules.solvers.krylov.gmres as gmres

class PampTestCase(unittest.TestCase):
	def runTest(self):
		#Setup propagator
		conf = pyprop.config.Load("config.ini")
		prop = pyprop.problem.Problem(conf)
		prop.SetupStep()
		
		#propagate with imaginary time
		for t in prop.Advance(10):
			E = prop.GetEnergy()
			print "t = %s, E = %s" % (t, E)

		#Make sure we got the right energy
		E = prop.GetEnergy()
		assert(abs(E - 0.5) < 0.01)
		

class PiramTestCase(unittest.TestCase):		
	def runTest(self):
		#Setup propagator
		conf = pyprop.config.Load("config.ini")
		prop = pyprop.problem.Problem(conf)
		prop.SetupStep()
		
		solver = piram.PiramSolver(prop)
		solver.Solve()
		
		E = sort(solver.GetEigenvalues())
		assert(abs(E[0] - 0.5) / E[0] < 0.01)
		assert(abs(E[1] - 1.5) / E[1] < 0.01)
		assert(abs(E[2] - 2.5) / E[2] < 0.01)
		


class GmresTestCase(unittest.TestCase):		
	def runTest(self):
		#Setup propagator
		conf = pyprop.config.Load("config.ini")
		prop = pyprop.problem.Problem(conf)
		prop.SetupStep()
		
		#Create Shift Invert solver
		shift = 2.5
		conf.GMRES.shift = shift
		shiftInvertSolver = gmres.GMRESShiftInvertSolver(prop)
		conf.Arpack.inverse_iterations = True
		conf.Arpack.matrix_vector_func = shiftInvertSolver.InverseIterations
		conf.Arpack.krylov_eigenvalue_shift = shift
		
		#Create eigenvalue solver
		solver = piram.PiramSolver(prop)
		solver.Solve()
		
		#Calculate the real eigenvalues from the shiftinverted eigenvalues
		shiftedE = solver.GetEigenvalues()
		E = sort(map(lambda Ei: (1.0 / Ei) + shift, shiftedE))
		
		assert(abs(E[0] - 1.5) / E[0] < 0.01 and "E = %s" % E[0])
		assert(abs(E[1] - 2.5) / E[1] < 0.01 and "E = %s" % E[1])
		assert(abs(E[2] - 3.5) / E[2] < 0.01 and "E = %s" % E[2])
		


if __name__ == "__main__":
	unittest.main()

