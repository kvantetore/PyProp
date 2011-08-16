#import system modules
import sys
import os

from numpy import double, abs, exp, pi, log10, sqrt
import pylab

#Make sure we use the correct pyprop library
sys.path.insert(1, os.path.abspath("../../.."))
import pyprop

import pyprop
if pyprop.IsRunningFromSource:
	sys.path.append(os.path.join(pyprop.BuildPath, "../examples", "tensor",
		"hydrogen_1d"))
import libhydrogen1d
from libhydrogen1d import KineticEnergyPotential_1, \
	RegularizedCoulombPotential_1, ManolopoulosAbsorber_1

import pyprop.config
import pyprop.problem

#referenced from config file
import pyprop.modules.discretizations.bspline as bspline
import pyprop.modules.discretizations.finitedifference as finitedifference
import pyprop.core as core
import pyprop.customgrid as customgrid
import pyprop.modules.solvers.krylov.pamp as pamp
import pyprop.propagator.basispropagator as basispropagator


def FindGroundstate():
	conf = pyprop.config.Load("config.ini")
	prop = pyprop.problem.Problem(conf)
	prop.SetupStep()

	for t in prop.Advance(5):
		print "t = ", t, " E = ", prop.GetEnergyExpectationValue()

	return prop


def Propagate():
	conf = pyprop.config.Load("config.ini")
	prop = pyprop.problem.Problem(conf)
	prop.SetupStep()
	prop.psi.Normalize()

	xgrid = prop.psi.GetRepresentation().GetLocalGrid(0)

	initialPsi = prop.psi.Copy()
	fig = pylab.figure()
	ax1 = fig.add_subplot(311)
	ax2 = fig.add_subplot(312)
	ax3 = fig.add_subplot(313)
	pylab.show()
	probList = []
	timeList = []
	probDensityList = []

	for t in prop.Advance(100):
		corr = abs(prop.psi.InnerProduct(initialPsi))
		prob = prop.psi.GetNorm()**2
		print "t = %f, N(t)**2 = %.15f, C(t) = %.15f" % (t, prob, corr)

		probDensity = abs(prop.psi.GetData())**2
		probDensityList.append(probDensity)
		ax1.plot(xgrid, log10(probDensity))
		ax2.plot(xgrid, abs(probDensity))
		pylab.draw()

		timeList.append(t)
		probList.append(prob)

		if prob < 1e-26:
			break

	ax3.plot(timeList, probList)
	pylab.draw()

	return prop, xgrid, probDensityList
