import sys
import os
import time

import pylab
import numpy

sys.path.insert(0, "./pyprop")
import pyprop
pyprop = reload(pyprop)

from libpotential import *

execfile("../common/benchmark.py")
execfile("../common/initialization.py")
execfile("../common/basic-propagation.py")
execfile("../common/grid-generation.py")

def SetupConfig(**args):
	conf = CommonSetupConfig(**args)

	if "differenceOrder" in args:
		differenceOrder = args["differenceOrder"]
		conf.SetValue("ElectronPropagator", "difference_order", differenceOrder)

	if "polarCount" in args:
		polarCount = args["polarCount"]
		conf.SetValue("ElectronRepresentation", "rank0", [0, 2*pi, polarCount])

	return conf

def Propagate(**args):
	initPsi = args.get("initPsi", None)
	
	args["imtime"] = False
	#args["potentialList"] = ["LaserPotential", "AbsorbingPotential"]
	args["duration"] = 40

	prop = SetupProblem(**args)
	if initPsi != None: 
		prop.psi.GetData()[:] = initPsi.GetData()
	else:
		initPsi = prop.psi.Copy()

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	#phi = prop.psi.GetRepresentation().GetLocalGrid(1)

	for t in prop.Advance(20):
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		norm = prop.psi.GetNorm()
		
		#pcolormesh(phi/pi, r, abs(prop.psi.GetData())**2, vmax=0.05)
		#imshow(abs(prop.psi.GetData())**2, vmax=0.05)

		print "t = %03.2f, corr(t) = %1.8f, N = %01.8f" % (t, corr, norm)

	return prop


def FindElectronicEnergyCurves(**args):
	args["silent"] = True

	nuclearSeparationList = \
		array([0.4000, 0.4500, 0.5000, 0.5500,\
 		0.6500, 0.7000, 0.7500, 0.8000, 0.9000, \
		1.0000, 1.1000, 1.2000, 1.3000, 1.3500, \
		1.3900, 1.4000, 1.4010, 1.4011, 1.4100, \
		1.4500, 1.5000, 1.6000, 1.7000, 1.8000, \
		1.9000, 2.0000, 2.1000, 2.2000, 2.3000, \
		2.4000, 2.5000, 2.6000, 2.7000, 2.8000, \
		2.9000, 3.0000, 3.1000, 3.2000, 3.3000, \
		3.4000, 3.5000, 3.6000, 3.7000, 3.8000, \
		3.9000, 4.0000, 4.1000, 4.2000, 4.3000, \
		4.4000, 4.5000, 4.6000, 4.7000, 4.8000, \
		4.9000, 5.0000, 5.1000, 5.2000, 5.3000, \
		5.4000, 5.5000, 5.6000, 5.7000, 5.8000, \
		5.9000, 6.0000, 6.1000, 6.2000, 6.3000, \
		6.4000, 6.5000, 6.6000, 6.7000, 6.8000, \
		6.9000, 7.0000, 7.2000, 7.4000, 7.6000, \
		7.8000, 8.0000, 8.2500, 8.5000, 9.0000, \
		9.5000, 10.0000], dtype=double)

	nuclearSeparationList = nuclearSeparationList[::5]

	energyList = []

	#Set up a temporary psi
	prop = SetupProblem(**args)
	tempPsi = prop.psi
	initialResidual = tempPsi.Copy()

	#Solve eigenstate for each nuclear separation
	estimatedTime = False
	for i, separation in enumerate(nuclearSeparationList):
		if estimatedTime:
			print "Estimated Time Remaining: %s" % formatTimedelta(estimatedTime)

		currentTime = - time.time()
		solver = FindEigenstates(nuclearSeparation=separation, initialResidual=initialResidual, **args)
		E = sort(solver.GetEigenvalues().copy())
		restartCount = solver.Solver.GetIterationCount()
		operationCount = solver.Solver.GetOperationCount()

		print "Found %i eigenvalues, used %i restart iterations and %i multiply hamiltonian operations" % (len(E), restartCount, operationCount)
		 
		energyList += [E]

		#Use a linear combination of the current eigenvectors as init residual to the next parameter value
		psi = solver.BaseProblem.psi
		initialResidual.GetData()[:] = 0
		for j, curE  in enumerate(E):
			solver.SetEigenvector(psi, j)
			initialResidual.GetData()[:] += psi.GetData()


		currentTime += time.time()
		estimatedTime = currentTime * (len(nuclearSeparationList) - i - 1)


		
	return nuclearSeparationList, energyList


def DetectPossibleCrossings(nuclearSeparationList, energyList, maximumDistance):
	"""
	Detect possible crossings in energyList. A possible crossing is where the
	slope of difference = E_i+1 - E_i changes sign from negative to positive
	(note that difference is always positive), and difference is less than 
	maximumDistance

	returns a list of tuples [(separationIndex, stateIndex), ...] containing the
	index in the separationList and the stateIndex where the possible crossing is 
	located
	"""

	energyArray = array(energyList)
	energyCount = energyArray.shape[1]
	paramCount = len(nuclearSeparationList)

	assert(paramCount == energyArray.shape[0])

	possibleCrossings = []

	#loop over all states, 
	for state in range(energyCount-1):
		#energyArray is sorted, so curCurve < nextCurve for all indices
		# => difference > 0
		curCurve = energyArray[:,state]
		nextCurve = energyArray[:, state+1]
		difference = nextCurve - curCurve 

		#Derivative of the difference
		differenceDeriv = diff(difference) / diff(nuclearSeparationList)
		for i in range(paramCount-2):
			if differenceDeriv[i] < 0 and differenceDeriv[i+1] > 0 and difference[i+1] < maximumDistance:
				possibleCrossings += [(i+1, state)]

	return possibleCrossings

try:
	import ctypesGsl
	import ctypesGsl.minim as minim
except:
	print "Warning: could not load ctypesGsl"

def ResolveCrossing(nuclearSeparationList, energyList, possibleCrossing, accuracy, maximumDistance, **args):
	"""
	Resolve a possible crossing as returned from detect possible crossings, 
	by searching for a minimum between the difference of the two energy curves. 
	if the minimum is > maximumDistance, it is considered an avoided crossing,
	and the method returns false.

	Returns true if it is a crossing, and false if it is not
	"""
	
	energyArray = array(energyList)
	parameterIndex, state = possibleCrossing

	def CalculateDifference(nuclearSeparation):
		if nuclearSeparation in nuclearSeparationList:
			idx = nuclearSeparationList.index(nuclearSeparation)
			return energyList[idx][state+1] - energyList[idx][state]
			
		#Calculate eigenvalues
		solver = FindEigenstates(nuclearSeparation=nuclearSeparation, silent=True, **args)
		E = sort(solver.GetEigenvalues().copy())

		#Store eigenvalues
		nuclearSeparationList.append(nuclearSeparation)
		energyList.append(E)
		
		#Calculate separation between energy curves
		difference = E[state+1] - E[state]

		return difference

	#Initialize minimization routine
	gslFunc = ctypesGsl.gsl_function(CalculateDifference)
	minimizer = minim.min_fminimizer(minim.min_fminimizer_brent, gslFunc)
	minimizer.init(nuclearSeparationList[parameterIndex], nuclearSeparationList[parameterIndex-1], nuclearSeparationList[parameterIndex+1])

	#Iterate until we know it is a crossing, or 
	converged = False
	while not converged:
		minimizer.iterate()
		difference = ctypesGsl.libgsl.gsl_min_fminimizer_f_minimum(minimizer.ptr)
		print "min difference = %f" % difference
		if difference < maximumDistance:
			converged = True
		if minimizer.test_interval(accuracy, accuracy):
			converged = True

	return difference < maximumDistance
	
		


