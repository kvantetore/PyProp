#Import system modules
import sys
import os
from pylab import *
from numpy import *
import tables

#Load pyprop
sys.path.insert(1, os.path.abspath("./pyprop"))
import pyprop
pyprop = reload(pyprop)
pyprop.ProjectNamespace = locals()

#Load the project module
from libpotential import *

#load the other python files in this project
import potential_data
import spline
execfile("constants.py")
execfile("delay_scan.py")
execfile("initialization.py")
execfile("serialization.py")
execfile("potential.py")

def IPython1Test(host, port):
	
	ctrl.execute("all", 'execfile("example.py")', block=True)
	r = r_[0:10]
	ctrl.scatterAll("delayList", r)
	print ctrl.execute("all", "print delayList", block=True)

def MakePhaseShiftPlot(**args):
	
	args["pulsePhase"] = 0
	t, c, r, rad = Propagate(**args)

	args["pulsePhase"] = pi/2.
	t, c2, r, rad = Propagate(**args)

	args["pulsePhase"] = pi/4.
	t, c3, r, rad = Propagate(**args)

	colors = "bgrcmykw"
	for i in range(len(colors)):
		plot(t, c[:,i], "%c-" % colors[i])
		plot(t, c2[:,i], "%c--" % colors[i])
		plot(t, c3[:,i], "%c." % colors[i])

	return t, c, c2

def Propagate(**args):
	args['config'] = "config.ini"
	args['imtime'] = False
	inputfile = GetInputFile(**args)
	prop = SetupProblem(**args)

	f = tables.openFile(inputfile, "r")
	try:
		#Setup initial state
		prop.psi.GetData()[:,0] = f.root.initial_state[:]
		prop.psi.GetData()[:,1] = 0
		initPsi = prop.psi.Copy()

		#Load eigenstates
		eigenstates = conj(f.root.eigenstates[:])
		f.close()
	finally:
		f.close()
	del f

	boundE, boundV = LoadBoundEigenstates(**args)
	contE1, contV1, contE2, contV2 = LoadContinuumEigenstates(**args)
	dE = average(diff(contE2))
	E = r_[-0.5:0:dE]
	
	if len(eigenstates.shape) != 2:
		raise Exception("Please implement for Rank!=2")
	eigenstateSize = len(eigenstates[0,:])
	weight = prop.psi.GetRepresentation().GetRepresentation(0).GetRange(0).Dx

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	timeList = []
	corrList = []
	radialData = []

	outputCount = 500
	if "outputCount" in args:
		outputCount = args["outputCount"]

	tPrev = 0
	for t in prop.Advance(outputCount):
		corrList += [abs(GetEigenstateCorrelations(prop, eigenstates))**2]
		radialData += [sum(abs(prop.psi.GetData())**2, axis=1)]
		
		norm = prop.psi.GetNorm()
		corr = abs(prop.psi.InnerProduct(initPsi))**2
		timeList.append(t/femtosec_to_au)
		if t-tPrev > 20*femtosec_to_au:
			print "t = %f, N = %f, Corr = %.17f" % (t/femtosec_to_au, norm, corr) 
			tPrev = t
		#plot(r, abs(prop.psi.GetData())**2)

	timeList.append(prop.PropagatedTime / femtosec_to_au)
	corrList += [abs(GetEigenstateCorrelations(prop, eigenstates))**2]

	norm = prop.psi.GetNorm()
	corr = abs(prop.psi.InnerProduct(initPsi))**2
	print "t = %f, N = %f, Corr = %.17f" % (t/femtosec_to_au, norm, corr) 

	energyDistrib1, energyDistrib2 = CalculateEnergyDistribution(prop.psi.GetData(), E, contE1, contV1, contE2, contV2)

	return array(timeList), array(corrList), E, array([energyDistrib1, energyDistrib2])



def GetEigenstateCorrelations(prop, eigenstates):
	"""
	Returns the correlation betweeen the propagated state and the eigenstates
	"""
	eigenstateSize = eigenstates.shape[1]
	weight = prop.psi.GetRepresentation().GetRepresentation(0).GetRange(0).Dx

	psiSlice = prop.psi.GetData()[0:eigenstateSize,0]

	return dot(eigenstates, psiSlice)*weight

	
def PlotElectricField(**args):
	conf = SetupConfig(**args)
	psi = pyprop.CreateWavefunction(conf)
	pot = pyprop.CreatePotentialFromSection(conf.ElectronicCoupling, "ElectronicCoupling", psi)
	t = r_[0:conf.Propagation.duration:abs(conf.Propagation.timestep)]
	field = array([pot.GetTimeValue(curT) for curT in t])
	plot(t/femtosec_to_au, field)

	return t, field


def PlotEigenstate(solver, index):
	prop = solver.BaseProblem
	solver.SetEigenvector(prop.psi, index)

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	data = prop.psi.GetData()
	plot(r, abs(data)**2 -0.6)
