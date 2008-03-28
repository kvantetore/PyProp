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
execfile("load_cmap.py")

def ConcatenateHDF5(inFileList, outFile, outputMode="w"):
	"""
	Concatenates several HDF5 files to one file.
	* inFileList: list of filenames to be concatenated
	* outFile: the file to which the concatenation should be written
	* outputMode: the mode in which outFile will be opened ("w" or "a")

	All non-group nodes are copied from the input files in the order in which they are
	specified. If a node already exists in outFile, it is skipped, and a warning is printed.

	If the same node exists in several of the input files, the node being in the first file 
	containing the duplicates will be used.
	"""

	def NodeExists(f, where, path):
		"""
		Checks wheter the node path exists in where in the file f
		"""
		nodeExists = False
		try:
			node = f.getNode(where, path)
			nodeExists = True
		except tables.NoSuchNodeError: pass
		return nodeExists

	def GetParentPath(path):
		"""
		Splits path into the parent path and the name
		"""
		path = path.strip()
		if path.endswith("/"):
			path = path[:-1]
		pathList = path.split("/")
		if len(pathList) == 1:
			return None, None
		if len(pathList) == 2:
			return "/", pathList[1]
		path = "/".join(pathList[:-1])
		name = pathList[-1]
		return path, name

	def GetOrCreateGroup(f, path):
		"""
		Gets path from the file f if it exists, or creates a group with that path
		"""
		path, name = GetParentPath(path)
		if path == None:
			return f.root
		node = None
		try:
			node = f.getNode(path, name)
			if not isinstance(node, tables.Group):
				raise Exception("Path %s/%s exists but it is not a group" (path,name))
		except tables.NoSuchNodeError: pass
		if node == None:
			node = f.createGroup(path, name, createparents=True)
		return node

	#Main method
	outFile = tables.openFile(outFile, outputMode)
	try:
		for fileName in inFileList:
			inFile = tables.openFile(fileName, "r")
			try:
				for node in inFile.walkNodes(inFile.root):
					if not(isinstance(node, tables.Group)):
						name = node._v_name
						path = node._v_parent._v_pathname
						if NodeExists(outFile, path, name):
							print "Node %s/%s already exists, skipping" % (path, name)
						else:
							#print "Copying %s/%s" % (path, name)
							outParent = GetOrCreateGroup(outFile, path)
							node.copy(outParent, name)
			finally:
				inFile.close()

	finally:
		outFile.close()

def IPython1Test(host, port):
	"""
	Very simple ipython1 test
	see http://ipython.scipy.org/moin/IPython1 for more info
	on ipython1
	"""
	ctrl.execute("all", 'execfile("example.py")', block=True)
	r = r_[0:10]
	ctrl.scatterAll("delayList", r)
	print ctrl.execute("all", "print delayList", block=True)

def MakePhaseShiftPlot(**args):
	"""
	Runs a problem with 3 different phases on the pulse,
	and makes
	"""
	
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

	conf = SetupConfig(**args)
	prop = SetupProblem(**args)

	pumpFrequency = conf.InitialCondition.pump_frequency
	pumpTransitionProbability = conf.InitialCondition.pump_probability
	pumpCount = conf.InitialCondition.pump_count
	pumpStateEnergy = GetInitialStateEnergy(inputfile)
	pumpTimes, pumpProb, pumpPhase = GetPump(pumpFrequency, pumpCount, pumpTransitionProbability, pumpStateEnergy)
	
	LoadInitialState(prop, **args)
	initPsi = prop.psi.Copy()

	boundE, boundV = LoadBoundEigenstates(**args)
	contE1, contV1, contE2, contV2 = LoadContinuumEigenstates(**args)
	dE = average(diff(contE2))
	E = r_[-0.5:0:dE]

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	dr = diff(r)[0]
	timeList = []
	corrList = []
	radialData = []


	outputCount = 500
	if "outputCount" in args:
		outputCount = args["outputCount"]

	def output():
		"""
		Performs one output from the propagation. We have many 
		advance loops, and creating an inner function like this is a 
		nice way of grouping the output functionality
		"""
		corrList.append(abs(dot(boundV, prop.psi.GetData()[:,0]))**2)
		radialData.append(prop.psi.GetData().copy())
		
		timeList.append(t/femtosec_to_au)

	#output initial state
	t = 0
	output()

	#
	fullDuration = prop.Duration
	outputFrequency = float(outputCount) / prop.Duration

	#Propagate to the end of the pump pulse
	prop.psi.GetData()[:] = 0
	for i in range(len(pumpTimes)):
		#Add a franck-condon with the correct phase wavepacket to our wavefunction.
		prop.psi.GetData()[:] += pumpProb[i] * pumpPhase[i] * initPsi.GetData()
		if i<len(pumpTimes)-1:
			prop.Duration = pumpTimes[i+1]
			curOutCount = outputFrequency * (pumpTimes[i+1] - pumpTimes[i])
			for t in prop.Advance(curOutCount):
				outputCount -= 1
				output()

	prop.psi.Normalize()

	if "removeStates" in args:
		if args["removeStates"] == "odd":
			states = r_[1:len(boundE):2]
		elif args["removeStates"] == "even":
			states = r_[0:len(boundE):2]
		else:
			raise Exception("Invalid removeStates %s" % args["removeStates"])
		
		data = prop.psi.GetData()[:,0].copy()
		for i in states:
			print "Removing %i" % i
			data -= dot(conj(boundV[i, :]), data) * boundV[i,:]/dr
		prop.psi.GetData()[:,0] = data

		plot(r, abs(prop.psi.GetData()[:,0].copy())**2)


	prop.Duration = fullDuration
	tPrev = 0
	for t in prop.Advance(outputCount):
		output()
		if t-tPrev > 20*femtosec_to_au:
			corr = abs(prop.psi.InnerProduct(initPsi))**2
			norm = prop.psi.GetNorm()
			print "t = %f, N = %f, Corr = %.17f" % (t/femtosec_to_au, norm, corr) 
			tPrev = t


	output()

	norm = prop.psi.GetNorm()
	corr = abs(prop.psi.InnerProduct(initPsi))**2
	print "t = %f, N = %f, Corr = %.17f" % (t/femtosec_to_au, norm, corr) 

	energyDistrib1, energyDistrib2 = CalculateEnergyDistribution(prop.psi.GetData(), E, contE1, contV1, contE2, contV2)

	return array(timeList), array(corrList), E, array([energyDistrib1, energyDistrib2]), array(r), array(radialData)


	
def PlotElectricField(**args):
	conf = SetupConfig(**args)
	psi = pyprop.CreateWavefunction(conf)
	pot = pyprop.CreatePotentialFromSection(conf.ProbePulsePotential, "ProbePulsePotential", psi)
	pot2 = pyprop.CreatePotentialFromSection(conf.ControlPulsePotential, "ProbePulsePotential", psi)
	t = r_[0:conf.Propagation.duration:abs(conf.Propagation.timestep)]
	field = array([pot.GetTimeValue(curT) + pot2.GetTimeValue(curT) for curT in t])
	plot(t/femtosec_to_au, field)

	return t, field


def PlotEigenstate(solver, index):
	prop = solver.BaseProblem
	solver.SetEigenvector(prop.psi, index)

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	data = prop.psi.GetData()
	plot(r, abs(data)**2 -0.6)
