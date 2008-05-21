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
execfile("figures.py")

try:
	import pyprop.plotting as myplot
except:
	print "pyprop.plotting not available"

	

def PlotTransitionMatrixElements(**args):
	matrix = real(GetTransitionMatrixElements(**args))
	figure()
	left = r_[0:matrix.shape[0]]
	bar(left-0.1125, diagonal(matrix), width=0.25)
	bar(left[:-1]+0.1125, diagonal(matrix,1), width=0.25, color="r")
	bar(left[1:]-0.25-0.1125, diagonal(matrix,-1), width=0.25, color="g")
	title("Matrix elements < i | c**2 | j > for |i-j| = 0, 1")
	xlabel("Vibrational states. Red is coupling to higher vibrational levels, Green is coupling to lower vib levels")
	ylabel("Coupling Strength")
	#axis((-0.5, 7.5, -.8, 2))

def GetTransitionMatrixElements(**args):
	args["silent"] = True
	args["configSilent"] = True
	prop = SetupProblem(**args)

	dr = prop.psi.GetRepresentation().GetLocalWeights(0)[0]
	coupling = GetElectronicCouplingBates(prop.psi, prop.Config.ProbePulsePotential)
	E, V = LoadRotatedBoundEigenstates(**args)

	#calculate matrix <v_i | c | v_ii>
	N = len(E)
	matrix = zeros((N, N), dtype=complex)
	for row in range(N):
		for col in range(N):
			matrix[row, col] = dot(V[row,:], coupling**2 * conj(V[col,:])) / dr

	return matrix

def MakePhaseDelayPlot(**args):
	args["silent"] = True
	args["configSilent"] = True
	args["molecule"] = "d2+"
	args["radialScaling"] = 2
	#args["pulseIntensity"] = 0.5e14
	args["pulseIntensity"] = 5e15
	args["pulseDelay"] = 300*femtosec_to_au
	args["duration"] = 305*femtosec_to_au
	args["outputCount"] = 1
	args["fastForward"] = 295*femtosec_to_au
	args["pulseShape"] = "square"

	#durationList = r_[0:7:0.2]
	#durationList = r_[0:1:0.05]
	durationList = r_[0:0.2:0.01]
	stateList = r_[0:7]
	theta = zeros((len(durationList), len(stateList)), dtype=double)
	for i, pulseDuration in enumerate(durationList):
		for j, initState in enumerate(stateList):
			args["pulseDuration"] = pulseDuration*femtosec_to_au
			args["initStates"] = [initState]
			theta[i,j] = GetPhaseDelay(**args)

	return durationList, stateList, theta

	
def GetPhaseDelay(**args):
	initStates=args["initStates"]
	if len(initStates) != 1:
		raise Exception("Please specify one init state")

	args["relativePhase"] = "follow"
	t, c, X, Y, r, psi = Propagate(**args)	

	phases = arctan2(c[-1,:].imag, c[-1,:].real)
	return phases[initStates[0]]


def MakeProbabilityDelayPlot(**args):
	args["silent"] = True
	args["configSilent"] = True
	args["molecule"] = "d2+"
	args["radialScaling"] = 2
	args["pulseIntensity"] = 0.5e14
	#args["pulseIntensity"] = 5e15
	args["pulseShape"] = "square"
	args["pulseDelay"] = 300*femtosec_to_au
	args["duration"] = 305*femtosec_to_au
	args["outputCount"] = 1
	args["fastForward"] = 295*femtosec_to_au

	#durationList = r_[0:7:0.5]
	durationList = r_[0:1:0.05]
	#durationList = r_[0:0.2:0.01]
	stateList = [4] #r_[0:7]
	probability = zeros((len(durationList), len(stateList)), dtype=double)
	for i, pulseDuration in enumerate(durationList):
		for j, initState in enumerate(stateList):
			args["pulseDuration"] = pulseDuration*femtosec_to_au
			args["initStates"] = [initState]
			probability[i,j] = GetPhaseDelay(**args)

	return durationList, stateList, probability



def GetSingleStateProbability(**args):
	initStates=args["initStates"]
	if len(initStates) != 1:
		raise Exception("Please specify one init state")

	args["relativePhase"] = "follow"
	t, c, X, Y, r, psi = Propagate(**args)	

	probability = abs(c[-1,:])**2
	return probability

def RunCheckerboardTest(**args):
	args["silent"] = True
	args["configSilent"] = True
	args["molecule"] = "d2+"
	args["radialScaling"] = 2
	args["pulseIntensity"] = 0.5e14
	args["pulseDuration"] = 5*sqrt(2)*femtosec_to_au
	args["duration"] = 350*femtosec_to_au
	args["outputCount"] = 200
	args["fastForward"] = 250*femtosec_to_au

	delayList = array([-1000, 293, 306])
	output = {}
	for delay in delayList:
		args["pulseDelay"] = delay*femtosec_to_au
		t, c, X, Y, r, psi = Propagate(**args)
		output["t"] = t
		output["corr_%i" % delay] = c
		output["psi_%i" % delay] = psi

	return output	



def CorrelationBarPlot(corr, ax=None, axisHeight=None, colors="b"):
	if ax == None:
		ax = gca()

	basisCount = len(corr)
	left = r_[:basisCount] - 0.4
	for i in range(basisCount):
		bar(left[i], corr[i], color=colors[i%len(colors)])

	if axisHeight == None:
		axisHeight = ax.axis()[3]
	axisWidth = len(corr) + 1

	yStart = 0
	yEnd = axisHeight
	xStart = -1
	xEnd = len(corr)
	ax.axis([xStart, xEnd, yStart, yEnd])


def CorrelationBarPlotPhase(corr, ax=None, axisHeight=None):
	if ax == None:
		ax = gca()
	
	basisCount = len(corr)
	left = r_[:basisCount] - 0.4
	bar(left, abs(corr)**2)

	if axisHeight == None:
		axisHeight = ax.axis()[3]
	axisWidth = len(corr) + 1

	winWidth = ax.get_window_extent().width()
	winHeight = ax.get_window_extent().height()

	scale = winWidth * axisHeight / (winHeight * axisWidth)

	yStart = - scale / (1 - scale)
	yEnd = axisHeight
	
	xStart = -1
	xEnd = len(corr)

	ax.axis([xStart, xEnd, yStart, yEnd])

	for i, curCorr in enumerate(corr[:-1]):
		if abs(corr[i])**2 > 1e-2:
			theta1 = arctan2(corr[i].imag, corr[i].real)
			theta2 = arctan2(corr[i+1].imag, corr[i+1].real)
			theta = theta1# - theta2
		
			xCenter = i# + 0.5
			yCenter = yStart / 2
			xy = [xCenter, yCenter]
		
			width = 0.7
			height = 2 * abs(yCenter) * 0.7
			ell = matplotlib.patches.Ellipse(xy=xy, width=width, height=height, facecolor="w")
			ax.add_artist(ell)
		
			dx = width/2. * cos(+pi/2 - theta)
			dy = height/2. * sin(+pi/2 - theta)
			lin = matplotlib.patches.Line2D(xdata=[xCenter, xCenter+dx], ydata=[yCenter, yCenter+dy])
			ax.add_artist(lin)
		
	draw_if_interactive()
	
	return ax.axis()
		


def MakeMovieFrame(conf, frameIndex, t, corr, r, rad, potData):
	fig = gcf()
	clf()

	probePulse = array([conf.ProbePulsePotential.time_function(conf.ProbePulsePotential, curT*femtosec_to_au) for curT in t])
	controlPulse = array([conf.ControlPulsePotential.time_function(conf.ControlPulsePotential, curT*femtosec_to_au) for curT in t])
	pulse = probePulse + controlPulse

	ax = fig.gca()
	ax.set_position([0.05,0.05,0.9,0.68])
	#ax.plot(r, potData)
	#ax.plot(r, 0.1*abs(rad[frameIndex,:,0])**2 + potData[:,0], "k")
	#ax.plot(r, 0.1*abs(rad[frameIndex,:,1])**2 + potData[:,1], "k")
	ax.plot(r, 0.1*abs(rad[frameIndex,:,0])**2 -0.6, "k")
	ax.plot(r, 0.1*abs(rad[frameIndex,:,1])**2 -0.2, "k")
	ax.axis([0,10,-0.65,-0.1])
	yticks([])

	ax = axes([0.52, 0.77, 0.43, 0.20])
	ax.plot(t, pulse, "b")
	ax.plot([t[frameIndex]], [pulse[frameIndex]], "r.", markersize=10)
	ax.axis([t[0], t[-1], 1.2*min(pulse), 1.2*max(pulse)])
	yticks([])

	basisCount = corr.shape[1]
	ax = axes([0.05, 0.77, 0.43, 0.20])
	curAxis = list(CorrelationBarPlotPhase(corr[frameIndex,:10], ax, 1.1*numpy.max(corr)))
	#ax.axis([-1, basisCount, 0, 1.1*numpy.max(corr)])
	yticks([])


def MakeMovie(**args):
	conf = SetupConfig(**args)
	prop = SetupProblem(**args)
	potData = prop.Propagator.PotentialList[0].Potential.GetPotentialData().copy()
	del prop

	interactive = rcParams["interactive"]
	rcParams["interactive"] = False

	generateFrames = args.get("generateFrames", True)
	generateMovie = args.get("generateMovie", True)
	isTest = args.get("isTest", False)

	if generateFrames:
		t, corr, X, Y, r, rad = Propagate(**args)
		
		figure(figsize=(8,8))
		for i, curT in enumerate(t):
			progressStr = "%#3i%s Complete" % (i*100 / float(len(t)), "%")
			sys.stdout.write(progressStr)
			sys.stdout.flush()
			sys.stdout.write("\b"*len(progressStr))
			MakeMovieFrame(conf, i, t, corr, r, rad, potData)
			if isTest:
				clf()
				MakeMovieFrame(conf, len(t)-1, t, corr, r, rad, potData)
				show()
				return
			savefig("movie/frame%05i.png" % i)

		close()

	if generateMovie:
		mymovie = myplot.MakeMovie()
		conf.Movie.Apply(mymovie)
		mymovie.CreateMovie()

	rcParams["interactive"] = interactive

def LoadDifferenceBasis(**args):
	boundE, boundV = LoadBoundEigenstates(**args)
	V = zeros((2*(len(boundE)/2), boundV.shape[1]), dtype=boundV.dtype)
	E = zeros((2*(len(boundE)/2)), dtype=boundE.dtype)
	for i in range(len(boundE)/2):
		V[2*i,:] = (boundV[2*i,:] - boundV[2*i+1,:]) / sqrt(2)
		V[2*i+1,:] = (boundV[2*i,:] + boundV[2*i+1,:]) / sqrt(2)
		E[2*i] = (boundE[2*i] + boundE[2*i+1]) / 2.0
		E[2*i+1] = (boundE[2*i] + boundE[2*i+1]) / 2.0

	return E, V

		


def MakeOddEvenMovie(**args):
	args['config'] = "config.ini"
	args['imtime'] = False
	inputfile = GetInputFile(**args)

	conf = SetupConfig(**args)
	prop = pyprop.Problem(conf)
	prop.SetupStep()

	LoadInitialState(prop, **args)
	initPsi = prop.psi.Copy()

	E, V = LoadDifferenceBasis(**args)

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	dr = diff(r)[0]

	outputCount = 500
	if "outputCount" in args:
		outputCount = args["outputCount"]

	curIndex = 0

	interactive = rcParams["interactive"]
	rcParams["interactive"] = False

	figure(figsize=(8,8))

	def output():
		corr = abs(dot(V, prop.psi.GetData()[:,0]))**2
		curtime = t/femtosec_to_au
		
		progressStr = "%#3i%s Complete" % (curIndex*100 / float(outputCount), "%")
		sys.stdout.write(progressStr)
		sys.stdout.flush()
		sys.stdout.write("\b"*len(progressStr))

		clf()
		CorrelationBarPlot(abs(corr)**2)
		axis([-1,20,0,0.35])
		title("t = %3.3f" % curtime)
		savefig("movie/frame%05i.png" % curIndex)

	
	#output initial state
	for t in prop.Advance(outputCount):
		output()
		curIndex += 1

	output()

	mymovie = myplot.MakeMovie()
	conf.Movie.Apply(mymovie)
	mymovie.CreateMovie()

	rcParams["interactive"] = interactive


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

	boundE, boundV = LoadRotatedBoundEigenstates(**args)
	contE1, contV1, contE2, contV2 = LoadContinuumEigenstates(**args)
	dE = average(diff(contE2))
	E = r_[-0.5:0:dE]

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	dr = diff(r)[0]

	initProj = dot(boundV, prop.psi.GetData()[:,0])
	initPhases = numpy.arctan2(initProj.imag, initProj.real)

	if "initStates" in args:
		initStates = args["initStates"]
		initPsi.GetData()[:] = 0

		initStatesWeight = args.get("initStatesWeight", "one")
		if initStatesWeight == "one":
			#use equal population in the selected states
			for state in initStates:
				initPsi.GetData()[:,0] += boundV[state,:]
				print "Input Energy = ", boundE[state]
	
				initPsi.Normalize()
		elif initStatesWeight == "fc":
			#use franck-condon population for the selected staes
			for state in initStates:
				initPsi.GetData()[:,0] += abs(initProj[state]) * boundV[state,:] / dr
				print "Input Energy = ", boundE[state]

		else:
			raise Exception("Invalid initStatesWeight '%s'" % initStatesWeight)
	

		prop.psi.GetData()[:] = initPsi.GetData()
			

	if "testEnergy" in args:
		initStates = args["initStates"]
		prop.psi.GetData()[:] = 0
		for state in initStates:
			prop.psi.GetData()[:,0] = boundV[state,:]
			print prop.GetEnergyExpectationValue() - boundE[state]

		return

	relativePhase = "initial"
	if "relativePhase" in args:
		relativePhase = args["relativePhase"]
		print "Using relative phase '%s'" % relativePhase

	timeList = []
	corrList = []
	radialData = []
	expectRList = []


	outputCount = 500
	if "outputCount" in args:
		outputCount = args["outputCount"]

	def output():
		"""
		Performs one output from the propagation. We have many 
		advance loops, and creating an inner function like this is a 
		nice way of grouping the output functionality
		"""
		if relativePhase == "initial":
			curPhases = initPhases
		elif relativePhase == "follow":
			curPhases = initPhases - boundE * t
		else:
			raise Exception("Unknown relativePhase '%s'" % relativePhase)
		corrList.append(dot(boundV, prop.psi.GetData()[:,0]) * exp(-1.0j * curPhases))
		radialData.append(prop.psi.GetData().copy())
		expectR = dr*sum( conj(prop.psi.GetData()[:,0]) * r * prop.psi.GetData()[:,0] ) 
		#only calc <R> for the bound potential
		#dr*sum( conj(prop.psi.GetData()[:,1]) * r * prop.psi.GetData()[:,1] )
		expectRList.append(expectR)
		
		timeList.append(t/femtosec_to_au)


	""" Don't use pump pulse

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
	"""

	if "removeStates" in args:
		if args["removeStates"] == "odd":
			states = r_[1:len(boundE):2]
		elif args["removeStates"] == "even":
			states = r_[0:len(boundE):2]
		elif iterable(args["removeStates"]):
			states = array(args["removeStates"], dtype=int)
		else:
			raise Exception("Invalid removeStates %s" % args["removeStates"])
		
		data = prop.psi.GetData()[:,0].copy()
		for i in states:
			print "Removing %i" % i
			data -= dot(conj(boundV[i, :]), data) * boundV[i,:]/dr
		prop.psi.GetData()[:,0] = data

		#plot(r, abs(prop.psi.GetData()[:,0].copy())**2)

	if "fastForward" in args:
		fastForward = args["fastForward"]
		proj = dot(boundV, prop.psi.GetData()[:,0])
		proj *= exp(-1.0j*fastForward*boundE)
		prop.psi.GetData()[:] = 0
		prop.psi.GetData()[:,0] = dot(conj(boundV.transpose()), proj) / dr
		prop.PropagatedTime = fastForward

	
	#prop.Duration = fullDuration
	tPrev = 0
	t = prop.PropagatedTime
	output()
	for t in prop.Advance(outputCount):
		output()
		if t-tPrev > 20*femtosec_to_au:
			norm = prop.psi.GetNorm()
			corr = sum(abs(dot(boundV, prop.psi.GetData()[:,0]))**2)
			print "t = %f, N = %f, Corr = %f" % (t/femtosec_to_au, norm, corr) 
			tPrev = t

	t = prop.PropagatedTime
	output()

	norm = prop.psi.GetNorm()
	corr = sum(abs(dot(boundV, prop.psi.GetData()[:,0]))**2)
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
