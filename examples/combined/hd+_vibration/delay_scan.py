import commands
import pyprop.utilities.submitpbs as submit

def GetScanDelayEnergyDistribution(**args):
	molecule = args["molecule"]
	outputfile = args["outputfile"]
	partitionCount = args["partitionCount"]
	threshold = -0.5
	inputfile = GetInputFile(**args)

	f = tables.openFile(inputfile, "r")
	try:
		E1 = f.root.complete.binding.eigenvalues[:]
		V1 = f.root.complete.binding.eigenstates[:]
		E2 = f.root.complete.unbinding.eigenvalues[:]
		V2 = f.root.complete.unbinding.eigenstates[:]
	finally:
		f.close()

	#We only project on certain states
	maxIndex1 = max(where(E1<0)[0]) + 1
	thresholdIndex = max(where(E1<=threshold)[0])
	E1 = E1[thresholdIndex:maxIndex1]
	V1 = conj(V1[:, thresholdIndex:maxIndex1].transpose())

	maxIndex2 = max(where(E2<0)[0]) + 1
	E2 = E2[:maxIndex2]
	V2 = conj(V2[:, :maxIndex2].transpose())

	#To add up the energies we must interpolate the values from 
	#the E1 and E2 grids to a common E-grid
	dE = average(diff(E2))
	E = r_[-0.5:0:dE]

	print "Using dE             = %.4f" % dE
	print "Using thresholdIndex = %i" % thresholdIndex


	projList = []
	timeList = []

	for i in range(partitionCount):
		filename = outputfile % i

		f = tables.openFile(filename, "r")
		try:
			for node in f.listNodes(f.root):
				name = node._v_name
				if name.startswith("delay_"):
					psi = None
					try:
						psi = node.wavefunction[-1,:,:]
						t = float(name[len("delay_"):])
					except:
						print "Could not process %s" % name
					if psi != None:
						proj = CalculateEnergyDistribution(psi, E, E1, V1, E2, V2) 
	 					projList.append(proj)
						timeList.append(t)

		finally:
			f.close()

	sortedIndex = argsort(timeList)
	time = array(timeList)[sortedIndex]
	proj= array(projList)[sortedIndex,:]

	return time, E, proj



def CalculateEnergyDistribution(psi, outputEnergies, E1, V1, E2, V2):
	proj1 = dot(V1, psi[:,0])
	proj2 = dot(V2, psi[:,1])

	interp1 = spline.Interpolator(E1, abs(proj1)**2)
	interp2 = spline.Interpolator(E2, abs(proj2)**2)

	return array([interp1.Evaluate(e) + interp2.Evaluate(e) for e in outputEnergies])

def MakeDelayScanPlots():
	interactive = rcParams["interactive"]
	rcParams["interactive"] = False

	def Plot(filename):
		t, c = GetScanDelayCorrelation(outputfile=filename, partitionCount=1)
		figure()
		ionization = 1 - sum(abs(c[:,:20])**2, axis=1)
		plot(t, ionization, label="Ionization")
		for i in range(10):
			if i < 5:
				label = "$n = %i$" % i
			else:
				label = "_nolegend_"
			plot(t, abs(c[:,i])**2, label=label)
		axis([0,50, 0,  0.8])
		legend()
		xlabel("Pulse Delay")
		ylabel("Proabability")

	Plot("outputfiles/hd+/delay_%i.h5")
	title("$HD^+$")
	
	Plot("outputfiles/h2+/delay_%i.h5")
	title("$H_2^+$")
	
	Plot("outputfiles/d2+/delay_%i.h5")
	title("$D_2^+$")
	
	Plot("outputfiles/hd+/nodipole_delay_%i.h5")
	title("$HD^+$ without static dipole moment")

	#Plot the difference between HD+ with and without the static dipole term
	figure()
	t, c = GetScanDelayCorrelation(outputfile="outputfiles/hd+/delay_%i.h5", partitionCount=1)
	t2, c2 = GetScanDelayCorrelation(outputfile="outputfiles/hd+/nodipole_delay_%i.h5", partitionCount=1)
	plot(t, abs(c)**2 - abs(c2)**2)
	title("Difference between HD+ with and without static dipole moment")
	xlabel("Pulse Delay")
	ylabel("Proabability difference")


	#restore previous interactive setting
	rcParams["interactive"] = interactive
	show()


def RepackDelayScan(**args):
	outputfile = args["outputfile"]
	partitionCount = args["partitionCount"]
	repackFile = args["repackFile"]

	print "Repacking correlation"
	time, proj = GetScanDelayCorrelation(**args)
	print "Repacking energy distribution"
	t, E, corr = GetScanDelayEnergyDistribution(**args)
	print "Repacking norm"
	t, norm = GetScanDelayNorm(**args)

	

	output = tables.openFile(repackFile, "w")
	try:
		SaveArray(output, "/", "pulse_delay", time)
		SaveArray(output, "/", "energy", E)

		SaveArray(output, "/", "final_correlation", abs(proj)**2)
		SaveArray(output, "/", "energy_distribution", corr)
		SaveArray(output, "/", "norm", norm)

	finally:
		output.close()

	
def GetScanDelayCorrelation(**args):
	outputfile = args["outputfile"]
	partitionCount = args["partitionCount"]

	projList = []
	timeList = []

	for i in range(partitionCount):
		filename = outputfile % i

		f = tables.openFile(filename, "r")
		try:
			for node in f.listNodes(f.root):
				name = node._v_name
				if name.startswith("delay_"):
					try:
						proj = node.eigenstateProjection[-1,:]
						t = float(name[len("delay_"):])
						projList.append(proj)
						timeList.append(t)
					except:
						print "Could not process %s" % name
		finally:
			f.close()

	sortedIndex = argsort(timeList)
	time = array(timeList)[sortedIndex]
	proj= array(projList)[sortedIndex]

	return time, proj

def GetScanDelayNorm(**args):
	outputfile = args["outputfile"]
	partitionCount = args["partitionCount"]

	projList = []
	timeList = []

	for i in range(partitionCount):
		filename = outputfile % i

		f = tables.openFile(filename, "r")
		try:
			for node in f.listNodes(f.root):
				name = node._v_name
				if name.startswith("delay_"):
					try:
						proj = node.norm[-1]
						t = float(name[len("delay_"):])
						projList.append(proj)
						timeList.append(t)
					except:
						print "Could not process %s" % name
		finally:
			f.close()

	sortedIndex = argsort(timeList)
	time = array(timeList)[sortedIndex]
	proj= array(projList)[sortedIndex]

	return time, proj

def SubmitDelayScan(**args):
	delayList = args["delayList"]
	outputfile = args["outputfile"]
	molecule = args["molecule"]
	partitionCount = args["partitionCount"]

	partitionSize = int(ceil(len(delayList)/float(partitionCount)))

	plist = []
	
	for i in range(partitionCount):
		args["delayList"] = delayList[i*partitionSize:(i+1)*partitionSize]
		args["outputfile"] = outputfile % i

		executable = './run_delay_scan.py %s > /dev/null' % commands.mkarg(repr(args))  
		plist.append(os.popen(executable))

	while len(plist) > 0:
		for p in plist:
			try: 
				print p.next()
			except StopIteration:
				print p.close()
				plist.remove(p)
				break

	

def RunDelayScan(**args):
	#Required parameters
	delayList = args["delayList"]
	outputfile = args["outputfile"]
	molecule = args["molecule"]

	print "USING OUTPUTFILE ", outputfile

	for delay in delayList:
		args["pulseDelay"] = delay*femtosec_to_au
		args["outputpath"] = "/delay_%i" % (delay)
		print "Propagating %s with pulse delay %ifs" % (molecule, delay)
		PropagateDelayScan(**args)

def PropagateDelayScan(**args):
	"""
	Propagates one run for the scan pulse delay experiment, and
	saves the result to disk.

	The following arguments are required
	molecule: The molecule to simulate "hd+" or "d2+" or "h2+"
	outputfile: The hdf5 file to save the results to (i.e. pulsedelay_30.h5)
	outputpath: The group inside outputfile where (i.e. "/delay_30fs")

	"""	
	args['config'] = "config.ini"
	args['imtime'] = False
	inputfile = GetInputFile(**args)
	outputfile = args["outputfile"]
	outputpath = args["outputpath"]

	conf = SetupConfig(**args)
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

	if len(eigenstates.shape) != 2:
		raise Exception("Please implement for Rank!=2")
	eigenstateSize = len(eigenstates[0,:])

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	timeList = []
	initCorrList = []
	corrList = []
	normList = []
	psiTimeList = []
	timeList = []

	pulseStart = conf.ElectronicCoupling.delay
	pulseDuration = conf.ElectronicCoupling.duration
	
	output = tables.openFile(outputfile, "a")
	try:
		RemoveNodeIfExists(output, outputpath, "wavefunction")
		atom = tables.ComplexAtom(16)
		shape = (0,) + prop.psi.GetData().shape
		psiList = output.createEArray(outputpath, "wavefunction", atom, shape, createparents=True)
		
		#0) Store the initial packet
		timeList.append(prop.PropagatedTime)
		corrList.append(GetEigenstateCorrelations(prop, eigenstates))
		initCorrList.append(prop.psi.InnerProduct(initPsi))
		normList.append(prop.psi.GetNorm())
		psiTimeList.append(prop.PropagatedTime)
		psiList.append(prop.psi.GetData().reshape((1,) + prop.psi.GetData().shape))
		
		#1) Propagate until the pulse starts
		prop.Duration = pulseStart - 2*pulseDuration
		for t in prop.Advance(minimum(20, prop.Duration)): 
			timeList.append(t)
			corrList.append(GetEigenstateCorrelations(prop, eigenstates))
			initCorrList.append(prop.psi.InnerProduct(initPsi))
			normList.append(prop.psi.GetNorm())
			psiTimeList.append(t)
			psiList.append(prop.psi.GetData().reshape((1,) + prop.psi.GetData().shape))
		
		#1a) Save the wavepacket to see how much has flowed out
		timeList.append(prop.PropagatedTime)
		corrList.append(GetEigenstateCorrelations(prop, eigenstates))
		initCorrList.append(prop.psi.InnerProduct(initPsi))
		normList.append(prop.psi.GetNorm())
		psiTimeList.append(prop.PropagatedTime)
		psiList.append(prop.psi.GetData().reshape((1,) + prop.psi.GetData().shape))
		
		#2) Propagate until the end of the pulse
		prop.Duration = pulseStart + 2*pulseDuration
		index = 0
		for t in prop.Advance(True):
			timeList.append(t)
			corrList.append(GetEigenstateCorrelations(prop, eigenstates))
			initCorrList.append(prop.psi.InnerProduct(initPsi))
			normList.append(prop.psi.GetNorm())
			#if index % 20 == 0:
			#	psiTimeList.append(t)
			#	psiList.append(prop.psi.GetData().reshape((1,) + prop.psi.GetData().shape))
			index+=1
		
		
		#2a) Save the wavepacket at the end 
		timeList.append(prop.PropagatedTime)
		corrList.append(GetEigenstateCorrelations(prop, eigenstates))
		initCorrList.append(prop.psi.InnerProduct(initPsi))
		normList.append(prop.psi.GetNorm())
		psiTimeList.append(prop.PropagatedTime)
		psiList.append(prop.psi.GetData().reshape((1,) + prop.psi.GetData().shape))
		
		SaveArray(output, outputpath, "time", array(timeList))
		SaveArray(output, outputpath, "eigenstateProjection", array(corrList))
		SaveArray(output, outputpath, "initstateProjection", array(initCorrList))
		SaveArray(output, outputpath, "norm", array(normList))
		SaveArray(output, outputpath, "timeWavefunction", array(corrList))
		psiList.close()

	finally:
		output.close()

	return array(timeList), array(corrList)

