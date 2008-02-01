import commands
import pyprop.utilities.submitpbs as submit

def PlotDelayScan(**args):
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
					projList.append([node.eigenstateProjection[-1,:]])
					timeList.append(float(name[len("delay_"):]))
		finally:
			f.close()

	time = array(timeList)
	proj= array(projList)

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

		executable = './run_delay_scan.py' + commands.mkarg(repr(args))  
		plist.append(os.popen(executable))

	for p in plist:
		print p.readlines()
		print p.close()
	

def RunDelayScan(**args):
	#Required parameters
	delayList = args["delayList"]
	outputfile = args["outputfile"]
	molecule = args["molecule"]

	print "USING OUTPUTFILE ", outputfile

	for delay in delayList:
		args["pulseDelay"] = delay*femtosec_to_au
		args["outputpath"] = "/delay_%i" % (delay)
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
			if index % 20 == 0:
				psiTimeList.append(t)
				psiList.append(prop.psi.GetData().reshape((1,) + prop.psi.GetData().shape))
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

