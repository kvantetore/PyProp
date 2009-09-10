def SubmitTPDIRun():
	"""
	Calculate two-photon double ionization
	"""
	configFile = "config_stabilization_freq_5.ini"
	timestep = 0.01
	#frequency = 1.7
	cycles = 20
	#pulse_duration = cycles * 2 * pi / frequency
	#pulse_duration = 2 * femtosec_to_au
	#duration = 1.5 * pulse_duration
	lmax = 5
	L = [0,1,2,3]
	#freqList = [1.9, 1.94, 1.98, 1.99]
	freqList = linspace(1.47, 1.99, 15)
	#freqList = [1.99]
	amplitude = field_from_intensity(1e13)
	statePostfix = "1s1s_1"

#	radialGrid = {\
#		'xsize' : 20, \
#		'xmax' : 150, \
#		'order' : 5, \
#		'bpstype' : 'exponentiallinear', \
#		'xpartition' : 6.28344983487, \
#		'gamma' : 2.30042679472}

#	radialGrid = {\
#		'xsize' : 22, \
#		'xmax' : 200, \
#		'order' : 5, \
#		'bpstype' : 'exponentiallinear', \
#		'xpartition' : 8.63130809996, \
#		'gamma' : 2.19549141339}

	radialGrid = {\
		'xsize' : 20, \
		'xmax' : 400, \
		'order' : 5, \
		'bpstype' : 'exponentiallinear', \
		'xpartition' : 10.0, \
		'gamma' : 3.0}


	args = {
		'account':'fysisk', \
		'procPerNode':1, \
		'procMemory':"4000M", \
		'config':configFile, \
		'storeInitialState':True, \
		'saveWavefunctionDuringPropagation':True, \
		'writeScript':False, \
		'useStoredPotentials': False, \
		'radialGrid':radialGrid, \
		#'multipoleCutoff = 5, \
		'lmax':lmax, \
		'L': L, \
		}

	radialPostfix = "_".join(GetRadialGridPostfix(**args))
	angularPostFix = "_".join(GetAngularGridPostfix(**args))

	#Set initial state to ground state or excited states
	if statePostfix == "1s1s_1":
		args["initialStateL"] = 0
		args["initialStateIndex"] = 0
	elif statePostfix == "1s2s_3":	
		args["initialStateL"] = 0
		args["initialStateIndex"] = 1
	elif statePostfix == "1s2s_1":	
		args["initialStateL"] = 0
		args["initialStateIndex"] = 2
	elif statePostfix == "1s2p_3":	
		args["initialStateL"] = 1
		args["initialStateIndex"] = 0
	elif statePostfix == "1s2p_1":	
		args["initialStateL"] = 1
		args["initialStateIndex"] = 1
	else:
		raise Exception("Unknown statePostfix %s" % (statePostfix))

	
	#outputDir = "tpdi/%s_%s_%s/" % (radialPostfix, angularPostFix, statePostfix)
	#if not os.path.exists(outputDir):
	#	print "Created output dir: %s" % outputDir
	#	os.makedirs(outputDir)

	#Check that boundstate files exists
	boundstateArgs = dict(args)
	boundstateArgs["L"] = args["initialStateL"]
	boundstateFile = GetBoundstatesFilename(**boundstateArgs)
	if not os.path.exists(boundstateFile):
		raise Exception("Required boundstate File %s does not exist. Run SetupBoundstates first." % boundstateFile)

	depend = None
	#if not os.path.exists(GetInitialStateFilename(**args)):
	#	dependJob = RunSubmitFullProcCount(RunStabilizationInputFile, walltime=timedelta(hours=0.5), **args)
	#	depend = "afterany:%s" % (dependJob,)

	for frequency in freqList:
		pulse_duration = cycles * 2 * pi / frequency
		duration = pulse_duration + 2 * femtosec_to_au
		intensityStr = GetIntensityStringFromAmplitude(amplitude)

		outputDir = "tpdi/freq_%s_%s_%s_%s/" % (frequency, radialPostfix, angularPostFix, statePostfix)
		if not os.path.exists(outputDir):
			print "Created output dir: %s" % outputDir
			os.makedirs(outputDir)

		name = outputDir + "tpdi_I_%s_freq_%s_dt_%1.e_T_%3.1f" % (intensityStr, frequency, timestep, duration)
		args['walltime'] = timedelta(hours=100)
		args['timestep'] = timestep
		args['duration'] = duration
		args['frequency'] = frequency
		args['amplitude'] = amplitude/frequency
		args['pulse_duration'] = pulse_duration
		args['outputCount'] = 100
		args['outputFilename'] = name
		args['findGroundstate'] = False
		args['saveWavefunctionDuringPropagation'] = True
		args['laserOff'] = False
		args['depend'] = depend
		RunSubmitFullProcCount(RunStabilization, **args)


def SubmitTPDIEigenvaluesJob():
	"""
	Calculate a number of the lowest eigenvalues of Helium for a range of L's
	"""
	lmax = 5
	radialGrid = {\
		'xsize' : 25, \
		'xmax' : 400, \
		'order' : 5, \
		'bpstype' : 'exponentiallinear', \
		'xpartition' : 10., \
		'gamma' : 3.}

	shift = [ 
		-2.903,					# L == 0 
		-2.3225840129999999,    # L == 1
		-2.05563627,            # L == 2
		-2.0306810854999999,    # L == 3
		-2.0133175794999998,    # L == 4
		-2.0096329019999999,    # L == 5
		]

	#for L in range(lmax+1):
	#for L in [0,1,2,3]:
	for L in [0]:
		RunSubmitFullProcCount(RunFindBoundstates, \
			account='fysisk',\
			procPerNode=1, \
			procMemory="4000M", \
			walltime=timedelta(hours=150, minutes=00), \
			config="config_tpdi_eigenvalues.ini", \
			writeScript=False, \
			radialGrid=radialGrid, \
			L = L, \
			lmax = lmax, \
			shift = shift[L], \
			eigenvalueCount = 10, \
			eigenvalueBasisSize = 60)
