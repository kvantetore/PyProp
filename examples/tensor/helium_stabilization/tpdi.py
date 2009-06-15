def SubmitTPDIRun():
	"""
	Calculate two-photon double ionization
	"""
	configFile = "config_stabilization_freq_5.ini"
	timestep = 0.01
	frequency = 1.65
	#cycles = 6
	#pulse_duration = cycles * 2 * pi / frequency
	pulse_duration = 2 * femtosec_to_au
	duration = 1.5 * pulse_duration
	lmax = 5
	L = [0,1,2,3]
	freqList = [1.65]
	amplitude = field_from_intensity(1e12)

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
		'xsize' : 24, \
		'xmax' : 400, \
		'order' : 5, \
		'bpstype' : 'exponentiallinear', \
		'xpartition' : 10.0, \
		'gamma' : 3.0}


	args = {
		'account':'fysisk', \
		'procPerNode':1, \
		'procMemory':"8000M", \
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
	#statePostfix = "groundstate"
	statePostfix = "1s1s"
	args["shift"] = -2.9
	args["eigenstateIndex"] = 0
	args["groundstateFilename"] = "output/initialstates/helium_%s_%s_%s.h5" % (statePostfix, radialPostfix, angularPostFix)
	
	outputDir = "tpdi/freq_%s_%s_%s_%s/" % (frequency, radialPostfix, angularPostFix, statePostfix)
	if not os.path.exists(outputDir):
		print "Created output dir: %s" % outputDir
		os.makedirs(outputDir)

	depend = None
	if not os.path.exists(GetInitialStateFilename(**args)):
		dependJob = RunSubmitFullProcCount(RunStabilizationInputFile, walltime=timedelta(hours=0.5), **args)
		depend = "afterany:%s" % (dependJob,)

	for frequency in freqList:
		name = outputDir + "tpdi_I_%i_freq_%s_dt_%1.e_T_%3.1f" % (amplitude, frequency, timestep, duration)
		args['walltime'] = timedelta(hours=50)
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

