import signal

#----------------------------------------------------------------------------------------------------
# Problem
#----------------------------------------------------------------------------------------------------
class Problem:
	"""
	This is the main class of pyprop. It reads a config object and sets up everything
	to allow propagation. See the examples folder for examples how to use this class.
	"""

	def __init__(self, config):
		self.Config = config
		try:
			#Enable redirect
			if hasattr(config.Propagation, "silent"):
				self.Silent = config.Propagation.silent
			else:
				self.Silent = False

			Redirect.Enable(self.Silent)
			
			print "Creating DistributionModel..."
			self.Distribution = CreateDistribution(config)
		
			print "Creating Representation..."
			self.Representation = CreateRepresentation(config, self.Distribution)
		
			print "Creating Wavefunction..."
			self.psi = CreateWavefunction(config, self.Representation)
		
			print "Creating Momentum Evaluator..."
			self.Propagator = CreatePropagator(config, self.psi)
		
			#apply propagation config
			config.Propagation.Apply(self)

			#Disable redirect
			Redirect.Disable()

		except:
			#Diasable redirect
			Redirect.Disable()
			raise
			

	def GetGrid(self):
		"""
		Helper function to construct the grid in each rank of the wavefunction

		This method is most useful for cartesian coordinates, and possibly other
		non-compressed grids. For compressed grids (like spherical) the grid returned
		here will typically just be the idices converted to double
		"""
		#set up grid
		repr = self.psi.GetRepresentation()
		grid = [repr.GetLocalGrid(i) for i in range(0, self.psi.GetRank())]
		return grid
		
		
	#Propagation-------------------------------------------------
	def SetupStep(self):
		"""
		Runs the nescescary setup routines to allow propagation.

		This function must be called before the first call to AdvanceStep()
		or Advance().
		"""
		try:
			#Enable redirect
			Redirect.Enable(self.Silent)
			
			print "Starting setup timestep..."
			print "    Setting up Propagator."
			if self.Propagator != None:
				self.Propagator.SetupStep(self.TimeStep )

			print "    Setting up initial wavefunction"
			self.SetupWavefunction()
			
			self.PropagatedTime = 0
			print "Setup timestep complete."
	
			#Disable redirect
			Redirect.Disable()

		except:
			#Diasable redirect
			Redirect.Disable()
			raise	

	def AdvanceStep(self):
		"""
		Advances the wavefunction one timestep.
	
		most of the work is done in the propagator object. This function
		is merley to keep a track of propagated time, and provide a simple 
		interface to the user.
		"""
		self.Propagator.AdvanceStep(self.PropagatedTime, self.TimeStep )
		if self.Propagator.RenormalizeActive:
			self.psi.Normalize()

		self.PropagatedTime += abs(self.TimeStep)

	def Advance(self, yieldCount, duration=0):
		"""
		Returns a generator for advancing the wavefunction a number of timesteps. 
		If duration is specified the wavefunction is propagated until propagated 
		time >= duration.

		if yieldCount is a number, it is the number of times this routine will 
		yield control back to the caller. 
		If is a boolean, yieldCount=True will make this function yield every timestep
		and yieldCount=False will make it yield one time. (at the last timestep (i think))
		"""
		if duration == 0:
			duration = self.Duration

		#Determine how often we should yield.
		if yieldCount.__class__ is bool:
			yieldStep = 1
		else:
			yieldStep = int((duration / abs(self.TimeStep)) / yieldCount)

		#Modify the interrupt signal handler
		try:
			InterruptHandler.UnRegister()
		except: pass
		InterruptHandler.Register()
	
		index = 0
		while self.PropagatedTime < duration:
			#next timestep	
			self.AdvanceStep()

			#check keyboard interrupt
			if InterruptHandler.IsInterrupted():
				InterruptHandler.ProcessInterrupt()
			
			index += 1
			if index % yieldStep == 0:
				yield self.PropagatedTime
			
		InterruptHandler.UnRegister()
				
	def GetEnergy(self):
		"""
		Advance one timestep without normalization and imaginary time.
		we can then find the energy by looking at the norm of the decaying
		wavefunction.

		This will only work when negaive imaginary time is used to propagate
		to the groundstate, and will only give a good estimate for the ground
		state energy when the wavefunction is well converged to the ground state.

		Because it uses the norm of the wavefunction to measure energy, it will 
		be very sensitive to differences in the shape of the wavefunction.
		"""
		
		if isreal(self.TimeStep):
			raise "Can only find energy for imaginary time propagation"
		
		renorm = self.Propagator.RenormalizeActive
		self.psi.Normalize()
		self.Propagator.RenormalizeActive = False
		self.AdvanceStep()
		self.Propagator.RenormalizeActive = renorm
		
		norm = self.psi.GetNorm()
		energy = - log(norm) / (2 * abs(self.TimeStep))
		return energy
	
	#Initialization-----------------------------------------------
	def ApplyConfigSection(self, configSection):
			self.TimeStep = configSection.timestep
			self.Duration = configSection.duration

	def SetupWavefunction(self):
		"""
		Initializes the wavefunction according to the initially specified configuration
		file.
		The supported initial condition types are:
		InitialConditionType.Function
			see SetupWavefunctionFunction()
			
		InitialConditionType.File
			See SetupWavefunctionFile()

		"""
		type = self.Config.InitialCondition.type
		if type == InitialConditionType.Function:
			self.SetupWavefunctionFunction(self.Config)
		elif type == InitialConditionType.File:
			self.SetupWavefunctionFile(self.Config)
		elif type == InitialConditionType.Class:
			self.SetupWavefunctionClass(self.Config, self.psi)
		else:
			raise "Invalid InitialConditionType: " + config.InitialCondition.type
			
	def SetupWavefunctionClass(self, config, psi):
		classname = config.InitialCondition.classname
		
	  	#Create globals
		glob = dict(sys.modules['__main__'].__dict__)
		glob.update(globals())	
	
		#try to
		evaluator = None
		try: evaluator = eval(classname + "()", glob)
		except: pass

		if evaluator == None:
			try: evaluator = eval(classname + "_" + str(psi.GetRank()) + "()", glob)
			except: pass

		if evaluator == None:
			raise "Invalid classname", classname

		config.Apply(evaluator)
		config.InitialCondition.Apply(evaluator)
		evaluator.SetupWavefunction(psi)

	def SetupWavefunctionFunction(self, config):
		"""
		Initializes the wavefunction from a function specified in the InitialCondition
		section of config. The function refrence config.InitialCondition.function 
		is evaulated in all grid points of the wavefunction, and the wavefunction is
		set up accordingly.

		for example. if this is a rank 2 problem, then this will initialize the wavefunction
		to a 2D gaussian wavepacket.
		
		def func(x,  conf):
			return exp(- x[0]**2 + x[1]**2])
		
		config.InitialCondition.function = func
		prop.SetupWavefunctionFunction(config)

		REMARK: This function should probably not be called on directly. Use SetupWavefunction()
		instead. That function will automatically determine the type of initial condition to be used.		
		"""
		func = config.InitialCondition.function
		conf = config.InitialCondition
		
		#TODO: We should get this class from Propagator in order to use a different
		#evaluator for a compressed (i.e. spherical) grid
		evalfunc = eval("core.SetWavefunctionFromGridFunction_" + str(self.psi.GetRank()))
		evalfunc(self.psi, func, conf)
	
	def SetupWavefunctionFile(self, config):
		"""
		Initializes the wavefunction from a file specified in the InitialCondition.filename
		according to the format InitialCondition.format. If InitialCondition.format is not 
		specified, this routine should automatically try to figure out which format it it, 
		but this is not yet implemented.

		See the functions LoadWavefunction*() SaveWavefunction*() for details on how to
		load and save wavefunctions

		REMARK: This function should probably not be called on directly. Use SetupWavefunction()
		instead. That function will automatically determine the type of initial condition to be used.		

		REMARK2: Currently it is not possible to interpolate a wavefunction between different grids. 
		this means that exactly the same grid used for saving the wavefunction must be used when loading
		it
		"""
		format   = config.InitialCondition.format
		filename = config.InitialCondition.filename
		
		if format == WavefunctionFileFormat.Ascii:
			self.LoadWavefunctionAscii(filename)
		elif format == WavefunctionFileFormat.Binary:
			self.LoadWavefunctionPickle(filename)
		elif format == WavefunctionFileFormat.HDF:
			datasetPath = str(config.InitialCondition.dataset)
			self.LoadWavefunctionHDF(filename, datasetPath)
		else:
			raise "Invalid file format: " + format
	
	#(de)serialization---------------------------------------------
	def LoadWavefunctionData(self, newdata):
		data = self.psi.GetData()
		if newdata.shape != data.shape:
			raise "Invalid shape on loaded wavefunction, got " + str(newdata.shape) + " expected " + str(data.shape)
		data[:] = newdata
		
	def SaveWavefunctionPickle(self, filename):
		serialization.SavePickleArray(filename, self.psi.GetData())
	
	def LoadWavefunctionPickle(self, filename):
		arr = serialization.LoadPickleArray(filename)
		self.LoadWavefunctionData(arr)

	def LoadWavefunctionHDF(self, filename, datasetPath):
		serialization.LoadWavefunctionHDF(filename, datasetPath, self)

	def SaveWavefunctionHDF(self, filename, datasetPath):
		serialization.SaveWavefunctionHDF(filename, datasetPath, self)

	def SaveWavefunctionAscii(self, filename):
		pylab.save(filename, self.psi.GetData())
	
	def LoadWavefunctionAscii(self, filename):
		arr = pylab.load(filename)
		self.LoadWavefunctionData(arr)
		
	

