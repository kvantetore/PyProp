import signal

class InterruptHandlerClass:
	"""
	KeyboardInterrupt handler, which is used in Problem.Advance, in order to not
	to put the problem in an invalid state when using ctrl+c during propagation
	"""
	def __init__(self):
		self.__interrupt = False
		self.__signum = -1
		self.__frame = -1
		self.__origHandler = None

	def Register(self):
		if self.__origHandler != None:
			raise "Already registered"

		self.__origHandler = signal.signal(signal.SIGINT, self.Handler)

	def UnRegister(self):
		if self.__origHandler == None:
			print "Not registered..."
			raise "Not registered"
			
		signal.signal(signal.SIGINT, signal.default_int_handler)
		self.__origHandler = None
		self.__signum = -1
		self.__frame = -1
		self.__interrupt = False
	
	def Handler(self, signum, frame):
		print "Got keyboard interrupt, will terminate next timestep"
		self.__signum = signum
		self.__frame = frame
		self.__interrupt = True

	def ProcessInterrupt(self):
		if not self.IsInterrupted():
			print "HM? Not interrupted"
			return

		hndlr = self.__origHandler
		signum = self.__signum
		frame = self.__frame
		self.UnRegister()
		hndlr(signum, frame)

	def IsInterrupted(self):
		return self.__interrupt
	
InterruptHandler = InterruptHandlerClass()

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

	def GetGrid(self):
		"""
		Helper function to construct the grid in each rank of the wavefunction

		This method is most useful for cartesian coordinates, and possibly other
		non-compressed grids. For compressed grids (like spherical) the grid returned
		here will typically just be the idices converted to double
		"""
		#set up grid
		repr = self.psi.GetRepresentation()
		grid = [repr.GetLocalGrid(self.psi, i) for i in range(0, self.psi.GetRank())]
		return grid
		
		
	#Propagation-------------------------------------------------
	def SetupStep(self):
		"""
		Runs the nescescary setup routines to allow propagation.

		This function must be called before the first call to AdvanceStep()
		or Advance().
		"""
		print "Starting setup timestep..."
		print "    Setting up Propagator."
		if self.Propagator != None:
			self.Propagator.SetupStep(self.TimeStep )

		print "    Setting up initial wavefunction"
		self.SetupWavefunction()
		
		self.PropagatedTime = 0
		print "Setup timestep complete."
		

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
			self.LoadWavefunction(filename)
		else:
			raise "Invalid file format: " + format
	
	def ChangeDistribution(self):
		"""
		TODO: This should probably be updated or removed... first i must figure out how
		to deal with distributed memory.
		"""
		self.psi.GetRepresentation().GetDistributedModel().ChangeRepresentation(self.psi)
		
	#(de)serialization---------------------------------------------
	def LoadWavefunctionData(self, newdata):
		data = self.psi.GetData()
		if newdata.shape != data.shape:
			raise "Invalid shape on loaded wavefunction, got " + str(newdata.shape) + " expected " + str(data.shape)
		data[:] = newdata
		
	def SaveWavefunction(self, filename):
		SavePickleArray(filename, self.psi.GetData())
	
	def LoadWavefunction(self, filename):
		arr = LoadPickleArray(filename)
		self.LoadWavefunctionData(arr)

	def SaveWavefunctionAscii(self, filename):
		SaveMatlabArray(filename, self.psi.GetData())
	
	def LoadWavefunctionAscii(self, filename):
		arr = LoadMatlabArray(filename)
		self.LoadWavefunctionData(arr)
		
	

