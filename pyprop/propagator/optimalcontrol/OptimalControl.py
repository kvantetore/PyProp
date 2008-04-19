
class OptimalControl:
	"""
	"""
	
	def __init__(self, prop):
		self.BaseProblem = prop
		self.PsiSize = prop.psi.GetData().size
		self.PotentialList = prop.Propagator.PotentialList
	
		#Get time step and propagation time from config
		self.TimeStep = prop.TimeStep.real
		self.PropagationTime = prop.Duration
		self.TimeGridSize = int(round(self.PropagationTime / self.TimeStep))
		self.TimeGridResolution = self.TimeStep

		self.ApplyConfigSection(prop.Config)
		self.SetupCommon()
		self.Setup()
		self.SetupPenaltyMatrix()


	def SetupCommon(self):
		#Keep the calculated J's (shows convergence)
		self.J = []
		self.Yield = []

		#Two matrices to hold previous backward and forward propagation solutions (for each time step)
		self.ForwardSolution = zeros((self.PsiSize, self.TimeGridSize), dtype = complex)
		self.BackwardSolution = zeros((self.PsiSize, self.TimeGridSize), dtype = complex)

		#The target state
		self.SetupTargetState()

		#Create a copy of the wavefunction to calculate H|psi>
		self.TempPsi = self.BaseProblem.psi.CopyDeep()
		self.TempPsi2 = self.BaseProblem.psi.CopyDeep()

		#Make a list of the control functions
		self.ControlFunctionList = []
		for controlFunc in self.ControlFunctionNamesList:
			for potential in self.PotentialList:
				if potential.Name == controlFunc:
					self.ControlFunctionList += [potential]

		#A vector to store the optimized control functions
		self.ControlVectors = ones((self.NumberOfControls, self.TimeGridSize), dtype = double)
		for a in range(self.NumberOfControls):
			self.ControlVectors[a,:] *= self.ControlFunctionList[a].ConfigSection.strength

	
	def Setup(self):
		raise NotImplementedError


	def Run(self):
		"""
		Run the Krotov algorithm, until desired yield or MaxIterations is reached.
		"""

		#
		# Main control loop
		# Iterate until converged or MaxIterations is reached
		#
		self.currIter = 0
		self.__YieldNotReached = True
		while (self.currIter < self.MaxIterations) & self.__YieldNotReached:
			
			# Forward propagation
			self.ForwardPropagation()

			# Project on target state and check yield.
			targetProjection = self.ComputeTargetProjection()
			currentYield = abs(targetProjection)**2
			self.Yield.append(currentYield)

			#Get value of cost functional
			curJ = self.ComputeCostFunctional(currentYield)
			self.J.append(curJ)
			norm = self.BaseProblem.psi.GetNorm()
			print "Yield = %.8f, J = %.8f, norm = %.15f" % (currentYield, curJ, norm)

			#Desired yield reached? -> Finished
			if( currentYield > self.YieldRequirement ):
				print "Yield reached!"
				self.__YieldNotReached = False
				continue

			# Backward propagation
			self.BackwardPropagation(targetProjection)

			self.currIter += 1

	
	def ComputeNewControlFunctions(self, timeGridIndex, t, direction):
		"""
		Updates the control function based
		"""
		raise NotImplementedError


	def ComputeCostFunctional(self, currentYield):
		"""
		"""
		raise NotImplementedError


	def ComputeTargetProjection(self):
		"""
		Project solution on target state
		"""
		#Note: Argument is complex conjugated
		return self.BaseProblem.psi.InnerProduct(self.TargetState)


	def ForwardPropagation(self):
		"""
		Propagate forward from 0 to T.
		"""
		
		print "Forward propagation, pass ", self.currIter

		self.SetupStep(Direction.Forward)
		self.InitializeControls(Direction.Forward)
		
		#Initial step
		if self.currIter > 0:
			self.ComputeNewControlFunctions(0, 0.0, Direction.Forward)
		#self.ForwardSolution[:, 0] = self.BaseProblem.psi.GetData()[:]

		for idx, t in enumerate(self.BaseProblem.Advance(self.TimeGridSize)):
			self.ForwardSolution[:, idx] = self.BaseProblem.psi.GetData()[:]
			if idx < self.TimeGridSize - 1 and self.currIter > 0:
				self.ComputeNewControlFunctions(idx+1, t, Direction.Forward)

	
	def BackwardPropagation(self, targetProjection):
		"""
		Propagate backwards from T to 0. Terminal condition is Psi(T) = <F|Psi_forward(T)> |F>,
		where |F> is the desired final state.
		"""

		print "Backward propagation, pass ", self.currIter

		self.SetupStep(Direction.Backward)
		self.InitializeControls(Direction.Backward)

		# Set the initial state
		self.BaseProblem.psi.GetData()[:] = targetProjection * self.TargetState.GetData()[:]
		
		self.ComputeNewControlFunctions(self.TimeGridSize - 1, self.PropagationTime, Direction.Backward)

		#Set control and store backward solution for t = T
		#self.BackwardSolution[:, self.TimeGridSize - 1] = self.BaseProblem.psi.GetData()[:]

		T = self.PropagationTime
		for idx, t in enumerate(self.BaseProblem.Advance(self.TimeGridSize, duration = T)):
			self.BackwardSolution[:, self.TimeGridSize - idx - 1] = self.BaseProblem.psi.GetData()[:]
			if idx < self.TimeGridSize - 2:
				self.ComputeNewControlFunctions(self.TimeGridSize - idx - 2, t, Direction.Backward)
			#if idx < self.TimeGridSize - 1:
			#	self.ControlFunction.ConfigSection.e0 = self.ControlVector[self.TimeGridSize - idx - 2]

	
	def SetupStep(self, direction):
		"""
		Set up propagator for forward or backward run
		"""

		if direction == Direction.Forward:
			self.BaseProblem.TimeStep = complex(self.TimeStep)
			self.BaseProblem.StartTime = 0.0
			self.BaseProblem.Duration = self.PropagationTime
			if hasattr(self.BaseProblem.Propagator, "OdeWrapper"):
				self.BaseProblem.Propagator.OdeWrapper.SetStartTime(0)
			self.BaseProblem.SetupStep()

		elif direction == Direction.Backward:
			self.BaseProblem.TimeStep = -complex(self.TimeStep)
			self.BaseProblem.StartTime = self.PropagationTime
			self.BaseProblem.Duration = 0.0
			if hasattr(self.BaseProblem.Propagator, "OdeWrapper"):
				self.BaseProblem.Propagator.OdeWrapper.SetStartTime(self.PropagationTime)
			self.BaseProblem.SetupStep()

	
	def SetupTargetState(self):
		"""
		Set up the desired target state
		"""
		self.TargetState = self.BaseProblem.psi.CopyDeep()
		self.TargetState.Clear()
		
		self.TargetState.GetData()[self.BaseProblem.Config.FinalState.states] \
			= self.BaseProblem.Config.FinalState.population[:]

		self.TargetState.Normalize()

	
	def InitializeControls(self, direction):
		for a in range(self.NumberOfControls):
			self.ControlFunctionList[a].ConfigSection.strength = self.ControlVectors[a, direction]

	
	def SetupPenaltyMatrix(self):
		raise NotImplementedError



class Direction:
	Forward = 0
	Backward = -1
