
class Krotov:
	"""
	This class implements Krotov's method, which is an optimization algorithm for quantum control.
	The user may specify several control functions which are used to guide the system from a
	specified initial to final state. 

	More precisely, given the TDSE:
	
	             d            /                        \ 
				---- Psi = -i | H_0 + sum( u_i(t)X_i ) | Psi
				 dt           \                        / 
	
	for a given H_0 (time independent) and controls (scalars) u_i, Krotov solves

            ______     |            ____
			\	 /	   |    ____   /    \ 
			 \	/	   |    ____   |    |
			  \/	\__/           \____/

	where J =  <F|Psi> - mu \sum( \int_0^T |u_i(t)|**2 dt ).
	"""

	def __init__(self, prop):
		self.BaseProblem = prop
		self.PsiSize = prop.psi.GetRepresentation().VectorSize

		# Shortcuts to the salient control parameters
		self.ControlCutoff = prop.Config.Krotov.control_cutoff
		self.EnergyPenalty = prop.Config.Krotov.energy_penalty
		self.MaxIterations = prop.Config.Krotov.max_iterations
		self.YieldRequirement = prop.Config.Krotov.yield_requirement
		self.TimeGridSize = prop.Config.Krotov.time_grid_size
		self.TimeStep = prop.Config.Propagation.timestep
		self.StoreSteps = int(round(self.BaseProblem.Duration / self.TimeStep / self.TimeGridSize))
		self.Omega = prop.Config.ControlFunction.omega

		#Keep the calculated J's (shows convergence)
		self.J = []

		#A matrix to hold previous backward propagation solutions (for each time step)
		self.BackwardSolution = zeros((self.PsiSize, self.TimeGridSize), dtype = complex)

		#The target state
		self.SetupTargetState()

		#Create a copy of the wavefunction to calculate H|psi>
		self.TempPsi = prop.GetTempPsi()
		self.TempPsi2 = prop.psi.Copy()

		#Find control function (potential)
		for potential in prop.Propagator.PotentialList:
			if(potential.Name == 'ControlFunction'):
				self.ControlFunction = potential

		#A vector to store the optimized control function
		self.ControlVector = zeros((self.TimeGridSize), dtype = float)
		self.ControlVector = ones(self.TimeGridSize, dtype = float) * self.ControlFunction.ConfigSection.e0
		#self.ControlVector = (random.rand(self.TimeGridSize) - 0.5) * 2.0
		#self.ControlVector = self.ControlVector / max(abs(self.ControlVector)) * self.ControlFunction.ConfigSection.e0


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

			# Interlude: project on target state and check yield
			#            If yield is reached, set flag and cycle
			targetProjection = self.ComputeTargetProjection()
			currentYield = abs(targetProjection)**2
			energyTerm = self.EnergyPenalty * sum(numpy.abs(self.ControlVector)**2) * self.TimeStep
			self.J += [currentYield - energyTerm]
			print "Yield = ", currentYield, " J = ", self.J[-1]
			if( currentYield > self.YieldRequirement ):
				print "Yield reached!"
				self.__YieldNotReached = False
				continue

			# Backward propagation
			self.BackwardPropagation(targetProjection)

			self.currIter += 1


	def ComputeCostFunctional(self):
		"""
		"""
		raise NotImplementedError, 'Not implemented yet'


	def ComputeNewControlFunction(self, timeGridIndex):
		"""
		Helper function, updates the control function based on backward propagation solution.
		"""
		self.ControlFunction.ConfigSection.e0 = 1.0
		self.TempPsi.GetData()[:] = 0.0
		self.ControlFunction.MultiplyPotential(self.TempPsi,0,0)
		self.TempPsi2.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
		newControl = conj(self.TempPsi2.InnerProduct(self.TempPsi))
		fieldValue = -imag(newControl) / (2.0 * self.EnergyPenalty)
		
		# Warn if dt * fieldValue is to great
		errorOrder = abs(self.TimeStep**2 * fieldValue)
		if(errorOrder > 1e-6):
			print "Warning: Omitted dexp terms not small (", errorOrder, ")!"
	
		# Check if control amplitude exceeds control cutoff, and store
		#self.ControlVector[timeGridIndex] = min(abs(self.ControlCutoff), abs(fieldValue)) * sign(fieldValue)
		self.ControlVector[timeGridIndex] = fieldValue


	def ComputeTargetProjection(self):
		"""
		"""
		return conj(self.TargetState.InnerProduct(self.BaseProblem.psi))


	def ForwardPropagation(self):
		"""
		Propagate forward from 0 to T.
		"""
		
		print "Forward propagation, pass ", self.currIter

		self.SetupStep(self.TimeStep)
		self.BaseProblem.SetupWavefunction()

		index = 0
		timeGridIndex = 0
		while self.BaseProblem.PropagatedTime < self.BaseProblem.Duration-1:
			
			# Time to do something?
			if (index % self.StoreSteps == 0):
				if not (self.currIter == 0):
					self.ComputeNewControlFunction(timeGridIndex)

				# Update control function
				self.ControlFunction.ConfigSection.e0 = self.ControlVector[timeGridIndex]

				timeGridIndex += 1
				
			self.BaseProblem.AdvanceStep()
			index += 1

			# Check for convergence here


	def BackwardPropagation(self, targetProjection):
		"""
		Propagate backwards from T to 0. Terminal condition is Psi(T) = <F|Psi_forward(T)> |F>,
		where |F> is the desired final state.
		"""

		print "Backward propagation, pass ", self.currIter

		self.SetupStep(-self.TimeStep)
		self.BaseProblem.PropagatedTime = self.BaseProblem.Duration

		# Set the initial state
		self.BaseProblem.psi.GetData()[:] = targetProjection * self.TargetState.GetData()[:]

		index = 0
		timeGridIndex = 0
		while self.BaseProblem.PropagatedTime > 0.0:

			# Use control function from forward propagation
			if (index % self.StoreSteps == 0):
				timeGridIndex += 1
				self.ControlFunction.ConfigSection.e0 = self.ControlVector[self.TimeGridSize - timeGridIndex]

				# Store backward solution for next iteration
				self.BackwardSolution[:, self.TimeGridSize - timeGridIndex] = self.BaseProblem.psi.GetData()[:]

			# Propagate one step
			self.BaseProblem.AdvanceStep()

			index += 1

	def SetupStep(self,timeStep):
		"""
		Helper function; set up propagator for forward or backward run
		"""
		self.BaseProblem.TimeStep = timeStep
		self.BaseProblem.PropagatedTime = 0.0

	
	def SetupTargetState(self):
		"""
		Set up the desired target state
		"""
		self.TargetState = self.BaseProblem.psi.Copy()
		self.TargetState.GetData()[:] = 0.0
		self.TargetState.GetData()[self.BaseProblem.Config.FinalState.states] = self.BaseProblem.Config.FinalState.population
