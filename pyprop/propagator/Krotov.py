
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

	    grad(J) = 0

	where J =  <F|Psi> - mu \sum( \int_0^T |u_i(t)|**2 dt ).
	"""

	def __init__(self, prop):
		self.BaseProblem = prop
		self.PsiSize = prop.psi.GetRepresentation().VectorSize

		#Get time step and propagation time from config
		self.TimeStep = prop.Config.Propagation.timestep
		self.PropagationTime = prop.Config.Propagation.duration

		# Shortcuts to the salient control parameters
		self.ControlCutoff = prop.Config.Krotov.control_cutoff
		self.EnergyPenalty = prop.Config.Krotov.energy_penalty
		self.MaxIterations = prop.Config.Krotov.max_iterations
		self.YieldRequirement = prop.Config.Krotov.yield_requirement
		self.TimeGridSize = int(self.PropagationTime / self.TimeStep)
		self.TimeGridResolution = self.TimeStep
		if hasattr(prop.Config.Krotov, "time_grid_size"):
			self.TimeGridSize = prop.Config.Krotov.time_grid_size
			self.TimeGridResolution = self.PropagationTime / self.TimeGridSize

		#Set control cutoff from time step
		err = 1e-6
		self.ControlCutoff = err / self.TimeStep**2

		#Check for the debug option
		self.Debug = False
		if hasattr(prop.Config.Krotov, "debug"):
			self.Debug = prop.Config.Krotov.debug
		
		#Keep the calculated J's (shows convergence)
		self.J = []
		self.Yield = []

		#A matrix to hold previous backward propagation solutions (for each time step)
		self.BackwardSolution = zeros((self.PsiSize, self.TimeGridSize), dtype = complex)

		#The target state
		self.SetupTargetState()

		#Create a copy of the wavefunction to calculate H|psi>
		self.TempPsi = prop.psi.CopyDeep()
		self.TempPsi2 = prop.psi.CopyDeep()

		#Find control function (potential)
		for potential in prop.Propagator.PotentialList:
			if(potential.Name == 'ControlFunction'):
				self.ControlFunction = potential

		#A vector to store the optimized control function
		self.ControlVector = ones(self.TimeGridSize, dtype = double) * self.ControlFunction.ConfigSection.e0


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
			energyTerm = self.EnergyPenalty * self.TimeGridResolution \
				* sum(numpy.abs(self.ControlVector)**2)
			curJ = currentYield - energyTerm
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

		#A final forward propagation
		#if self.__YieldNotReached:
		#	self.ForwardPropagation()


	def ComputeNewControlFunction(self, timeGridIndex, t):
		"""
		Updates the control function based on backward propagation solution:

		    E(t) = - 1 / mu * imag(<etta|X|psi>)

		where E is the new control, mu is the energy penalty, |etta> is the
		backward solution and X the control matrix.
		"""
		#Compute new control function. First the control matrix
		#is multiplied on the current forward solution, and then
		#the result is projected onto the backward solution. 
		self.ControlFunction.ConfigSection.e0 = 1.0
		self.TempPsi.Clear()
		self.TempPsi2.Clear()
		self.ControlFunction.MultiplyPotential(self.TempPsi,0,0)
		self.TempPsi2.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
		newControl = self.TempPsi2.InnerProduct(self.TempPsi)
		fieldValue = -imag(newControl) / (2.0 * self.EnergyPenalty)
		
		# Warn if dt**2 * fieldValue is to great
		errorOrder = abs(self.TimeStep**2 * fieldValue)
		if(errorOrder > 1e-6 and self.Debug == True):
			print "Warning: Omitted dexp terms not small (", errorOrder, ")!"
	
		# Check if control amplitude exceeds control cutoff, and store
		#fieldValue = min(abs(self.ControlCutoff), abs(fieldValue)) * sign(fieldValue)

		#Update control function. We apply a window function to
		#force control to be zero at beginning and end. Since we
		#are potentially making a _lot_ of integration steps in
		#some cases, the current time as reported by the propagator
		#may drift off the mark by several orders of magnitude. To
		#avoid the sin()-window going negative, we take the max()
		#of it with 0.
		#T = self.BaseProblem.Duration
		#mask = max(sin(pi * min(t,T) / (T - self.TimeStep)), 0.0)
		mask = 1.0
		self.ControlVector[timeGridIndex] = sqrt(mask) * fieldValue


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
		
		#Initial step
		if self.currIter > 0:
			self.ComputeNewControlFunction(0,0.0)
		self.ControlFunction.ConfigSection.e0 = self.ControlVector[0]

		for idx, t in enumerate(self.BaseProblem.Advance(self.TimeGridSize)):
			if idx < self.TimeGridSize - 1:
				if self.currIter > 0:
					self.ComputeNewControlFunction(idx + 1, t)
				self.ControlFunction.ConfigSection.e0 = self.ControlVector[idx + 1]

	
	def BackwardPropagation(self, targetProjection):
		"""
		Propagate backwards from T to 0. Terminal condition is Psi(T) = <F|Psi_forward(T)> |F>,
		where |F> is the desired final state.
		"""

		print "Backward propagation, pass ", self.currIter

		self.SetupStep(Direction.Backward)

		# Set the initial state
		self.BaseProblem.psi.GetData()[:] = targetProjection * self.TargetState.GetData()[:]
		
		#Set control and store backward solution for t = T
		self.ControlFunction.ConfigSection.e0 = self.ControlVector[self.TimeGridSize - 1]
		#self.BackwardSolution[:, self.TimeGridSize - 1] = self.BaseProblem.psi.GetData()[:]

		T = self.PropagationTime
		for idx, t in enumerate(self.BaseProblem.Advance(self.TimeGridSize, duration = T)):
			self.BackwardSolution[:, self.TimeGridSize - idx - 1] = self.BaseProblem.psi.GetData()[:]
			if idx < self.TimeGridSize - 1:
				self.ControlFunction.ConfigSection.e0 = self.ControlVector[self.TimeGridSize - idx - 2]


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



class Direction:
	Forward = 0
	Backward = 1
