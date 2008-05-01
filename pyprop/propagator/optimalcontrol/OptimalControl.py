import warnings

class OptimalControl:
	"""
	Some description here.

	PenaltyMatrix:

	    The PenaltyMatrix (A) provides a cost measure for a given control. It can be built do penalize
		a wide variety of pulse aspects, such as frequency, fluency (overall energy penalty), pulse
		shape (windoved energy penalty), etc. However, since it is a quite large matrix (if it is full),
		the number of time intervals squared, we need to keep it sparse. In practise, the types of 
		penalties can be divided into two types: diagonal in time-space and diagonal in frequency-space.
		Since the penalty is computed from u^T A u, where the controls u are time-dependent, we use
		FFT to diagonlise A for frequency penalties: u^T A u = u^T IFFT{ A' {FFT{u}} }, and we only need
		to store the diagonal matrix A'.
	"""
	
	def __init__(self, prop):
		self.BaseProblem = prop
		self.Psi = prop.psi
		self.PsiSize = prop.psi.GetData().size
		self.PotentialList = prop.Propagator.PotentialList
		if hasattr(self.BaseProblem.Propagator, "BasePropagator"):
			self.PotentialList = self.BaseProblem.Propagator.BasePropagator.PotentialList
	
		#Get time step and propagation time from baseproblem
		self.TimeStep = prop.TimeStep.real
		self.PropagationTime = prop.Duration
		self.TimeGridSize = int(round(self.PropagationTime / self.TimeStep))
		self.TimeGridResolution = self.TimeStep
		self.PerturbControl = 0


	def ApplyConfigSection(self, config):
		"""
		"""
	
		self.Config = config.Config

		#Get the list of control function potentials
		self.ControlFunctionNamesList = config.control_functions
		self.NumberOfControls = len(self.ControlFunctionNamesList)

		# Shortcuts to the salient control parameters
		self.ControlCutoff = config.control_cutoff
		self.EnergyPenalty = config.energy_penalty
		self.MaxIterations = config.max_iterations
		self.YieldRequirement = config.yield_requirement
		if hasattr(config, "time_grid_size"):
			self.TimeGridSize = config.time_grid_size
			self.TimeGridResolution = self.PropagationTime / self.TimeGridSize

		#Set control cutoff from time step
		err = 1e-6
		self.ControlCutoff = err / self.TimeStep**2

		#Check for the debug option
		self.Debug = False
		if hasattr(config, "debug"):
			self.Debug = config.debug

		#Do backward updates?
		self.UpdateBackward = hasattr(config, "update_backwards") and config.update_backwards or False

		self.PenaltyMatrixIsDiagonal = True

		self.PerturbControl = config.__dict__.get("perturb_control", 0)
		print self.PerturbControl


	def Setup(self):
		self.SetupCommon()
		self.SetupPenaltyMatrix()
		#self.SetupInitialControl()


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

		#Check that we actually got some control functions
		if len(self.ControlFunctionList) == 0:
			raise Exception("Did not find any control functions!")

		#A vector to store the optimized control functions
		self.ControlVectors = ones((self.NumberOfControls, self.TimeGridSize), dtype = double)
		for a in range(self.NumberOfControls):
			self.ControlVectors[a,:] *= self.ControlFunctionList[a].ConfigSection.strength

	
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
		Updates the control function based on backward propagation solution:

		    E(t) = - 1 / mu * imag(<etta|X1 + X2 + O(h) |psi>)

		where E is the new control, mu is the energy penalty, |etta> is the
		backward solution and X the control matrix. The number of higher order
		terms depend on the algorithm at hand (Krotov, Zhu-Rabitz or Degani)
		"""

		#Check if we should skip backward update
		if direction == Direction.Backward and not self.UpdateBackward:
			for a in range(self.NumberOfControls):
				self.ControlFunctionList[a].ConfigSection.strength = self.ControlVectors[a, timeGridIndex]

		else:

			self.SetupVectorB(timeGridIndex, direction)
			self.SetupMatrixM(timeGridIndex, direction)
			newControls = linalg.solve(self.M, self.b)

			if self.PerturbControl > 0:	
				newControls[0] += [(rand() -0.5) * self.PerturbControl]

			#Update controls
			self.UpdateControls(newControls, timeGridIndex)
	

	def ComputeCostFunctional(self, currentYield):
		penalty = 0
		for a in range(self.NumberOfControls):
			penalty += self.TimeGridResolution * numpy.dot(self.ControlVectors[a,:],  \
				self.PenaltyMatrix[:] * self.ControlVectors[a,:])
		return currentYield - penalty


	def ComputeTargetProjection(self):
		"""
		Project solution on target state
		"""
		#Note: Argument is complex conjugated
		return self.BaseProblem.psi.InnerProduct(self.TargetState)


	def UpdateControls(self, newControls, timeGridIndex):
		"""
		Update controls from list newControls
		"""
		for a in range(self.NumberOfControls):
			self.ControlVectors[a, timeGridIndex] = newControls[a]
			self.ControlFunctionList[a].ConfigSection.strength = newControls[a]


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

		for idx, t in enumerate(self.BaseProblem.Advance(self.TimeGridSize)):
			self.ForwardSolution[:, idx] = self.BaseProblem.psi.GetData()[:]
			if idx < self.TimeGridSize - 1 and self.currIter > 0:
				self.ComputeNewControlFunctions(idx+1, t, Direction.Forward)
			else:
				for a in range(self.NumberOfControls):
					self.ControlFunctionList[a].ConfigSection.strength = self.ControlVectors[a, idx]
	
	
	def BackwardPropagation(self, targetProjection):
		"""
		Propagate backwards from T to 0. Terminal condition is Psi(T) = <F|Psi_forward(T)> |F>,
		where |F> is the desired final state.
		"""

		print "Backward propagation, pass ", self.currIter

		self.SetupStep(Direction.Backward)
		#self.InitializeControls(Direction.Backward)

		# Set the initial state
		self.BaseProblem.psi.GetData()[:] = targetProjection * self.TargetState.GetData()[:]
		#self.BaseProblem.psi.GetData()[:] = self.ForwardSolution[:,-1]
		
		#Update controls
		self.ComputeNewControlFunctions(self.TimeGridSize - 1, self.PropagationTime, Direction.Backward)

		for idx, t in enumerate(self.BaseProblem.Advance(self.TimeGridSize)):
			self.BackwardSolution[:, self.TimeGridSize - 1 - idx] = self.BaseProblem.psi.GetData()[:]
			if idx < self.TimeGridSize - 1:
				self.ComputeNewControlFunctions(self.TimeGridSize - idx - 2, t, Direction.Backward)

	
	def SetupStep(self, direction):
		"""
		Set up propagator for forward or backward run
		"""

		if direction == Direction.Forward:
			self.BaseProblem.RestartPropagation(complex(self.TimeStep), 0.0, self.PropagationTime)

		elif direction == Direction.Backward:
			self.BaseProblem.RestartPropagation(-complex(self.TimeStep), self.PropagationTime, self.PropagationTime)

		#Apply initial condition
		self.BaseProblem.SetupWavefunction()

	
	def SetupPenaltyMatrix(self):
		self.PenaltyMatrix = numpy.ones(self.TimeGridSize) * self.EnergyPenalty
		#self.PenaltyMatrix[:] *= self.EnergyPenalty * self.TimeGridResolution
	

	def SetupTargetState(self):
		"""
		Set up the desired target state
		"""
		self.TargetState = self.BaseProblem.psi.CopyDeep()
		self.TargetState.Clear()

		if self.BaseProblem.Config.FinalState.type == "vector":
			
			self.TargetState.GetData()[self.BaseProblem.Config.FinalState.states] \
				= self.BaseProblem.Config.FinalState.population[:]

		elif self.BaseProblem.Config.FinalState.type == "function":
			grid = self.Psi.GetRepresentation().GetLocalGrid(0)
			func = self.BaseProblem.Config.FinalState.grid_function
			self.TargetState.GetData()[:] = func(self.BaseProblem.Config.FinalState, grid)

		else:
			raise Exception("Unknown target specification")
		
		#Normalize target state
		self.TargetState.Normalize()

	
	def SetupInitialControl(self):
		"""
		"""
		#f = None
		#if self.:
		#	f = self.
		#else:
		f = lambda omega, t, T : sin(pi * t / T)**2 * cos(omega * t)
		
		T = self.PropagationTime
		omega = 0.1
		for a in range(self.NumberOfControls):
			self.ControlVectors[a, :] = [f(omega, t, T) for t in linspace(0, T, self.TimeGridSize)]
			self.ControlVectors[a,:] *= self.ControlFunctionList[a].ConfigSection.strength

	
	def InitializeControls(self, direction):
		for a in range(self.NumberOfControls):
			self.ControlFunctionList[a].ConfigSection.strength = self.ControlVectors[a, direction]


	def MultiplyCommutatorAB(self, A, B, psi, tmpPsi, outPsi):
		"""
		[A,B] * |psi> = AB * |psi> - BA * |psi>. Result is returned
		in outPsi.
		"""

		tmpPsi.Clear()
		outPsi.Clear()

		#Calculate -BA * |psi> store in outPsi
		A.MultiplyPotential(psi, tmpPsi, 0, 0)
		B.MultiplyPotential(tmpPsi, outPsi, 0, 0)
		outPsi.GetData()[:] *= -1.0

		#Calculate AB * |psi> store in outPsi
		tmpPsi.Clear()
		B.MultiplyPotential(psi, tmpPsi, 0, 0)
		A.MultiplyPotential(tmpPsi, outPsi, 0, 0)


	def MultiplyCommutatorAAB(self, A, B, psi, tmpPsi1, tmpPsi2, outPsi):
		"""
		Multiply the commutator [A,[A,B]]] on psi. Result is returned in outPsi.
		All buffers are destroyed.

		Straightforward expansion of the commutator would require 12 matrix-vector 
		multiplications. However, we write:

		    P1 = A |psi>
		    P2 = B |psi>

		giving

		    [A,[A,B]] = AAB - 2*ABA +  BAA = AA*P2 - 2*AB*P1 + BA*P1

		where we now must do only 8 matrix-vector multiplications.
		"""
		#Clear buffers
		tmpPsi1.Clear()
		tmpPsi2.Clear()
		outPsi.Clear()

		#Calculate P1
		A.MultiplyPotential(psi, tmpPsi1, 0, 0)

		#Calculate BA*P1 and store in outPsi
		A.MultiplyPotential(tmpPsi1, outPsi, 0, 0)
		B.MultiplyPotential(outPsi, tmpPsi2, 0, 0)
		outPsi.GetData()[:] = tmpPsi2.GetData()[:]

		#Calculate -2*AB*P1 and add to outPsi
		tmpPsi2.Clear()
		B.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		A.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] -= tmpPsi1.GetData()[:]

		#Calculate P2
		tmpPsi1.Clear()
		B.MultiplyPotential(psi, tmpPsi1, 0, 0)

		#Calculate AA*P2
		tmpPsi2.Clear()
		A.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		A.MultiplyPotential(tmpPsi2, tmpPsi1,0, 0)
		outPsi.GetData()[:] += tmpPsi1.GetData()[:]


	def MultiplyCommutatorABA(self, A, B, psi, tmpPsi1, tmpPsi2, outPsi):
		"""
		Multiply the commutator [A,[B,A]]] on psi. Result is returned in outPsi.
		All buffers are destroyed. Since [A,[B,A]] = -[A,[A,B]], we just call
		on MultiplyCommutatorAAB and multiply the result by -1
		"""
		self.MultiplyCommutatorAAB(A, B, psi, tmpPsi1, tmpPsi2, outPsi)
		outPsi.GetData()[:] *= -1


	def MultiplyCommutatorABC(self, A, B, C, psi, tmpPsi1, tmpPsi2, outPsi):
		"""
		Multiply the commutator [A,[B,C]] on psi. Result is returned in outPsi.

		[A,[B,C]] = [A,BC-CB] = ABC - ACB - BCA + CBA 
		"""

		def ClearBuffers():
			tmpPsi1.Clear()
			tmpPsi2.Clear()
			outPsi.Clear()

		#Clear buffers
		ClearBuffers()

		#Calculate CBA and store in outPsi
		A.MultiplyPotential(psi, tmpPsi1, 0, 0)
		B.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		C.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] = tmpPsi1.GetData()

		#Calculate -BCA and add to outPsi
		ClearBuffers()
		A.MultiplyPotential(psi, tmpPsi1, 0, 0)
		C.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		B.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] -= tmpPsi1.GetData()[:]

		#Calculate -ACB 
		ClearBuffers()
		B.MultiplyPotential(psi, tmpPsi1, 0, 0)
		C.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		A.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] -= tmpPsi1.GetData()[:]

		#Calculate ABC
		ClearBuffers()
		C.MultiplyPotential(psi, tmpPsi1, 0, 0)
		B.MultiplyPotential(tmpPsi1, tmpPsi2, 0, 0)
		tmpPsi1.Clear()
		A.MultiplyPotential(tmpPsi2, tmpPsi1, 0, 0)
		outPsi.GetData()[:] = tmpPsi1.GetData()[:]


class Direction:
	Forward = 0
	Backward = -1
