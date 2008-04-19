from numpy import linalg

class Degani(OptimalControl):
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

	def ApplyConfigSection(self, config):
		"""
		"""
		self.Config = config.Degani

		#Get the list of control function potentialss
		self.ControlFunctionNamesList = self.Config.control_functions
		self.NumberOfControls = len(self.ControlFunctionNamesList)

		# Shortcuts to the salient control parameters
		self.ControlCutoff = self.Config.control_cutoff
		self.EnergyPenalty = self.Config.energy_penalty
		self.MaxIterations = self.Config.max_iterations
		self.YieldRequirement = self.Config.yield_requirement
		if hasattr(self.Config, "time_grid_size"):
			self.TimeGridSize = self.Config.time_grid_size
			self.TimeGridResolution = self.PropagationTime / self.TimeGridSize

		#Set control cutoff from time step
		err = 1e-6
		self.ControlCutoff = err / self.TimeStep**2

		#Check for the debug option
		self.Debug = False
		if hasattr(self.Config, "debug"):
			self.Debug = self.Config.debug

		#Do backward updates?
		self.UpdateBackward = hasattr(self.Config, "update_backwards") and self.Config.update_backwards or False

		self.PenaltyMatrixIsDiagonal = True


	def Setup(self):
		self.M = numpy.zeros((self.NumberOfControls, self.NumberOfControls), dtype=double)
		self.b = numpy.zeros(self.NumberOfControls, dtype=double)

		self.TempPsi3 = self.BaseProblem.psi.CopyDeep()

		#Find h0 potential
		for potential in self.PotentialList:
			if potential.Name == self.Config.h0[0]:
				self.H_0 = potential

	def ComputeNewControlFunctions(self, timeGridIndex, t, direction):
		"""
		Updates the control function based on backward propagation solution:

		    E(t) = - 1 / mu * imag(<etta|X1 + X2 + ...|psi>)

		where E is the new control, mu is the energy penalty, |etta> is the
		backward solution and X the control matrix.
		"""

		if direction == Direction.Backward and not self.UpdateBackward:
			for a in range(self.NumberOfControls):
				self.ControlFunctionList[a].ConfigSection.strength = self.ControlVectors[a, timeGridIndex]
			return

		self.SetupVectorB(timeGridIndex, direction)
		self.SetupMatrixM(timeGridIndex, direction)

		#newControls = linalg.solve(self.M, self.b)
		newControls = self.b[0] / self.M[0,:]

		#Update controls
		for a in range(self.NumberOfControls):
			self.ControlVectors[a, timeGridIndex] = newControls[a]
			self.ControlFunctionList[a].ConfigSection.strength = newControls[a]

	
	def ComputeCostFunctional(self, currentYield):
		penalty = 0
		for a in range(self.NumberOfControls):
			#penalty += numpy.dot(numpy.transpose(self.ControlVectors[a,:]), numpy.dot(self.PenaltyMatrix, self.ControlVectors[a,:]))
			penalty += numpy.dot(numpy.transpose(self.ControlVectors[a,:]), self.PenaltyMatrix * self.ControlVectors[a,:])
		return currentYield - penalty


	def ComputeTargetProjection(self):
		"""
		Project solution on target state
		"""

		#Note: Argument is complex conjugated
		return self.BaseProblem.psi.InnerProduct(self.TargetState)


	def SetupPenaltyMatrix(self):
		self.PenaltyMatrix = numpy.ones(self.TimeGridSize)
		self.PenaltyMatrix[:] *= self.EnergyPenalty * self.TimeGridResolution

	
	def SetupVectorB(self, timeGridIndex, direction):

		commutatorScaling = 1j * self.TimeStep / 2.0
		#self.TempPsi3.GetData()[:] = self.BaseProblem.psi.GetData()[:]
		
		for a in range(self.NumberOfControls):
			self.ControlFunctionList[a].ConfigSection.strength = 1.0

			#We must compute (X_a - ih/2 * [H_0, X_a]) * |psi>
			self.TempPsi2.Clear()
			
			#First we compute X_a * |psi>. This involves three steps:
			#
			#    1. Clear TmpPsi
			#    2. Perform X_a * |psi> and store in TmpPsi
			#    3. Store TmpPsi in TmpPsi2
			#
			self.TempPsi.Clear()
			self.ControlFunctionList[a].MultiplyPotential(self.TempPsi,0,0)
			#self.TempPsi2.GetData()[:] = self.TempPsi.GetData()[:]

			#Then we do X_a * H_0 * |psi>. The steps involved are:
			#
			#    1. Clear TmpPsi
			#    2. Perform H_0 * |psi> and store in TmpPsi
			#    3. Set |psi> = TmpPsi
			#    4. Clear TmpPsi again
			#    5. Perform X_a * |psi> and store in TmpPsi
			#    6. Store TmpPsi in TmpPsi2
			#self.TempPsi.Clear()
			#self.H_0.MultiplyPotential(self.TempPsi, 0, 0)
			#self.BaseProblem.psi.GetData()[:] = self.TempPsi.GetData()[:]
			#self.TempPsi.Clear()
			#self.ControlFunctionList[a].MultiplyPotential(self.TempPsi, 0, 0)
			#self.TempPsi2.GetData()[:] += commutatorScaling * self.TempPsi.GetData()[:]

			#H_0 * X_a * |psi>
			#self.BaseProblem.psi.GetData()[:] = self.TempPsi3.GetData()[:]
			#self.TempPsi.Clear()
			#self.ControlFunctionList[a].MultiplyPotential(self.TempPsi, 0, 0)
			#self.BaseProblem.psi.GetData()[:] = self.TempPsi.GetData()[:]
			#self.TempPsi.Clear()
			#self.H_0.MultiplyPotential(self.TempPsi, 0, 0)
			#self.TempPsi2.GetData()[:] -= commutatorScaling * self.TempPsi.GetData()[:]

			#Now inner product
			if direction == Direction.Forward:
				self.TempPsi2.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
				self.b[a] = -numpy.imag(self.TempPsi2.InnerProduct(self.TempPsi))
			else:
				self.TempPsi2.GetData()[:] = self.ForwardSolution[:, timeGridIndex]
				self.b[a] = -numpy.imag(self.TempPsi.InnerProduct(self.TempPsi2))


			#Last, the penalty matrix
			if not self.PenaltyMatrixIsDiagonal:
				pass #do something here

			#Reset |psi>
			#self.BaseProblem.psi.GetData()[:] = self.TempPsi3.GetData()[:]

	def SetupMatrixM(self, timeGridIndex, direction):

		self.M[:] = 0
		for a in range(self.NumberOfControls):
			for b in range(self.NumberOfControls):
				if a == b:
					self.M[a,b] = self.PenaltyMatrix[timeGridIndex] / self.TimeStep
	#			else:
	#				#X_a * X_b * |psi>
	#				self.TempPsi.Clear()
	#				self.H_0.MultiplyPotential(self.TempPsi, 0, 0)
	#				self.ControlFunctionList[a].MultiplyPotential(self.TempPsi,0,0)
	#				self.TempPsi2.GetData()[:] = commutatorScaling * self.TempPsi.GetData()[:]



