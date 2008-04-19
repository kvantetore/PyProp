
class Krotov(OptimalControl):
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
		conf = config.Krotov
		#Get the list of control function potentialss
		self.ControlFunctionNamesList = conf.control_functions
		self.NumberOfControls = len(self.ControlFunctionNamesList)

		# Shortcuts to the salient control parameters
		self.ControlCutoff = conf.control_cutoff
		self.EnergyPenalty = conf.energy_penalty
		self.MaxIterations = conf.max_iterations
		self.YieldRequirement = conf.yield_requirement
		if hasattr(conf, "time_grid_size"):
			self.TimeGridSize = conf.time_grid_size
			self.TimeGridResolution = self.PropagationTime / self.TimeGridSize

		#Set control cutoff from time step
		err = 1e-6
		self.ControlCutoff = err / self.TimeStep**2

		#Check for the debug option
		self.Debug = False
		if hasattr(conf, "debug"):
			self.Debug = conf.debug

		#Do backward updates?
		self.UpdateBackward = hasattr(conf, "update_backwards") and conf.update_backwards or False

	def Setup(self):
		pass
	

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

		#Compute new control function. First the control matrix
		#is multiplied on the current forward solution, and then
		#the result is projected onto the backward solution.
		for a in range(self.NumberOfControls):
			self.ControlFunctionList[a].ConfigSection.strength = 1.0
			self.TempPsi.Clear()
			self.TempPsi2.Clear()
			self.ControlFunctionList[a].MultiplyPotential(self.TempPsi,0,0)
			self.TempPsi2.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
			newControl = self.TempPsi2.InnerProduct(self.TempPsi)
			fieldValue = -imag(newControl) / (2.0 * self.EnergyPenalty)
			
			# Check if control amplitude exceeds control cutoff, and store
			#fieldValue = min(abs(self.ControlCutoff), abs(fieldValue)) * sign(fieldValue)

			#Update control function. We apply a window function to
			#force control to be zero at beginning and end. Since we
			#are potentially making a _lot_ of integration steps in
			#some cases, the current time as reported by the propagator
			#may drift off the mark by several orders of magnitude. To
			#avoid the sin()-window going negative, we take the max()
			#of it with 0.
			T = self.BaseProblem.Duration
			mask = max(sin(pi * min(t,T) / (T - self.TimeStep)), 0.0)
			self.ControlVectors[a, timeGridIndex] = sqrt(mask) * fieldValue

			#Update controls
			self.ControlFunctionList[a].ConfigSection.strength = self.ControlVectors[a, timeGridIndex]

			# Warn if dt**2 * fieldValue is to great
			errorOrder = abs(self.TimeStep**2 * fieldValue)
			if(errorOrder > 1e-6 and self.Debug == True):
				print "Warning: Omitted dexp terms not small (", errorOrder, ")!"


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


