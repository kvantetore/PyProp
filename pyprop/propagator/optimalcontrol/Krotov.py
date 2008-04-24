from numpy import linalg

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
		self.Config = config.Krotov

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

	
	def SetupVectorB(self, timeGridIndex, direction):

		for a in range(self.NumberOfControls):
			self.ControlFunctionList[a].ConfigSection.strength = 1.0

			#We must compute X_a * |psi>
			self.TempPsi.Clear()
			self.TempPsi2.Clear()
			self.ControlFunctionList[a].MultiplyPotential(self.Psi, self.TempPsi, 0, 0)

			#Now inner product <X_a * psi | etta>
			if direction == Direction.Forward:
				self.TempPsi2.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
				self.b[a] = -numpy.imag(self.TempPsi2.InnerProduct(self.TempPsi))
			else:
				self.TempPsi2.GetData()[:] = self.ForwardSolution[:, timeGridIndex]
				self.b[a] = -numpy.imag(self.TempPsi.InnerProduct(self.TempPsi2))

			#Last, the penalty matrix
			if not self.PenaltyMatrixIsDiagonal:
				pass #do something here


	def SetupMatrixM(self, timeGridIndex, direction):

		self.M[:] = 0
		for a in range(self.NumberOfControls):
			for b in range(self.NumberOfControls):
				if a == b:
					self.M[a,b] = self.PenaltyMatrix[timeGridIndex] / self.TimeGridResolution
	#			else:
	#				#X_a * X_b * |psi>
	#				self.TempPsi.Clear()
	#				self.H_0.MultiplyPotential(self.TempPsi, 0, 0)
	#				self.ControlFunctionList[a].MultiplyPotential(self.TempPsi,0,0)
	#				self.TempPsi2.GetData()[:] = commutatorScaling * self.TempPsi.GetData()[:]



