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

		#Get the list of control function potentials
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

		self.ControlsAreCommuting = True

	def Setup(self):
		self.M = numpy.zeros((self.NumberOfControls, self.NumberOfControls), dtype=double)
		self.b = numpy.zeros(self.NumberOfControls, dtype=double)

		self.TempPsi3 = self.BaseProblem.psi.CopyDeep()
		self.TempPsi4 = self.BaseProblem.psi.CopyDeep()

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
		
		#Should we skip update on backward propagation?
		if direction == Direction.Backward and not self.UpdateBackward:
			self.UpdateControls(self.ControlVectors[:,timeGridIndex])
		
		#Set up linear system of equations for the new controls
		self.SetupVectorB(timeGridIndex, direction)
		self.SetupMatrixM(timeGridIndex, direction)

		#Solve lin. system to find new controls
		newControls = linalg.solve(self.M, self.b)

		#Update controls
		self.UpdateControls(newControls)

	
	def UpdateControls(self, newControls):
		"""
		Update controls from list newControls
		"""
		for a in range(self.NumberOfControls):
			self.ControlVectors[a, timeGridIndex] = newControls[a]
			self.ControlFunctionList[a].ConfigSection.strength = newControls[a]


	def SetupVectorB(self, timeGridIndex, direction):
		"""
		Case 1: Controls are commuting
		Case 2: Controls are not commuting
		"""

		#Calculate the constants to be multiplied on the two
		#commutators that appear in the expression for b
		if direction == Direction.Forward:
			commutatorScaling1 = 1j * self.TimeStep / 2.0
			commutatorScaling2 = -self.TimeStep**2 / 6.0
		else:
			commutatorScaling1 = -1j * self.TimeStep / 2.0
			commutatorScaling2 = -self.TimeStep**2 / 6.0
		
		#Loop over all controls
		for a in range(self.NumberOfControls):
			self.ControlFunctionList[a].ConfigSection.strength = 1.0

			#CASE 1: Controls are commuting
			if self.ControlsAreCommuting:
				#X * |psi>
				self.TempPsi.Clear()
				self.ControlFunctionList[a].MultiplyPotential(self.Psi, self.TempPsi,0,0)
				self.TempPsi4.GetData()[:] = self.TempPsi.GetData()[:]

				#ih/2 * [H_0, X] * |psi>
				self.MultiplyCommutatorAB(self.H_0, self.ControlFunctionList[a], self.Psi, self.TempPsi, self.TempPsi2)
				self.TempPsi4.GetData()[:] +=  commutatorScaling1 * self.TempPsi2.GetData()[:]

				# h/6 * [H_0, [H_0, X]] * |psi>
				self.MultiplyCommutatorAAB(self.H_0, self.ControlFunctionList[a], self.Psi, self.TempPsi, self.TempPsi2, self.TempPsi3)
				self.TempPsi4.GetData()[:] += commutatorScaling2 * self.TempPsi3.GetData()[:]

			#CASE 2: Controls are not commuting
			else:
				pass #for now

			#Now perform innerproduct
			if direction == Direction.Forward:
				self.TempPsi.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
				self.b[a] = -numpy.imag(self.TempPsi.InnerProduct(self.TempPsi4))
			else:
				self.TempPsi.GetData()[:] = self.ForwardSolution[:, timeGridIndex]
				self.b[a] = -numpy.imag(self.TempPsi4.InnerProduct(self.TempPsi))

			#Last, the penalty matrix
			if not self.PenaltyMatrixIsDiagonal:
				pass #do something here


	def SetupMatrixM(self, timeGridIndex, direction):
		"""
		Case 1: Only one control
		Case 2: Controls are commuting
		Case 3: Controls are not commuting
		"""

		self.M[:] = 0

		#Penalty matrix factor
		energyPenalty = self.PenaltyMatrix[timeGridIndex] / self.TimeStep

		commutatorScaling = -self.TimeStep**2 / 6.0
		
		#CASE 1: One control
		if self.NumberOfControls == 1:
			#Commutator part h**2/6 [A,[B,A]]
			self.MultiplyCommutatorABA(self.ControlFunctionList[0], self.H_0, self.Psi, self.TempPsi, self.TempPsi2, self.TempPsi3)

			matrixElementM = 0
			if direction == Direction.Forward:
				self.TempPsi.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
				matrixElementM -= numpy.imag(commutatorScaling * self.TempPsi.InnerProduct(self.TempPsi3))
			else:
				self.TempPsi.GetData()[:] = self.ForwardSolution[:, timeGridIndex]
				matrixElementM -= numpy.imag(commutatorScaling * self.TempPsi3.InnerProduct(self.TempPsi))

			self.M[0,0] = energyPenalty + matrixElementM

		#CASE 2: Controls are commuting
		elif self.ControlsAreCommuting:
			for a in range(self.NumberOfControls):
				for b in range(self.NumberOfControls):

					#Diagonal commutator part h**2/6 [X,[H_0,Y]]
					if a == b:
						self.MultiplyCommutatorABA(self.ControlFunctionList[0], self.H_0, self.Psi, self.TempPsi, self.TempPsi2, self.TempPsi3)
						
						#Energy penalty
						self.M[a,b] += energyPenalty
					
					#Off-diagonal commutator part h**2/6 [X,[H_0,X]]
					else:
						self.MultiplyCommutatorABC(self.ControlFunctionList[a], self.H_0, self.ControlFunctionList[b], self.Psi, \
							self.TempPsi, self.TempPsi2, self.TempPsi3)
				
					#Project on forward or backward solution
					if direction == Direction.Forward:
						self.TempPsi.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
						self.M[a,b] -= numpy.imag(commutatorScaling * self.TempPsi.InnerProduct(self.TempPsi3))
					else:
						self.TempPsi.GetData()[:] = self.ForwardSolution[:, timeGridIndex]
						self.M[a,b] -= numpy.imag(commutatorScaling * self.TempPsi3.InnerProduct(self.TempPsi))

		#CASE 3: Controls are not commuting
		else:
			pass #for now



