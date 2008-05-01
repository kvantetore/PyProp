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

	BASE = OptimalControl

	def ApplyConfigSection(self, config):
		self.BASE.ApplyConfigSection(self, config)

		#Are controls commuting?
		self.ControlsAreCommuting = True
		if hasattr(self.Config, "controls_are_commuting"):
			self.ControlsAreCommuting = self.Control.controls_are_commuting

		#Find h0 potential
		for potential in self.PotentialList:
			if potential.Name == config.h0[0]:
				self.H_0 = potential

	def Setup(self):
		self.BASE.Setup(self)
		
		self.M = numpy.zeros((self.NumberOfControls, self.NumberOfControls), dtype=double)
		self.b = numpy.zeros(self.NumberOfControls, dtype=double)

		self.TempPsi3 = self.BaseProblem.psi.CopyDeep()
		self.TempPsi4 = self.BaseProblem.psi.CopyDeep()


	def SetupVectorB(self, timeGridIndex, direction):
		"""
		Case 1: Controls are commuting
		Case 2: Controls are not commuting
		"""

		#Calculate the constants to be multiplied on the two
		#commutators that appear in the expression for b
		commutatorScaling2 = -self.TimeStep**2 / 6.0
		if direction == Direction.Forward:
			commutatorScaling1 = 1j * self.TimeStep / 2.0
		else:
			commutatorScaling1 = -1j * self.TimeStep / 2.0
			commutatorScaling2 = -self.TimeStep**2 / 6.0
		
		#Loop over all controls
		for a in range(self.NumberOfControls):
			self.ControlFunctionList[a].ConfigSection.strength = 1.0
			
			#
			#CASE 1: Controls are commuting
			#
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

			#
			#CASE 2: Controls are not commuting
			#
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
		
		#Clear matrix M
		self.M[:] = 0

		#Penalty matrix factor
		energyPenalty = self.PenaltyMatrix[timeGridIndex] / self.TimeGridResolution

		commutatorScaling = -self.TimeStep**2 / 6.0
	
		#
		#CASE 1: One control
		#
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

		#
		#CASE 2: Controls are commuting
		#
		elif self.ControlsAreCommuting:
			for a in range(self.NumberOfControls):
				for b in range(self.NumberOfControls):

					#Diagonal commutator part h**2/6 [X,[H_0,Y]]
					if a == b:
						self.MultiplyCommutatorABA(self.ControlFunctionList[0], self.H_0, \
							self.Psi, self.TempPsi, self.TempPsi2, self.TempPsi3)
						
						#Energy penalty
						self.M[a,b] += energyPenalty
					
					#Off-diagonal commutator part h**2/6 [X,[H_0,X]]
					else:
						self.MultiplyCommutatorABC(self.ControlFunctionList[a], self.H_0, self.ControlFunctionList[b], \
							self.Psi, self.TempPsi, self.TempPsi2)
				
					#Project on forward or backward solution
					if direction == Direction.Forward:
						self.TempPsi.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
						self.M[a,b] -= numpy.imag(commutatorScaling * self.TempPsi.InnerProduct(self.TempPsi3))
					else:
						self.TempPsi.GetData()[:] = self.ForwardSolution[:, timeGridIndex]
						self.M[a,b] -= numpy.imag(commutatorScaling * self.TempPsi3.InnerProduct(self.TempPsi))

		#
		#CASE 3: Controls are not commuting
		#
		else:
			for a in range(self.NumberOfControls):
				for b in range(self.NumberOfControls):

					#Diagonal commutator part h**2/6 [X,[H_0,Y]]
					if a == b:
						#Energy penalty
						self.M[a,b] += energyPenalty
					
					#Off-diagonal commutator part h/2 [X_b,X_a]
					else:
						self.MultiplyCommutatorAB(self.ControlFunctionList[b], self.ControlFunctionList[a], \
							self.Psi, self.TempPsi, self.TempPsi2, self.TempPsi3)
				
					#Project on forward or backward solution
					if direction == Direction.Forward:
						self.TempPsi.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
						self.M[a,b] -= numpy.imag(commutatorScaling * self.TempPsi.InnerProduct(self.TempPsi3))
					else:
						self.TempPsi.GetData()[:] = self.ForwardSolution[:, timeGridIndex]
						self.M[a,b] -= numpy.imag(commutatorScaling * self.TempPsi3.InnerProduct(self.TempPsi))





