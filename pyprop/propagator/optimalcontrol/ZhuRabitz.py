from numpy import linalg

class ZhuRabitz(OptimalControl):
	"""
	This class implements Zhu and Rabitz's modified Krotov method, which is an optimization
	algorithm for quantum control. The user may specify several control functions which are 
	used to guide the system from a specified initial to final state. 

	More precisely, given the TDSE:
	
	     d            /                        \ 
	    ---- Psi = -i | H_0 + sum( u_i(t)X_i ) | Psi
	     dt           \                        / 
	
	for a given H_0 (time independent) and controls (scalars) u_i, Krotov solves

	    grad(J) = 0

	where J =  <Psi|P*Psi> - mu \sum( \int_0^T |u_i(t)|**2 dt ), and P is a projection operator
	"""

	BASE = OptimalControl
	
	def Setup(self):
		self.BASE.Setup(self)

		self.M = numpy.zeros((self.NumberOfControls, self.NumberOfControls), dtype=double)
		self.b = numpy.zeros(self.NumberOfControls, dtype=double)

		self.TempPsi3 = self.BaseProblem.psi.CopyDeep()

		#Find h0 potential
		for potential in self.PotentialList:
			if potential.Name == self.Config.ZhuRabitz.h0[0]:
				self.H_0 = potential


	def SetupVectorB(self, timeGridIndex, direction):
		"""
		Compute b vector:
		
		     b =  -Im(<psi| (X_a +/- ih/2 * [H_0, X_a]) |etta>)
		"""

		if direction == Direction.Forward:
			commutatorScaling = 1j * self.TimeGridResolution / 2.0
		else:
			commutatorScaling = -1j * self.TimeGridResolution / 2.0

		for a in range(self.NumberOfControls):

			#Set potential strength to 1 (operator multiplication only)
			self.ControlFunctionList[a].ConfigSection.strength = 1.0

			#First we compute X_a * |psi>.
			self.TempPsi.Clear()
			self.TempPsi2.Clear()
			self.ControlFunctionList[a].MultiplyPotential(self.Psi, self.TempPsi, 0, 0)
			self.TempPsi2.GetData()[:] = self.TempPsi.GetData()[:]

			#Then we do  H_0 * X_a * |psi>.
			self.MultiplyCommutatorAB(self.H_0, self.ControlFunctionList[a], self.Psi, self.TempPsi, self.TempPsi3)
			self.TempPsi2.GetData()[:] +=  commutatorScaling * self.TempPsi3.GetData()[:]

			#Now inner product
			if direction == Direction.Forward:
				self.TempPsi.GetData()[:] = self.BackwardSolution[:, timeGridIndex]
				self.b[a] = -numpy.imag(self.TempPsi.InnerProduct(self.TempPsi2))
			else:
				self.TempPsi.GetData()[:] = self.ForwardSolution[:, timeGridIndex]
				self.b[a] = -numpy.imag(self.TempPsi2.InnerProduct(self.TempPsi))

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



