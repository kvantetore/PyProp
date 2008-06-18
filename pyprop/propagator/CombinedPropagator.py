
class CombinedPropagator(PropagatorBase):
	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()

	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)
	
		#Create all sub propagators
		self.SubPropagators = []
		for i in range(self.Rank):
			sectionName = config.Propagation.Get("propagator" + str(i))
			if sectionName != None:
				section = config.GetSection(sectionName)

				#Create sub-propagator
				prop = section.propagator(self.psi, i)
				print "Propagator for rank %i is %s" % (i, prop)
			
				#Apply config to sub propagator
				config.Apply(prop)
				section.Apply(prop)
				
				self.SubPropagators.append(prop)
			
	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):

		for prop in self.SubPropagators:
			#parallelization
			if prop.TransformRank == self.Rank-1 and not IsSingleProc():
				self.SetupTranspose()
				self.Transpose(1, self.psi)
			prop.SetupStep(dt/2.)

		self.SetupPotential(dt)

		for prop in reversed(self.SubPropagators):
			prop.SetupStepConjugate(dt/2.)
			#parallelization
			if prop.TransformRank == self.Rank-1 and not IsSingleProc():
				self.Transpose(2, self.psi)

		if not IsSingleProc():
			distr = self.psi.GetRepresentation().GetDistributedModel().GetDistribution()
			if len(distr) != 1:
				raise "Invalid distribution length %i. Distribution (%s) is invalid" % (len(distr), distr)
			if distr[0] != self.Rank-1:
				raise "Distribution (%s) is invalid" % (distr)


	def MultiplyHamiltonian(self, destPsi, t, dt):
		#first halfstep
		for prop in self.SubPropagators:
			#parallelization
			if prop.TransformRank == self.Rank-1:
				self.Transpose(1, self.psi)
				self.Transpose(1, destPsi)

			#advance step
			prop.MultiplyHamiltonian(destPsi, t, dt/2.)

		#the potential is in the middle, so it does one full step
		self.MultiplyPotential(self.psi, destPsi, t, dt)

		#second halfstep
		for prop in reversed(self.SubPropagators):
			#advance step
			prop.MultiplyHamiltonianConjugate(destPsi, t, dt/2.)
			#parallelization
			if prop.TransformRank == self.Rank-1:
				self.Transpose(2, self.psi)
				self.Transpose(2, destPsi)
	
	def CalculatePotentialExpectationValue(self, tmpPsi, potential, t, dt):
		#first halfstep
		for prop in self.SubPropagators:
			#parallelization
			if prop.TransformRank == self.Rank-1:
				self.Transpose(1, self.psi)
				self.Transpose(1, tmpPsi)

			#transform step
			prop.InverseTransform(self.psi)
			prop.InverseTransform(tmpPsi)

		#multiply potential
		tmpPsi.Clear()
		potential.MultiplyPotential(self.psi, tmpPsi, t, dt)

		#second halfstep
		for prop in reversed(self.SubPropagators):
			#transform step
			prop.ForwardTransform(self.psi)
			prop.ForwardTransform(tmpPsi)

			#parallelization
			if prop.TransformRank == self.Rank-1:
				self.Transpose(2, self.psi)
				self.Transpose(2, tmpPsi)

		#inner product yields 
		return prop.psi.InnerProduct(tmpPsi)

	def AdvanceStep(self, t, dt):
		#Advance one step using strang splitting

		#first halfstep
		for prop in self.SubPropagators:
			#parallelization
			if prop.TransformRank == self.Rank-1:
				self.Transpose(1, self.psi)
			#advance step
			prop.AdvanceStep(t, dt/2.)

		#the potential is in the middle, so it does one full step
		self.ApplyPotential(t, dt)

		#second halfstep
		for prop in reversed(self.SubPropagators):
			#advance step
			prop.AdvanceStepConjugate(t, dt/2.)
			#parallelization
			if prop.TransformRank == self.Rank-1:
				self.Transpose(2, self.psi)

	#Transpose
	def SetupTranspose(self):
		#get transpose
		distrModel = self.psi.GetRepresentation().GetDistributedModel()
		if not distrModel.IsSingleProc():
			self.Distribution1 = distrModel.GetDistribution().copy()
			self.Distribution2 = GetAnotherDistribution(self.Distribution1, self.Rank)
			if len(self.Distribution1) > 1: 
				raise "Does not support more than 1D proc grid"

			transpose = distrModel.GetTranspose()
			#Setup shape	
			fullShape = self.psi.GetRepresentation().GetFullShape()
			distribShape = transpose.CreateDistributedShape(fullShape, self.Distribution2)
		
			self.DistributedShape1 = array(self.psi.GetData().shape)
			self.DistributedShape2 = distribShape

	def Transpose(self, stage, psi):
		distrModel = psi.GetRepresentation().GetDistributedModel()
		if not distrModel.IsSingleProc():
			if stage == 1:
				newDistrib = self.Distribution2
				newShape = self.DistributedShape2
				#Get transpose buffer
				transposeBufferName = psi.GetAvailableDataBufferName(newShape)
				if transposeBufferName == -1:
					transposeBufferName = psi.AllocateData(newShape)
				#Change Distribution
				distrModel.ChangeDistribution(psi, newDistrib, transposeBufferName)
				for i in range(self.Rank):
					psi.GetRepresentation().GetRepresentation(i).GetDistributedModel().SetDistribution(newDistrib)

			elif stage == 2:
				newDistrib = self.Distribution1
				newShape = self.DistributedShape1
				#Get transpose buffer
				transposeBufferName = psi.GetAvailableDataBufferName(newShape)
				if transposeBufferName == -1:
					transposeBufferName = psi.AllocateData(newShape)
				#Change Distribution
				distrModel.ChangeDistribution(psi, newDistrib, transposeBufferName)
				for i in range(self.Rank):
					psi.GetRepresentation().GetRepresentation(i).GetDistributedModel().SetDistribution(newDistrib)

			else:
				raise "Invalid stage %i" % stage
		
			for i in range(self.Rank):
				subDistribModel = self.psi.GetRepresentation().GetRepresentation(i).GetDistributedModel()
				subDistribModel.SetDistribution(newDistrib)


	def GetBasisFunction(self, rank, basisIndex):
		prop = self.SubPropagators[rank]
		return prop.GetBasisFunction(rank, basisIndex)


	
