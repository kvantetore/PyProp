
class CombinedPropagator(PropagatorBase):
	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.Rank = psi.GetRank()

	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)
	
		#Create all but spherical propagators
		self.SubPropagators = []
		for i in range(self.Rank):
			sectionName = config.Propagation.Get("propagator" + str(i))
			section = config.GetSection(sectionName)

			#Create sub-propagator
			prop = section.propagator(self.psi, i)
			print "Propagator for rank %i is %s" % (i, prop)
		
			#Apply config to sub propagator
			section.Apply(prop)
			
			self.SubPropagators.append(prop)
			
	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):
		for prop in self.SubPropagators:
			if prop.TransformRank == self.Rank-1:
				self.SetupTranspose()
				self.Transpose(1)
			prop.SetupStep(dt/2.)
				

		self.SetupPotential(dt)

		propagatorsReversed = list(self.SubPropagators)
		propagatorsReversed.reverse()
		for prop in propagatorsReversed:
			prop.SetupStepConjugate(dt/2.)
			if prop.TransformRank == self.Rank-1:
				self.Transpose(2)


	def AdvanceStep(self, t, dt):
		for prop in self.SubPropagators:
			if prop.TransformRank == self.Rank-1:
				self.Transpose(1)
			prop.AdvanceStep(t, dt/2.)

		self.ApplyPotential(t, dt)

		propagatorsReversed = list(self.SubPropagators)
		propagatorsReversed.reverse()
		for prop in propagatorsReversed:
			prop.AdvanceStepConjugate(t, dt/2.)
			if prop.TransformRank == self.Rank-1:
				self.Transpose(2)

	#Transpose
	def SetupTranspose(self):
		#get transpose
		distrModel = self.psi.GetRepresentation().GetDistributedModel()
		if not distrModel.IsSingleProc():
			self.Distribution1 = distrModel.GetDistribution().copy()
			if len(self.Distribution1) != 1: 
				raise "Does not support more than 1D proc grid"
			transpose = distrModel.GetTranspose()
			#Setup shape	
			fullShape = self.psi.GetRepresentation().GetFullShape()
			self.Distribution2 = array([self.Rank-2])
			distribShape = transpose.CreateDistributedShape(fullShape, self.Distribution2)
			#allocate wavefunction
			self.TransposeBuffer1 = self.psi.GetActiveBufferName()
			self.TransposeBuffer2 = self.psi.AllocateData(distribShape)

	def Transpose(self, stage):
		distrModel = self.psi.GetRepresentation().GetDistributedModel()
		if not distrModel.IsSingleProc():
			if stage == 1:
				distrModel.ChangeDistribution(self.psi, self.Distribution2, self.TransposeBuffer2)
			elif stage == 2:
				distrModel.ChangeDistribution(self.psi, self.Distribution1, self.TransposeBuffer1)
			else:
				raise "Invalid stage %i" % stage


	
