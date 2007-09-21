
#----------------------------------------------------------------------------------------------------
# Cartesian FFT Evaluator
#----------------------------------------------------------------------------------------------------
class CartesianPropagator(PropagatorBase):
	FFT_FORWARD = -1
	FFT_BACKWARD = 1
	
	def __init__(self, psi):
		PropagatorBase.__init__(self, psi)

		rank = psi.GetRank()
		self.FFTTransform = CreateInstanceRank("core.CartesianFourierTransform", rank)

		self.RepresentationMapping = dict()

	def ApplyConfig(self, config): 
		PropagatorBase.ApplyConfig(self, config)
		
	def ApplyConfigSection(self, configSection): 
		PropagatorBase.ApplyConfigSection(self, configSection)
		self.Mass = 1.0
		if hasattr(configSection, 'mass'):
			self.Mass = configSection.mass

	def SetupStep(self, dt):
		#representation mappings		
		startRepr = self.psi.GetRepresentation()
		fftRepr = self.FFTTransform.CreateFourierRepresentation(startRepr)
		self.RepresentationMapping[startRepr] = fftRepr
		self.RepresentationMapping[fftRepr] = startRepr

		# set up potential
		self.SetupPotential(dt/2.)

		self.SetupTranspose()
		
		# set up kinetic energy
		self.TransformForward(self.psi)
		self.SetupKineticPotential(dt)
		self.TransformInverse(self.psi)

	def AdvanceStep(self, t, dt):
		self.ApplyPotential(t, dt/2.)
		self.AdvanceKineticEnergy(t, dt)
		self.ApplyPotential(t, dt/2.)

	def MultiplyHamiltonian(self, destPsi, t, dt):
		#self.MultiplyPotential(destPsi, t, dt/2.)
		self.MultiplyKineticEnergy(destPsi, t, dt)
		self.MultiplyPotential(destPsi, t, dt/2.)

	def TransformForward(self, psi):
		# transform into fourier space
		if IsSingleProc():
			self.FFTTransform.ForwardTransform(psi)
		else:
			for rank in range(1, psi.GetRank()):
				self.TransformRank(rank, self.FFT_FORWARD, psi)
			self.Transpose(1, psi)
			self.TransformRank(0, self.FFT_FORWARD, psi)

		self.ChangeRepresentation(psi)

	def TransformInverse(self, psi):
		# transform back into real space
		if IsSingleProc():
			self.FFTTransform.InverseTransform(psi)
		else:
			self.TransformRank(0, self.FFT_BACKWARD, psi)
			self.Transpose(2, psi)
			for rank in range(1, psi.GetRank()):
				self.TransformRank(rank, self.FFT_BACKWARD, psi)
			self.TransformNormalize(psi)
		self.ChangeRepresentation(psi)


	def AdvanceKineticEnergy(self, t, dt):
		# transform into fourier space
		self.TransformForward(self.psi)
		
		# apply kinetic energy potential
		self.KineticPotential.AdvanceStep(t, dt) 
		
		# transform back into real space
		self.TransformInverse(self.psi)

	def MultiplyKineticEnergy(self, destPsi, t, dt):
		# transform into fourier space
		self.TransformForward(self.psi)
		self.TransformForward(destPsi)
		
		# apply kinetic energy potential
		self.KineticPotential.MultiplyPotential(destPsi, t, dt) 
		
		# transform back into real space
		self.TransformInverse(destPsi)
		self.TransformInverse(self.psi)

	def SetupKineticPotential(self, dt):
		#create config
		class staticEnergyConf(Section):
			def __init__(self, type, classname):
				self.type = type
				self.classname = classname
		conf = staticEnergyConf(PotentialType.Static, "CartesianKineticEnergyPotential")
		conf.mass = self.Mass

		#create potential 
		pot = CreatePotentialFromSection(conf, "KineticEnergy", self.psi)
		pot.SetupStep(dt)
		self.KineticPotential = pot
		
	#Transpose
	def SetupTranspose(self):
		#get transpose
		distrModel = self.psi.GetRepresentation().GetDistributedModel()
		if not distrModel.IsSingleProc():
			self.Distribution1 = distrModel.GetDistribution().copy()
			if len(self.Distribution1) > 1: 
				raise "Does not support more than 1D proc grid"
			transpose = distrModel.GetTranspose()
			#Setup shape	
			fullShape = self.psi.GetRepresentation().GetFullShape()
			self.Distribution2 = array([1])
			distribShape = transpose.CreateDistributedShape(fullShape, self.Distribution2)
			#allocate wavefunction
			self.TransposeBuffer1 = self.psi.GetActiveBufferName()
			self.TransposeBuffer2 = self.psi.AllocateData(distribShape)

	def Transpose(self, stage, psi):
		distrModel = psi.GetRepresentation().GetDistributedModel()
		if not distrModel.IsSingleProc():
			if stage == 1:
				distrModel.ChangeDistribution(psi, self.Distribution2, self.TransposeBuffer2)
			elif stage == 2:
				distrModel.ChangeDistribution(psi, self.Distribution1, self.TransposeBuffer1)
			else:
				raise "Invalid stage %i" % stage
		
	#Transforms one rank of psi between grid space and fourier space
	def TransformRank(self, rank, direction, psi):
		self.FFTTransform.TransformRank(psi, rank, direction)
		
	def TransformNormalize(self, psi):
		self.FFTTransform.Renormalize(psi)
		
	def ChangeRepresentation(self, psi, direction=1):
		oldRepr = psi.GetRepresentation()

		#forward change of representation
		if direction == 1:
			if not self.RepresentationMapping.has_key(oldRepr):
				raise "Unknown representation ", oldRepr
				
			newRepr = self.RepresentationMapping[oldRepr]
			psi.SetRepresentation(newRepr)
		
		#Reverse change of representation
		if direction == -1:
			for key, value in self.RepresentationMapping.items():
				if value == oldRepr:
					psi.SetRepresentation(key)
					return
		
			raise "Unknown representation", oldRepr


