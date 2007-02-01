
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

		self.KineticPotential = CreateInstanceRank("core.StaticPotential", rank)
		self.SplittingOrder = 2
		self.RepresentationMapping = dict()

	def ApplyConfig(self, config): 
		PropagatorBase.ApplyConfig(self, config)
		
	def ApplyConfigSection(self, configSection): 
		PropagatorBase.ApplyConfigSection(self, configSection)

	def SetupStep(self, dt):
		#representation mappings		
		startRepr = self.psi.GetRepresentation()
		fftRepr = self.FFTTransform.CreateFourierRepresentation(startRepr)
		self.RepresentationMapping[startRepr] = fftRepr
		self.RepresentationMapping[fftRepr] = startRepr

		#Set up potentials
		if self.SplittingOrder == 1:
			self.SetupPotential(dt/2.)
		
		elif self.SplittingOrder == 2:
			self.SetupPotential(dt/2.)
	
		else:
			raise "Invalid splitting order " + str(self.SplittingOrder)
		
		
		# transform into fourier space
		print "        Setting up forward transform"
		if IsSingleProc():
			self.FFTTransform.ForwardTransform(self.psi)
		else:
			rank = self.TransformA(self.FFT_FORWARD)
			self.TransformDistribution()
			self.TransformB(rank, self.FFT_FORWARD)
		self.ChangeRepresentation()
				
		# set up potential
		self.SetupKineticPotential(dt)
		
		# transform back into real space
		print "        Setting up inverse transform"
		if IsSingleProc():
			self.FFTTransform.InverseTransform(self.psi)
		else:
			self.TransformB(rank, self.FFT_BACKWARD)
			self.TransformDistribution()
			self.TransformA(self.FFT_BACKWARD)
			self.TransformNormalize()
		self.ChangeRepresentation()
		
	def AdvanceStep(self, t, dt):
		#Apply potential
		if self.SplittingOrder == 1:
			#first order splitting
			self.ApplyPotential(t, dt)
			self.AdvanceKineticEnergy()
		
		elif self.SplittingOrder == 2:
			#strang splitting
			self.ApplyPotential(t, dt/2.)
			self.AdvanceKineticEnergy()
			self.ApplyPotential(t, dt/2.)

		else:
			raise "Invalid splitting order " + str(self.SplittingOrder)
	
		

	def AdvanceKineticEnergy(self):
		# transform into fourier space
		if IsSingleProc():
			self.FFTTransform.ForwardTransform(self.psi)
		else:
			rank = self.TransformA(self.FFT_FORWARD)
			self.TransformDistribution()
			self.TransformB(rank, self.FFT_FORWARD)
		#change representation
		self.ChangeRepresentation()
		
		# apply kinetic energy potential
		self.KineticPotential.ApplyPotential(self.psi)
		
		# transform back into real space
		if IsSingleProc():
			self.FFTTransform.InverseTransform(self.psi)
		else:
			self.TransformB(rank, self.FFT_BACKWARD)
			self.TransformDistribution()
			self.TransformA(self.FFT_BACKWARD)
			self.TransformNormalize()
		self.ChangeRepresentation()

	def SetupKineticPotential(self, dt):
		#create config
		class staticEnergyConf(Section):
			def __init__(self, type, classname):
				self.type = type
				self.classname = classname
		conf = staticEnergyConf(PotentialType.Static, "CartesianKineticEnergyPotential")

		#create potential 
		pot = CreatePotentialFromSection(conf, "KineticEnergy", self.psi)
		pot.SetupStep(dt)
		self.KineticPotential = pot.Potential
		
	def IsDistributedRank(self, rank):
		return self.psi.GetRepresentation().GetDistributedModel().IsDistributedRank(self.psi, rank)
		
	def TransformA(self, direction):
		rank = self.psi.GetRepresentation().GetDistributedModel().GetDistributedRank(self.psi)
		self.FFTTransform.TransformExceptDistributedRank(self.psi, direction)
		return rank
			
	def TransformRepresentation(self):
		self.psi.GetRepresentation().GetDistributedModel().ChangeRepresentation(self.psi)
	
	def TransformB(self, rank, direction):
		self.FFTTransform.TransformRank(self.psi, rank, direction)

	#TransformB is not a very nice name, use this instead...
	def TransformRank(self, rank, direction):
		self.TransformB(rank, direction)
		
	def TransformNormalize(self):
		self.FFTTransform.Renormalize(self.psi)
		
	def ChangeRepresentation(self, direction=1):
		oldRepr = self.psi.GetRepresentation()

		#forward change of representation
		if direction == 1:
			if not self.RepresentationMapping.has_key(oldRepr):
				raise "Unknown representation ", oldRepr
				
			newRepr = self.RepresentationMapping[oldRepr]
			self.psi.SetRepresentation(newRepr)
		
		#Reverse change of representation
		if direction == -1:
			for key, value in self.RepresentationMapping.items():
				if value == oldRepr:
					self.psi.SetRepresentation(key)
					return
		
			raise "Unknown representation", oldRepr
