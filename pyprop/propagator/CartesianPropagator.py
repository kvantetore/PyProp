
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
		self.__WavefunctionData = {}

	def ApplyConfig(self, config): 
		PropagatorBase.ApplyConfig(self, config)
		
	def ApplyConfigSection(self, configSection): 
		PropagatorBase.ApplyConfigSection(self, configSection)
		self.Mass = 1.0
		if hasattr(configSection, 'mass'):
			self.Mass = configSection.mass

	def SetupStep(self, dt):
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
			distribRanks = psi.GetRepresentation().GetDistributedModel().GetDistribution().copy()
			for rank in range(0, psi.GetRank()):
				if rank not in distribRanks:
					self.TransformRank(rank, self.FFT_FORWARD, psi)
			self.Transpose(1, psi)
			for rank in distribRanks:
				self.TransformRank(rank, self.FFT_FORWARD, psi)

		self.SetFourierRepresentation(psi)

	def TransformInverse(self, psi):
		# transform back into real space
		if IsSingleProc():
			self.FFTTransform.InverseTransform(psi)
		else:
			distribRanks = psi.GetRepresentation().GetDistributedModel().GetDistribution().copy()
			for rank in range(0, psi.GetRank()):
				if rank not in distribRanks:
					self.TransformRank(rank, self.FFT_BACKWARD, psi)
			self.Transpose(2, psi)
			for rank in distribRanks:
				self.TransformRank(rank, self.FFT_BACKWARD, psi)
			self.TransformNormalize(psi)

		self.SetGridRepresentation(psi)


	def AdvanceKineticEnergy(self, t, dt):
		# transform into fourier space
		self.TransformForward(self.psi)
		
		# apply kinetic energy potential
		self.KineticPotential.AdvanceStep(t, dt) 
		
		# transform back into real space
		self.TransformInverse(self.psi)

	def MultiplyKineticEnergy(self, destPsi, t, dt):
		# transform into fourier space
		if ProcId == 0:
			pass
#			print "1 distribtuion (psi) = ", self.psi.GetRepresentation().GetDistributedModel().GetDistribution()
#			print "1 distribtuion (tempPsi) = ", self.psi.GetRepresentation().GetDistributedModel().GetDistribution()
		self.TransformForward(self.psi)
		if ProcId == 0:
			pass
#			print "2 distribtuion (psi) = ", self.psi.GetRepresentation().GetDistributedModel().GetDistribution()
#			print "2 distribtuion (tempPsi) = ", self.psi.GetRepresentation().GetDistributedModel().GetDistribution()
		self.TransformForward(destPsi)
		if ProcId == 0:
			pass
#			print "3 distribtuion (psi) = ", self.psi.GetRepresentation().GetDistributedModel().GetDistribution()
#			print "3 distribtuion (tempPsi) = ", self.psi.GetRepresentation().GetDistributedModel().GetDistribution()
	
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
			self.Distribution2 = array([0])
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
		self.FFTTransform.TransformRank(psi, int(rank), direction)
		
	def TransformNormalize(self, psi):
		self.FFTTransform.Renormalize(psi)

	def __GetPropagatorData(self, psi):
		class PropagatorData(object):
			def __init__(self):
				self.GridRepresentation = None
				self.FourierRepresentation = None

		if not psi in self.__WavefunctionData:
			propagatorData = PropagatorData()
			self.__WavefunctionData[psi] = propagatorData

		return self.__WavefunctionData[psi]

	def SetGridRepresentation(self, psi):
		data = self.__GetPropagatorData(psi)
		if data.GridRepresentation == None:
			raise RuntimeError("Grid Representation is not set. Can not continue")
	
		repr = data.GridRepresentation
		#keep the current distributed model
		distrib = psi.GetRepresentation().GetDistributedModel()
		repr.SetDistributedModel(distrib)
		psi.SetRepresentation(repr)

	def SetFourierRepresentation(self, psi):
		data = self.__GetPropagatorData(psi)
		if data.FourierRepresentation == None:
			gridRepr = psi.GetRepresentation()
			fourierRepr = self.FFTTransform.CreateFourierRepresentation(gridRepr)
			data.GridRepresentation = gridRepr
			data.FourierRepresentation = fourierRepr
	
		repr = data.FourierRepresentation
		#keep the current distributed model
		distrib = psi.GetRepresentation().GetDistributedModel()
		repr.SetDistributedModel(distrib)
		psi.SetRepresentation(repr)


