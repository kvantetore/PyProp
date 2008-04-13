def GetAnotherDistribution(distrib, rank):
	if rank <= len(distrib):
		raise Exception("Can not have distribution with length (%i) >= rank (%i)" % (rank, len(distrib)))

	ranks = [j for j in r_[0:rank] if j not in distrib]
	return array([min(ranks)], dtype=int)

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
			self.Distribution2 = GetAnotherDistribution(self.Distribution1, self.psi.GetRank())
			if len(self.Distribution1) > 1: 
				raise "Does not support more than 1D proc grid"

			transpose = distrModel.GetTranspose()
			#Setup shape	
			fullShape = self.psi.GetRepresentation().GetFullShape()
			print fullShape, self.Distribution2
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

			elif stage == 2:
				newDistrib = self.Distribution1
				newShape = self.DistributedShape1
				#Get transpose buffer
				transposeBufferName = psi.GetAvailableDataBufferName(newShape)
				if transposeBufferName == -1:
					transposeBufferName = psi.AllocateData(newShape)
				#Change Distribution
				distrModel.ChangeDistribution(psi, newDistrib, transposeBufferName)
			
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

	def GetFourierRepresentation(self, psi):
		data = self.__GetPropagatorData(psi)
		if data.FourierRepresentation == None:
			gridRepr = psi.GetRepresentation()
			fourierRepr = self.FFTTransform.CreateFourierRepresentation(gridRepr)
			data.GridRepresentation = gridRepr
			data.FourierRepresentation = fourierRepr
		return data.FourierRepresentation
		
	def SetFourierRepresentation(self, psi):
		repr = self.GetFourierRepresentation(psi)
		
		#keep the current distributed model
		distrib = psi.GetRepresentation().GetDistributedModel()
		repr.SetDistributedModel(distrib)

		psi.SetRepresentation(repr)


