def GetAnotherDistribution(distrib, rank):
	if rank <= len(distrib):
		raise Exception("Can not have distribution with length (%i) >= rank (%i)" % (rank, len(distrib)))

	ranks = [j for j in r_[0:rank] if j not in distrib]
	return array([min(ranks)], dtype=int)

def GetAnotherDistribution2(distrib, curRank, rank):
	if rank <= len(distrib):
		raise Exception("Can not have distribution with length (%i) >= rank (%i)" % (rank, len(distrib)))

	#find out which procRank curRank is distributed on
	curProcRank = find(array(distrib) == curRank)[0]

	#Ranks available for distribution
	availableRanks = [j for j in r_[0:rank] if j not in distrib]

	#The new distrib is similar to the old distrib, only with the distribution 
	#curProcRank changed
	newDistrib = array(distrib)
	newDistrib[curProcRank] = min(availableRanks)

	return newDistrib


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

	def MultiplyHamiltonian(self, srcPsi, destPsi, t, dt):
		#self.MultiplyPotential(destPsi, t, dt/2.)
		self.MultiplyKineticEnergy(srcPsi, destPsi, t, dt)
		self.MultiplyPotential(srcPsi, destPsi, t, dt/2.)

	def TransformForward(self, psi):
		# transform into fourier space
		if IsSingleProc():
			self.FFTTransform.ForwardTransform(psi)
		else:
			stage = 0
			for curRank in self.TransformOrder:
				distribRanks = psi.GetRepresentation().GetDistributedModel().GetDistribution().copy()
				if curRank in distribRanks:
					stage += 1
					self.Transpose(stage, psi)
				self.TransformRank(curRank, self.FFT_FORWARD, psi)

		self.SetFourierRepresentation(psi)

	def TransformInverse(self, psi):
		# transform back into real space
		if IsSingleProc():
			self.FFTTransform.InverseTransform(psi)
		else:
			stage = len(self.DistribList)-1
			for curRank in reversed(self.TransformOrder):
				distribRanks = psi.GetRepresentation().GetDistributedModel().GetDistribution().copy()
				if curRank in distribRanks:
					stage -= 1
					self.Transpose(stage, psi)
				self.TransformRank(curRank, self.FFT_BACKWARD, psi)
			self.TransformNormalize(psi)

		self.SetGridRepresentation(psi)


	def AdvanceKineticEnergy(self, t, dt):
		# transform into fourier space
		self.TransformForward(self.psi)
		
		# apply kinetic energy potential
		self.KineticPotential.AdvanceStep(t, dt) 
		
		# transform back into real space
		self.TransformInverse(self.psi)

	def MultiplyKineticEnergy(self, srcPsi, destPsi, t, dt):
		# transform into fourier space
		self.TransformForward(srcPsi)
		self.TransformForward(destPsi)
	
		# apply kinetic energy potential
		self.KineticPotential.MultiplyPotential(srcPsi, destPsi, t, dt) 
		
		# transform back into real space
		self.TransformInverse(destPsi)
		self.TransformInverse(srcPsi)

	def SetupKineticPotential(self, dt):
		#create config
		class staticEnergyConf(Section):
			def __init__(self, type, classname):
				self.type = type
				self.classname = classname
		conf = staticEnergyConf(PotentialType.Static, "CartesianKineticEnergyPotential")
		conf.mass = self.Mass
		conf.storage_model = StaticStorageModel.StorageValue

		#create potential 
		pot = CreatePotentialFromSection(conf, "KineticEnergy", self.psi)
		pot.SetupStep(dt)
		self.KineticPotential = pot
		
	#Transpose
	def SetupTranspose(self):
		#get transpose
		distrModel = self.psi.GetRepresentation().GetDistributedModel()

		if not distrModel.IsSingleProc():
			distribs = []
			transformOrder = []

			initDistrib = distrModel.GetDistribution().copy()
			fullShape = self.psi.GetRepresentation().GetFullShape()

			#the first transforms are the ones that are not distributed
			transformOrder = [r for r in r_[:self.psi.GetRank()] if not r in initDistrib]
			remainingRanks = filter(lambda r: r not in transformOrder, r_[:self.psi.GetRank()])

			#the rest are done in ascending order
			distribs.append( (initDistrib, array(self.psi.GetData().shape)) )
			for i, curRank in enumerate(remainingRanks):
				#Get the current distribution	
				startDistrib, startShape = distribs[i] 

				#Find the next distribution
				finalDistrib = GetAnotherDistribution2(startDistrib, curRank, self.psi.GetRank())
		
				#Find shape of the next distribution
				transpose = distrModel.GetTranspose()
				finalShape = transpose.CreateDistributedShape(fullShape, finalDistrib)

				#Update DistributionList
				distribs.append( (finalDistrib, finalShape) )
				transformOrder.append(curRank)

			self.DistribList = distribs	
			self.TransformOrder = transformOrder


	def Transpose(self, stage, psi):
		distrModel = psi.GetRepresentation().GetDistributedModel()
		if not distrModel.IsSingleProc():
			if 0 <= stage < len(self.DistribList):
				newDistrib, newShape = self.DistribList[stage]
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


