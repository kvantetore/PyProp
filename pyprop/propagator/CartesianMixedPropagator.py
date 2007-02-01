
class CartesianMixedPropagator(CartesianPropagator):
	__Base = CartesianPropagator

	def __init__(self, psi):
		self.__Base.__init__(self, psi)
		self.MixedPotentials = dict()
		del self.KineticPotential
		
	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)
		self.Config = config

	def ApplyConfigSection(self, configSection):
		self.__Base.ApplyConfigSection(self, configSection)


	def SetupStep(self, dt):
		#Setup representations
		print "    Setting up representations"
		gridRepr = self.psi.GetRepresentation()
		reprList = list()
		reprList.append(gridRepr)
		for rank in r_[0:self.psi.GetRank()]:
			print "        Setting up fourier representation ", rank
			fftReprRank = self.FFTTransform.CreateFourierRepresentation(gridRepr, rank)
			self.RepresentationMapping[reprList[rank]] = fftReprRank
			reprList.append(fftReprRank)
		self.RepresentationMapping[reprList[len(reprList)-1]] = gridRepr
		del reprList

		print "    Setting up potential"
		self.SetupPotential(dt)

		print "    Setting up mixed potentials"
		for rank in r_[0:self.psi.GetRank()]:
			self.ChangeTransformedRank(rank)

			#Set up the mixed potential dependant on p_{rank} and 
			#optionally, all grid coordinates except from x_{rank}
			repr = self.psi.GetRepresentation()
			potential = CreatePotential(self.Config, "MixedPotential_" + str(rank), self.psi)
			self.MixedPotentials[repr] = potential
			potential.SetupStep(dt)
	
		#Transform the last rank back to grid coordinated
		self.ChangeRepresentation()
		self.TransformRank(self.psi.GetRank()-1, self.FFT_BACKWARD)


	def AdvanceStep(self, t, dt):
		#Advancing potential
		self.ApplyPotential(t, dt)
		
		#Advancing mixed potentials
		for rank in r_[0:self.psi.GetRank()]:
			self.ChangeTransformedRank(rank)

			#Apply potential
			repr = self.psi.GetRepresentation()
			if not self.MixedPotentials.has_key(repr):
				raise "Missing mixed potential for rank " + str(rank)
			potential = self.MixedPotentials[repr]
			#print "Applying potential ", potential.Name
			potential.AdvanceStep(t, dt)

		#Transform the last rank back to grid coordinated
		self.ChangeRepresentation(1)
		self.TransformRank(self.psi.GetRank()-1, self.FFT_BACKWARD)	
		
		self.TransformNormalize()


	def ChangeTransformedRank(self, rank, direction=1):
		#Make sure that this rank is not currently distributed
		if self.IsDistributedRank(rank):
			self.TransformDistribution()

		#update the wavefunction with this representation
		self.ChangeRepresentation(direction)

		#Transform back the previous fourier transformed rank
		if rank - direction*1 >= 0:
			self.TransformRank(rank-direction*1, self.FFT_BACKWARD)
		#Transform this rank into fourier space
		self.TransformRank(rank, self.FFT_FORWARD)
