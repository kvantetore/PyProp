
class FourierPropagatorBase:
	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

		rank = psi.GetRank()	

		#instantiate transforms		
		self.FFTTransform = CreateInstanceRank("core.CartesianFourierTransform", 1)
		self.FourierTransform = CreateInstanceRank("core.RadialTransform", rank)

	def ApplyConfigSection(self, configSection):
		self.Mass = configSection.mass

		#Create all other fourier Potentials
		config = configSection.Config
		if hasattr(configSection, 'fourier_potentials'):
			potentialNames = configSection.fourier_potentials
			self.FourierPotentials = [CreatePotential(config, name, self.psi) for name in potentialNames]
		else:
			self.FourierPotentials = []

		#Create default fourier potential
		self.FourierPotentials.append(self.CreateDefaultFourierPotential())

		#Default is to not set origin to zero when imageinary time is used
		self.ForceOriginZero = False

	
	def SetupStep(self, dt):
		#setup representation
		self.GridRepresentation = self.psi.GetRepresentation().GetRepresentation(self.TransformRank)
		self.FourierRepresentation = self.CreateFourierRepr(self.GridRepresentation)

		# transform radial into fourier space
		print "	       Setup Forward Transform"
		self.TransformForward(self.psi)
		
		# set up potential
		print "	       Setup Fourier Potentials"
		for pot in self.FourierPotentials:
			pot.SetupStep(dt)

		# transform radial into grid space
		print "	       Setup Inverse Transform"
		self.TransformInverse(self.psi)

		#Get index of origin
		if self.ForceOriginZero:
			self.OriginIndex = self.GetOriginIndex()

	
	def AdvanceStep(self, t, dt):
		# transform the radial into fourier space
		self.TransformForward(self.psi)
	
		# apply radial kinetic energy potential
		for pot in self.FourierPotentials:
			pot.AdvanceStep(t, dt)
		
		# transform back into real space
		self.TransformInverse(self.psi)

		#For imaginary time propagation, set origin to 0 to keep symmetry
		if 1.0j * imag(dt) == dt and self.ForceOriginZero:
			self.SetValue(self.OriginIndex, 0)

	def MultiplyHamiltonian(self, dstPsi, t, dt):
		"""
		Calculates dstPsi += - 1/(2m) d^2/dx^2 psi
		plus any other fourier potentials which is defined.
		"""
		self.TransformForward(self.psi)
		self.TransformForward(dstPsi)
		
		for pot in self.FourierPotentials:
			pot.MultiplyPotential(self.psi, dstPsi, t, dt)

		self.TransformInverse(dstPsi)
		self.TransformInverse(self.psi)
	
	def SetupStepConjugate(self, dt):
		pass

	def AdvanceStepConjugate(self, t, dt):
		self.AdvanceStep(t, dt)

	def MultiplyHamiltonianConjugate(self, dstPsi, t, dt):
		pass
			
	def GetOriginIndex(self):
		"""
		Tries to determine what the index of the origin is
		"""
		grid = self.psi.GetRepresentation().GetLocalGrid(self.TransformRank)
		idxList = where(abs(grid) <  1e-14)[0]
		if len(idxList) != 1:
			print grid
			raise "hm!?", idxList
		else:
			idx = idxList[0]
			print "Found origin r_i = 0 at i = " + str(idx)
		return idx

	def SetValue(self, rankIndex, value):
		"""
		Sets the wavefunction to value for all points in the grid at the given
		radial index
		"""
		data = self.psi.GetData()

		#Create a slice object for all ranks except TransformRank
		index = [slice(None) for i in range(self.psi.GetRank())]
		index[self.TransformRank] = rankIndex

		data[index] = value

	def TransformForward(self, psi):
		"""
		Transforms psi from Grid space to Fourier space, and updates the representation
		"""
		self.FourierTransform.ForwardTransform(psi, self.TransformRank)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.FourierRepresentation)

	def TransformInverse(self, psi):
		"""
		Transforms psi from Fourier space to Grid space, and updates the representation
		"""
		self.FourierTransform.InverseTransform(psi, self.TransformRank)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.GridRepresentation)
	
	def CreateFourierRepr(self, repr):
		"""
		Creates a Fourier representation given a GridRepresentation
		"""
		fftrepr = self.FFTTransform.CreateFourierRepresentation(repr)
		return fftrepr

	def GetBasisFunction(self, rank, basisIndex):
		"""
		Return a fourier basis function
		"""
		fourierRepr = self.FourierRepresentation
		k = fourierRepr.GetGlobalGrid(rank)[basisIndex]
		x = self.psi.GetRepresentation().GetLocalGrid(rank)
		return exp( 1.0j * k * x )

	# Helper class for kinetic energy potentials
	class StaticEnergyConf(Section):
		def __init__(self, type, classname):
			self.type = type
			self.classname = classname
	
