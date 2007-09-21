
class CartesianRadialPropagator:
	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

		rank = psi.GetRank()	

		#instantiate transforms		
		self.FFTTransform = CreateInstanceRank("core.CartesianFourierTransform", 1)
		self.RadialTransform = CreateInstanceRank("core.RadialTransform", rank)
		self.RadialKineticPotential = None

	def ApplyConfigSection(self, configSection):
		self.Mass = configSection.mass

		self.ForceOriginZero = True
		if hasattr(configSection, "force_origin_zero"):
			self.ForceOriginZero = configSection.force_origin_zero

	def SetupStep(self, dt):
		#setup representation
		self.RadialRepresentation = self.psi.GetRepresentation().GetRepresentation(self.TransformRank)
		self.RadialFourierRepresentation = self.CreateRadialFourierRepr(self.RadialRepresentation)

		# transform radial into fourier space
		print "	       Setup Forward Radial Transform"
		self.TransformRadialForward(self.psi)
		
		# set up potential
		print "	       Setup Radial Kinetic Potential"
		self.SetupRadialKineticPotential(dt)

		# transform radial into grid space
		print "	       Setup Inverse Radial Transform"
		self.TransformRadialInverse(self.psi)

		#Get index of origin
		if self.ForceOriginZero:
			self.OriginIndex = self.GetOriginIndex()

	def AdvanceStep(self, t, dt):
		# transform the radial into fourier space
		self.TransformRadialForward(self.psi)
	
		# apply radial kinetic energy potential
		self.RadialKineticPotential.AdvanceStep(t, dt)
		
		# transform back into real space
		self.TransformRadialInverse(self.psi)

		#For imaginary time propagation, set origin to 0 to keep symmetry
		if 1.0j * imag(dt) == dt and self.ForceOriginZero:
			self.SetValue(self.OriginIndex, 0)

	def MultiplyHamiltonian(self, dstPsi, t, dt):
		# Calculates dstPsi += - 1/(2m) d^2/dx^2 psi
		self.TransformRadialForward(self.psi)
		self.TransformRadialForward(dstPsi)
		self.RadialKineticPotential.MultiplyPotential(dstPsi, t, dt)
		self.TransformRadialInverse(dstPsi)
		self.TransformRadialInverse(self.psi)
	
	def SetupStepConjugate(self, dt):
		pass

	def AdvanceStepConjugate(self, t, dt):
		self.AdvanceStep(t, dt)

	def MultiplyHamiltonianConjugate(self, dstPsi, t, dt):
		pass
			
	def GetOriginIndex(self):
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

	def TransformRadialForward(self, psi):
		self.RadialTransform.ForwardTransform(psi, self.TransformRank)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RadialFourierRepresentation)

	def TransformRadialInverse(self, psi):
		self.RadialTransform.InverseTransform(psi, self.TransformRank)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RadialRepresentation)
	
	# Helper class for kinetic energy potentials
	class StaticEnergyConf(Section):
		def __init__(self, type, classname):
			self.type = type
			self.classname = classname
	
	def SetupRadialKineticPotential(self, dt):
		radialConf = self.StaticEnergyConf(PotentialType.Static, "core.RadialKineticEnergyPotential")
		radialConf.mass = self.Mass
		radialConf.radial_rank = self.TransformRank
		radialPot = CreatePotentialFromSection(radialConf, "RadialKineticEnergy", self.psi)
		radialPot.SetupStep(dt)
		self.RadialKineticPotential = radialPot
	
	def CreateRadialFourierRepr(self, repr):
		fftrepr = self.FFTTransform.CreateFourierRepresentation(repr)
		return fftrepr


