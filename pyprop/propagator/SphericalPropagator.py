#---------------------------------------------------------------------
# SphericalHarmonic Momentum Evaluator
#---------------------------------------------------------------------

class RadialPropagator:
	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

	def SetupStep(self, dt): pass
	def AdvanceStep(self, t, dt): pass

class TransformedRadialPropagator(RadialPropagator):
	__Base = RadialPropagator
	
	def __init__(self, psi, transformRank):
		self.__Base.__init__(self, psi, transformRank)

		rank = psi.GetRank()	
		self.Propagator = CreateInstanceRank("core.TransformedGridPropagator", rank)
		self.TransformRank = transformRank

	def ApplyConfigSection(self, configSection):
		configSection.Apply(self.Propagator)

	def SetupStep(self, dt):
		param = self.psi.GetRepresentation().GetRepresentation(self.TransformRank).Range.Param
		
		self.Propagator.Setup(param, dt, self.psi, self.TransformRank)

	def AdvanceStep(self, t, dt):
		self.Propagator.AdvanceStep(self.psi)
	

class CartesianRadialPropagator(RadialPropagator):
	__Base = RadialPropagator
	
	def __init__(self, psi, transformRank):
		self.__Base.__init__(self, psi, transformRank)

		rank = psi.GetRank()	

		#instantiate transforms		
		self.FFTTransform = CreateInstanceRank("core.CartesianFourierTransform", 1)
		self.RadialTransform = CreateInstanceRank("core.RadialTransform", rank)
		self.RadialKineticPotential = None

	def ApplyConfigSection(self, configSection):
		self.Mass = configSection.mass

	def SetupStep(self, dt):
		#setup representation
		self.RadialRepresentation = self.psi.GetRepresentation().GetRepresentation(self.TransformRank)
		self.RadialFourierRepresentation = self.CreateRadialFourierRepr(self.RadialRepresentation)

		# transform radial into fourier space
		print "	       Setup Forward Radial Transform"
		self.TransformRadialForward()
		
		# set up potential
		print "	       Setup Radial Kinetic Potential"
		self.SetupRadialKineticPotential(dt)

		# transform radial into grid space
		print "	       Setup Inverse Radial Transform"
		self.TransformRadialInverse()

		#Get index of origin
		self.OriginIndex = self.GetOriginIndex()

	def AdvanceStep(self, t, dt):
		# transform the radial into fourier space
		self.TransformRadialForward()
	
		# apply radial kinetic energy potential
		self.RadialKineticPotential.AdvanceStep(t, dt)
		
		# transform back into real space
		self.TransformRadialInverse()

		#For imaginary time propagation, set origin to 0 to keep symmetry
		if 1.0j * imag(dt) == dt and self.OriginIndex != None:
			self.SetValue(self.OriginIndex, 0)
			
			
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

	def TransformRadialForward(self):
		self.RadialTransform.ForwardTransform(self.psi, self.TransformRank)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RadialFourierRepresentation)

	def TransformRadialInverse(self):
		self.RadialTransform.InverseTransform(self.psi, self.TransformRank)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RadialRepresentation)
	
	# Helper class for kinetic energy potentials
	class StaticEnergyConf(Section):
		def __init__(self, type, classname):
			self.type = type
			self.classname = classname
	
	def SetupRadialKineticPotential(self, dt):
		radialConf = self.StaticEnergyConf(PotentialType.Static, "core.DiagonalRadialPotential")
		radialConf.mass = self.Mass
		radialConf.radial_rank = self.TransformRank
		radialPot = CreatePotentialFromSection(radialConf, "RadialKineticEnergy", self.psi)
		radialPot.SetupStep(dt)
		self.RadialKineticPotential = radialPot
	
	def CreateRadialFourierRepr(self, repr):
		fftrepr = self.FFTTransform.CreateFourierRepresentation(repr)
		return fftrepr


class SphericalPropagator(PropagatorBase):
	__Base = PropagatorBase
	
	def __init__(self, psi):
		self.__Base.__init__(self, psi)

		#Create spherical transform
		self.Rank = psi.GetRank()
		self.SphericalTransform = CreateInstanceRank("core.SphericalTransform", self.Rank);

		#instantiate potentials
		self.AngularKineticPotential = None 

	def ApplyConfig(self, config):
		self.__Base.ApplyConfig(self, config)
	
		#Create all but spherical propagators
		self.SubPropagators = []
		for i in range(self.Rank-1):
			sectionName = config.Propagation.Get("propagator" + str(i))
			section = config.GetSection(sectionName)

			prop = section.propagator(self.psi, i)
			print "Propagator for rank %i is %s" % (i, prop)
			
			section.Apply(prop)

			self.SubPropagators.append(prop)
			
	def ApplyConfigSection(self, configSection): 
		self.__Base.ApplyConfigSection(self, configSection)

		self.AngularMass = configSection.mass
	
	def SetupStep(self, dt):
		print "        Setup SphericalTransform"
		self.SphericalTransform.SetupStep(self.psi);
	
		#get representations
		print "        Setup representations"
		self.LmRepresentation = self.psi.GetRepresentation().GetAngularRepresentation()
		self.AngularRepresentation = self.SphericalTransform.CreateAngularRepresentation()
		self.AngularRepresentation.SetDistributedModel(self.LmRepresentation.GetDistributedModel())

		#transform into grid space
		print "	       Setup Inverse Spherical Transform"
		self.TransformSphericalInverse()
		
		#setup potential
		print "        Setup Potential"
		self.SetupPotential(dt)
		
		# transform into spherical harmonics 
		print "        Setup Forward Spherical Transform"
		self.TransformSphericalForward()
	
		# set up potential
		print "        Setup Angular Kinetic Potential"
		self.SetupAngularKineticPotential(dt)

		# setup radial propagator
		print "	       Setup Radial Propagator"
		for prop in self.SubPropagators:
			prop.SetupStep(dt)
		
	
	def AdvanceStep(self, t, dt):
		# transform back into spherical (theta,phi) space
		self.TransformSphericalInverse()

		#apply grid potential
		self.ApplyPotential(t, dt)
	
		# transform into spherical harmonic (l,m) space
		self.TransformSphericalForward()
	
		#apply potential
		self.AngularKineticPotential.AdvanceStep(t, dt)
		
		#radial propagator
		for prop in self.SubPropagators:
			prop.AdvanceStep(t, dt)

	#Transforms:		
	def TransformSphericalForward(self):
		self.SphericalTransform.ForwardTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.Rank-1, self.LmRepresentation)

	def TransformSphericalInverse(self):
		self.SphericalTransform.InverseTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.Rank-1, self.AngularRepresentation)
		
	# Helper class for kinetic energy potentials
	class StaticEnergyConf(Section):
		def __init__(self, type, classname):
			self.type = type
			self.classname = classname
	
	def SetupAngularKineticPotential(self, dt):
		angularConf = self.StaticEnergyConf(PotentialType.Static, "core.DiagonalAngularPotential")
		angularConf.mass = self.AngularMass
		angularConf.radial_rank = self.Rank-2 #the last non-spherical rank
		print "using mass ", self.AngularMass
		angularPot = CreatePotentialFromSection(angularConf, "AngularKineticEnergy", self.psi)
		print "Create potential..."
		angularPot.SetupStep(dt)
		print "done..."
		self.AngularKineticPotential = angularPot

	def IsDistributedRank(self, rank):
		return self.psi.GetRepresentation().GetDistributedModel().IsDistributedRank(self.psi, rank)
		

