#---------------------------------------------------------------------
# SphericalHarmonic Propagator
#---------------------------------------------------------------------

class FullSphericalPropagator:
	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

		if transformRank != psi.GetRank() - 1:
			raise "SphericalTransform can only be used on the last rank"

		#Create spherical transform
		self.Rank = psi.GetRank()
		self.SphericalTransform = CreateInstanceRank("core.SphericalTransform", self.Rank);


	def ApplyConfigSection(self, config):
		self.Mass = config.mass
		self.RadialRank = config.radial_rank


	#Setup
	def SetupStep(self, dt):
		# Set up Transform
		self.SphericalTransform.SetupStep(self.psi);
		
		# Set up Representations 
		print "        Setup representations"
		self.LmRepresentation = self.psi.GetRepresentation().GetAngularRepresentation()
		self.AngularRepresentation = self.SphericalTransform.CreateAngularRepresentation()
		self.AngularRepresentation.SetDistributedModel(self.LmRepresentation.GetDistributedModel())

		# Set up angular kinetic energy potential
		print "        Setup Angular Kinetic Potential"
		self.SetupAngularKineticPotential(dt)

		#transform into grid space
		print "	       Setup Inverse Spherical Transform"
		self.TransformSphericalInverse()
	
	def SetupStepConjugate(self, dt):
		#transform into spherical harmonic space
		print "	       Setup Inverse Spherical Transform"
		self.TransformSphericalForward()
		#Angular kinetic energy potential is already set up

	
	#Advance
	def AdvanceStep(self, t, dt):
		#apply angular kinetic potential
		self.AngularKineticPotential.AdvanceStep(t, dt)
		# transform back into spherical (theta,phi) space
		self.TransformSphericalInverse()

	def AdvanceStepConjugate(self, t, dt):
		# transform into spherical harmonic (l,m) space
		self.TransformSphericalForward()
		#apply potential
		self.AngularKineticPotential.AdvanceStep(t, dt)
			

	#Transforms:		
	def TransformSphericalForward(self):
		self.SphericalTransform.ForwardTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.LmRepresentation)

	def TransformSphericalInverse(self):
		self.SphericalTransform.InverseTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.AngularRepresentation)
		
	# Helper class for kinetic energy potentials
	class StaticEnergyConf(Section):
		def __init__(self, type, classname):
			self.type = type
			self.classname = classname
	
	def SetupAngularKineticPotential(self, dt):
		angularConf = self.StaticEnergyConf(PotentialType.Static, "core.AngularKineticEnergyPotential")
		angularConf.mass = self.Mass
		angularConf.radial_rank = self.RadialRank 
		angularPot = CreatePotentialFromSection(angularConf, "AngularKineticEnergy", self.psi)
		angularPot.SetupStep(dt)
		self.AngularKineticPotential = angularPot



