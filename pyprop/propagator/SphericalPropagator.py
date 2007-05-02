#---------------------------------------------------------------------
# SphericalHarmonic Propagator
#---------------------------------------------------------------------

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

		# Transpose
		self.SetupTranspose()
		self.Transpose(1)

		# setup radial propagator
		print "	       Setup Radial Propagator"
		for prop in self.SubPropagators:
			prop.SetupStep(dt)

		#transpose back
		self.Transpose(2)

	def AdvanceStep(self, t, dt):
		# transform back into spherical (theta,phi) space
		self.TransformSphericalInverse()

		#apply grid potential
		self.ApplyPotential(t, dt)
	
		# transform into spherical harmonic (l,m) space
		self.TransformSphericalForward()
	
		#apply potential
		self.AngularKineticPotential.AdvanceStep(t, dt)
	
		#transpose	
		self.Transpose(1)
		
		#radial propagator
		for prop in self.SubPropagators:
			prop.AdvanceStep(t, dt)

		#transpose back
		self.Transpose(2)

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
			self.Distribution2 = array([self.Rank-1])
			distribShape = transpose.CreateDistributedShape(fullShape, self.Distribution2)
			#allocate wavefunction
			self.TransposeBuffer1 = self.psi.GetActiveBufferName()
			self.TransposeBuffer2 = self.psi.AllocateData(distribShape)

	def Transpose(self, stage):
		distrModel = self.psi.GetRepresentation().GetDistributedModel()
		if not distrModel.IsSingleProc():
			if stage == 1:
				distrModel.ChangeDistribution(self.psi, self.Distribution2, self.TransposeBuffer2)
			elif stage == 2:
				distrModel.ChangeDistribution(self.psi, self.Distribution1, self.TransposeBuffer1)
			else:
				raise "Invalid stage %i" % stage

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
		angularConf = self.StaticEnergyConf(PotentialType.Static, "core.AngularKineticEnergyPotential")
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
		

