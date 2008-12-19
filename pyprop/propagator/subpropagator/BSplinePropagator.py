class BSplinePropagator(SubPropagatorBase):
	__BASE = SubPropagatorBase

	def __init__(self, psi, transformRank):
		self.__BASE.__init__(self, psi, transformRank)

		self.Propagator = CreateInstanceRank("core.BSplinePropagator", psi.GetRank())
		self.TransformRank = transformRank

		self.Transform = CreateInstanceRank("core.BSplineTransform", psi.GetRank())

		self.RepresentationBSpline = psi.GetRepresentation().GetRepresentation(transformRank)
		self.BSplineObject = self.RepresentationBSpline.GetBSplineObject()


	def ApplyConfigSection(self, configSection):
		if hasattr(configSection, "propagation_algorithm"):
			if configSection.propagation_algorithm == 1:
				rank1repr = self.psi.GetRepresentation().GetRepresentation(1)
				if not isinstance(rank1repr, core.ReducedSphericalHarmonicRepresentation):
					raise Exception("Specified propagation_algorithm only works with ReducedSphericalRepresentation!")

		if hasattr(configSection, "centrifugal_potential"):
			if not hasattr(configSection, "angular_rank"):
				raise Exception("When centrifugal potential is present, angular rank must be specified!")

		configSection.Apply(self.Propagator)
		self.ConfigSection = configSection


	def SetupStep(self, dt):
	
		# Set up b-spline grid representation
		self.RepresentationGrid = core.BSplineGridRepresentation()
		self.RepresentationGrid.SetBaseRank(self.TransformRank)
		self.RepresentationGrid.SetupRepresentation(self.BSplineObject)
		self.RepresentationGrid.SetDistributedModel(self.RepresentationBSpline.GetDistributedModel())

		# Set up transform
		print "        Setup Forward Transform"
		self.Transform.SetupStep(self.psi, self.BSplineObject, self.TransformRank)
		print "        Setup Inverse Transform"
		self.InverseTransform(self.psi)


	def AdvanceStep(self, t, dt):
		self.Propagator.AdvanceStep(self.psi)
		self.InverseTransform(self.psi)


	def MultiplyHamiltonian(self, srcPsi, dstPsi, t, dt):
		self.Propagator.MultiplyHamiltonian(srcPsi, dstPsi)
		self.InverseTransform(srcPsi)
		self.InverseTransform(dstPsi)


	def SetupStepConjugate(self, dt):
		# Allocate arrays for b-spline potential and centrifugal potential
		globalGridSize = self.BSplineObject.GetQuadratureGridGlobal().size
		potentialVector = zeros(globalGridSize, dtype=double)

		# If present, set up centrifugal potential
		if hasattr(self.ConfigSection, "centrifugal_potential"):
			centrifugalPotentialName = self.ConfigSection.Get("centrifugal_potential")
			centrifugalPotentialSection = self.ConfigSection.Config.GetSection(centrifugalPotentialName)
			centrifugalPotential = CreatePotentialFromSection(centrifugalPotentialSection, centrifugalPotentialName, self.psi)
			centrifugalPotential.SetupStep(dt)
			centrifugalVector = zeros(globalGridSize, dtype=double)
			centrifugalVector[:] = centrifugalPotential.Evaluator.GetPotential(self.psi, dt, 0).real[:]
			self.Propagator.SetupCentrifugalPotential(centrifugalVector, self.psi)
		
		
		# If an additonal b-spline potential is present, we set it up, 
		# otherwise just kinetic potential
		if hasattr(self.ConfigSection, "potential"):
			potentialName = self.ConfigSection.Get("potential")
			potentialSection = self.ConfigSection.Config.GetSection(potentialName)
			potential = CreatePotentialFromSection(potentialSection, potentialName, self.psi)
			potential.SetupStep(dt)
			potentialVector[:] = potential.Evaluator.GetPotential(self.psi, dt, 0).real[:]

			# Set up propagator w/potential
			self.ForwardTransform(self.psi)
			self.Propagator.Setup(dt, self.psi, self.BSplineObject, potentialVector, self.TransformRank)

		else:

			# Set up propagator without potentials
			self.ForwardTransform(self.psi)
			self.Propagator.Setup(dt, self.psi, self.BSplineObject, self.TransformRank)

	
	def AdvanceStepConjugate(self, t, dt):
		self.ForwardTransform(self.psi)
		self.Propagator.AdvanceStep(self.psi)

	
	def MultiplyHamiltonianConjugate(self, srcPsi, dstPsi, t, dt):
		self.ForwardTransform(srcPsi)
		self.ForwardTransform(dstPsi)


	def ForwardTransform(self, psi):
		assert(not psi.GetRepresentation().GetDistributedModel().IsDistributedRank(self.TransformRank))
		self.Transform.ForwardTransform(psi)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationBSpline)

	def InverseTransform(self, psi):
		assert(not psi.GetRepresentation().GetDistributedModel().IsDistributedRank(self.TransformRank))
		self.Transform.InverseTransform(psi)
		psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationGrid)
