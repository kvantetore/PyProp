class BSplinePropagator:
	def __init__(self, psi, transformRank):
		self.psi = psi
		self.TransformRank = transformRank

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
					raise Exception("Specified projection algorithm only works with ReducedSphericalRepresentation!")

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
		self.Transform.InverseTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationGrid)


	def AdvanceStep(self, t, dt):
		self.Propagator.AdvanceStep(self.psi)
		self.InverseTransform()


	def MultiplyHamiltonian(self, dstPsi, t, dt):
		self.Propagator.MultiplyHamiltonian(self.psi, dstPsi)
		#self.Potential.MultiplyPotential(dstPsi, t, dt)
		self.Transform.InverseTransform(self.psi)
		self.Transform.InverseTransform(dstPsi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationGrid)
		dstPsi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationGrid)


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
			self.Transform.ForwardTransform(self.psi)
			self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationBSpline)
			assert(self.RepresentationBSpline == self.psi.GetRepresentation().GetRepresentation(self.TransformRank))
			self.Propagator.Setup(dt, self.psi, self.BSplineObject, potentialVector, self.TransformRank)

		else:

			# Set up propagator without potentials
			self.Transform.ForwardTransform(self.psi)
			self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationBSpline)
			assert(self.RepresentationBSpline == self.psi.GetRepresentation().GetRepresentation(self.TransformRank))
			self.Propagator.Setup(dt, self.psi, self.BSplineObject, self.TransformRank)

	
	def AdvanceStepConjugate(self, t, dt):
		self.ForwardTransform()
		self.Propagator.AdvanceStep(self.psi)

	
	def MultiplyHamiltonianConjugate(self, dstPsi, t, dt):
		self.Transform.ForwardTransform(self.psi)
		self.Transform.ForwardTransform(dstPsi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationBSpline)
		dstPsi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationBSpline)


	def ForwardTransform(self):
		self.Transform.ForwardTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationBSpline)

	def InverseTransform(self):
		self.Transform.InverseTransform(self.psi)
		self.psi.GetRepresentation().SetRepresentation(self.TransformRank, self.RepresentationGrid)
