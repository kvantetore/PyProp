
class GMRESShiftInvertSolver:
	"""
	A solver for using GMRES + some specified preconditioner to solve
	the linear system
	(S**-1 H - I * shift) psi_sol = psi_rhs                          (1)

	In orer to avold forming S**-1 H, this is rewritten as 
	(H - S * shift) psi_sol  = S psi_rhs,

	With left preconditioner this is
	M**-1 (H - S * shift) psi_sol  = M**-1 S psi_rhs

	which is the system we solve in this routine.

	It can be used by PiramSolver to speed up the search for eigenvalues
	significantly, especially if interior eigenvalues are sought

	In order for GMRES to work properly, it is almost alwyas nescessary to
	precondition it. The preconditioner is specified by the precondition config
	value in the GMRES config section, and should partially invert the matrix 
	(H - S*shift) 

	"""

	def __init__(self, prop):
		self.BaseProblem = prop
		self.Rank = prop.psi.GetRank()
		self.Config = self.BaseProblem.Config.GMRES
		self.psi = prop.psi
		self.TempPsi = prop.psi.Copy()
		self.TempPsi2 = prop.psi.Copy()
		self.Shift = self.Config.shift

		#Set up GMRES
		self.Solver = CreateInstanceRank("core.krylov_GmresWrapper", self.Rank)
		prop.Config.GMRES.Apply(self.Solver)
		self.Solver.Setup(self.psi)

		#Setup preconditioner
		config = prop.Config
		preconditionerName = self.Config.preconditioner
		if preconditionerName:
			preconditionerSection = config.GetSection(preconditionerName)
			preconditioner = preconditionerSection.type(self.psi)
			preconditionerSection.Apply(preconditioner)
			self.Preconditioner = preconditioner
			self.Preconditioner.SetHamiltonianScaling(1.0)
			self.Preconditioner.SetOverlapScaling(-self.Shift)
			self.Preconditioner.Setup(self.BaseProblem.Propagator)
		else:
			self.Preconditioner = None


	def __callback(self, srcPsi, dstPsi):
		"""
		This is the callback for the GMRES solver
		Here we do the following steps:
		
		   1) dstPsi = (H - shift * S) * srcPsi
		   2) dstPsi = M**-1 * dstPsi
		"""

		repr = dstPsi.GetRepresentation()
		self.TempPsi2.GetData()[:] = srcPsi.GetData()[:]
		
		#Multiply H * srcPsi, put in dstPsi (dstPsi is zeroed!)
		multiplyHam = self.BaseProblem.Propagator.BasePropagator.MultiplyHamiltonianNoOverlap
		multiplyHam(srcPsi, dstPsi, 0, 0)

		#Multiply shift * S * srcPsi
		repr.MultiplyOverlap(self.TempPsi2)
		self.TempPsi2.GetData()[:] *= -self.Shift
		dstPsi.GetData()[:] += self.TempPsi2.GetData()[:]

		#Solve left preconditioner in-place, dstPsi = M**-1 * dstPsi
		if self.Preconditioner:
			self.Preconditioner.Solve(dstPsi)


	def InverseIterations(self, srcPsi, destPsi, t, dt):
		"""
		InverseIterations is called from Piram as an alternative
		to MultiplyHamiltonian. On exit destPsi is set to 
		the solution of the linear system (1) with srcPsi as 
		right hand side
		"""

		repr = destPsi.GetRepresentation()

		#Multiply overlap, c = S * b
		self.TempPsi.GetData()[:] = srcPsi.GetData()[:]
		repr.MultiplyOverlap(self.TempPsi)

		#Preconditioner, d = M**-1 * c = M**-1 * S * b
		if self.Preconditioner:
			self.Preconditioner.Solve(self.TempPsi)

		#Solve for (H - shift * S)
		self.Solver.Solve(self.__callback, self.TempPsi, destPsi, False)


