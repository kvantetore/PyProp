
class ArpackSolver:

	def __init__(self, prop):
		self.BaseProblem = prop
		self.Rank = prop.psi.GetRank()

		#Create a copy of the wavefunction to calculate H|psi>
		self.TempPsi = prop.psi.CopyDeep()

		#Set up Arpack Solver
		self.Solver = CreateInstanceRank("core.krylov_ArpackPropagator", self.Rank)
		prop.Config.Arpack.Apply(self.Solver)
		self.Solver.Setup(prop.psi)

	def Solve(self):
		psi = self.BaseProblem.psi;
		tempPsi = self.TempPsi

		#Run the Arnoldi iterations
		self.Solver.Solve(self.__MatVecCallback, psi, tempPsi)

	def __MatVecCallback(self, psi, tempPsi):
		tempPsi.GetData()[:] = 0
		self.BaseProblem.Propagator.MultiplyHamiltonian(tempPsi, 0, 0)
	


