
def GetRadialMatricesCompressedCol(pot, psi):
	"""
	Sets up a list of radial matrices in compressed column 
	format from a tensor potential.

	The matrix is assumed to be diagonal in the angular rank, and one
	radial matrix is set up for each index in the angular rank.
	"""

	localShape = psi.GetData().shape

	#do parallelizaion on angular rank
	repr = psi.GetRepresentation()
	if repr.GetDistributedModel().IsDistributedRank(1) or repr.GetDistributedModel().IsDistributedRank(2):
		raise Exception("Only angular rank can be distributed")

	#Check that the angular basis pairs is diagonal
	for row, col in zip(pot.BasisPairs[0][:,0], pot.BasisPairs[0][:,1]):
		if row != col:
			raise Exception("%i != %i, angular rank must be diagonal for potential %s" % (row, col, pot.Name))

	#stats for each radial matrix
	matrixSize = localShape[1] * localShape[2]
	nonzeroCount = pot.BasisPairs[1].shape[0] * pot.BasisPairs[2].shape[0]

	#Create a list of radial index pairs in the full matrix
	matrixPairs = zeros((nonzeroCount, 2), dtype=int)
	matrixIndex = 0
	for i, (row0, col0) in enumerate(pot.BasisPairs[1]):
		rowIndex0 = (row0 * localShape[2])
		colIndex0 = (col0 * localShape[2]) 
		for j, (row1, col1) in enumerate(pot.BasisPairs[2]):
			rowIndex = rowIndex0 + row1
			colIndex = colIndex0 + col1

			matrixPairs[matrixIndex, 0] = rowIndex
			matrixPairs[matrixIndex, 1] = colIndex
			matrixIndex += 1

	#Sort indexpairs by column
	sortIndices = argsort(matrixPairs[:,1])

	#create colStart- and row-indices required by the col-compressed format
	# colStartIndices is the indices in the compressed matrix where a col starts
	# rowIndices is the row of every element in the compressed matrix
	colStartIndices = zeros(matrixSize+1, dtype=int32)
	rowIndices = zeros(nonzeroCount, dtype=int32)
	prevCol = -1
	for i, sortIdx in enumerate(sortIndices):
		row, col = matrixPairs[sortIdx, :]
		rowIndices[i] = row

		if prevCol != col:
			if prevCol != col -1:
				raise Exception("What? Not a col-compressed format?")
			colStartIndices[col] = i
			prevCol = col
	colStartIndices[-1] = nonzeroCount

	#Create a radial matrix in compressed col format for each angular index
	radialMatrices = []
	localAngularCount = localShape[0]
	for curAngularIndex in range(localAngularCount):
		potentialSlice = pot.PotentialData[curAngularIndex, :, :].flat
		radialMatrix = zeros(nonzeroCount, dtype=complex)
		
		for i, sortIdx in enumerate(sortIndices):
			radialMatrix[i] = potentialSlice[sortIdx]

		radialMatrices.append( radialMatrix )

	return rowIndices, colStartIndices, radialMatrices

				
class RadialTwoElectronPreconditioner:
	"""
	Preconditioner for GMRES for solving systems of the type

	(1)     (a H + b S) x = y  ->  x

	where H is the Hamiltonian, S is the overlap matrix
	and a and b are complex scaling factors (scalingH and scalingS).

	The preconditioner is given a set of tensor potentials which are 
	diagonal in the angular rank (also called radial potentials). 
	These potentials should approximate the Hamiltonian as well as 
	possible, as the system (1) is solved exactly for the given radial
	potentials
	"""

	def __init__(self, psi):
		self.Rank = psi.GetRank()
		self.psi = psi

	def ApplyConfigSection(self, conf):
		self.OverlapSection = conf.Config.GetSection("OverlapPotential")
		self.PotentialSections = [conf.Config.GetSection(s) for s in conf.potential_evaluation]

	def SetHamiltonianScaling(self, scalingH):
		self.HamiltonianScaling = scalingH

	def SetOverlapScaling(self, scalingS):
		self.OverlapScaling = scalingS

	def GetHamiltonianScaling(self):
		return self.HamiltonianScaling

	def GetOverlapScaling(self):
		return self.OverlapScaling

	def Setup(self, prop):
		"""
		Set up a tensor potential for overlap potential and all other potentials
		and add them together, assuming they have the same layout
		ending up with a potential containing S + scalingH * (P1 + P2 + ...)

		The radial part of this potential is then converted to compressed col storage
		and factorized.
		"""

		#Setup overlap potential
		tensorPotential = prop.BasePropagator.GeneratePotential(self.OverlapSection)
		tensorPotential.PotentialData *= self.GetOverlapScaling()

		#Add all potentials to solver
		scalingH = self.GetHamiltonianScaling()
		for conf in self.PotentialSections:
			#Setup potential in basis
			potential = prop.BasePropagator.GeneratePotential(conf)
			if not tensorPotential.CanConsolidate(potential):
				raise Exception("Cannot consolidate potential %s with overlap-potential" % (potential.Name))
		
			#Add potential
			potential.PotentialData *= scalingH
			tensorPotential.PotentialData += potential.PotentialData
			del potential

	
		#Setup radial matrices in CSC format
		tensorPotential.SetupStep(0.0)
		row, colStart, radialMatrices = GetRadialMatricesCompressedCol(tensorPotential, self.psi)

		shape = self.psi.GetRepresentation().GetFullShape()
		matrixSize = shape[1] * shape[2]

		#factorize each matrix
		radialSolvers = []
		for mat in radialMatrices:
			#M = scipy.sparse.csc_matrix((mat, row, colStart))
			#solve = scipy.linsolve.factorized(M)
			#radialSolvers.append(solve)

			solve = SuperLUSolver_2()
			solve.Setup(int(matrixSize), mat, row, colStart)
			radialSolvers.append(solve)

		self.RadialSolvers = radialSolvers

	def Solve(self, psi):
		data = psi.GetData()
		
		angularCount = data.shape[0]
		if angularCount != len(self.RadialSolvers):
			raise Exception("Invalid Angular Count")

		for angularIndex, solve in enumerate(self.RadialSolvers):
			#data[angularIndex,:,:].flat[:] = solve(data[angularIndex,:,:].flatten())
			solve.Solve(data[angularIndex, :, :])
		


