
#------------------------------------------------------------------------------------
#        Functions for converting tensor potentials to other formats
#------------------------------------------------------------------------------------

def SetupPotentialMatrix(prop, whichPotentials):
	"""
	sets up a "big-matrix", a dense matrix in normal storage format
	from some of the tensor potentials of a Problem object

	whichPotentials is a list of integers specifying the index
	of the potentials in the potential list of the prop object
	to include in the big matrix.
	"""

	print "Setting up potential matrix..."
	matrixSize = prop.psi.GetData().size
	
	#Allocate the potential matrix
	print "    Allocating potential matrix of size [%i, %i]  ~%.0f MB" \
		% (matrixSize, matrixSize, matrixSize**2 * 16 / 1024.**2)
	BigMatrix = zeros((matrixSize, matrixSize), dtype="complex")

	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)

		for idxL, idxR, i, j, k in TensorPotentialIndexMap(prop.psi.GetData().shape, potential):
			BigMatrix[idxL, idxR] += potential.PotentialData[i, j, k]

	return BigMatrix


def SetupBigMatrixReal(prop, whichPotentials):
	print "Setting up potential matrix..."
	matrixSize = prop.psi.GetData().size
	psiShape = prop.psi.GetData().shape
	
	#Allocate the hamilton matrix
	print "    Allocating potential matrix of size [%i, %i]  ~%.0f MB" \
		% (matrixSize, matrixSize, matrixSize**2 * 8 / 1024.**2)
	BigMatrix = zeros((matrixSize, matrixSize), dtype="double")

	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)

		basisPairs0 = potential.BasisPairs[0]
		basisPairs1 = potential.BasisPairs[1]
		basisPairs2 = potential.BasisPairs[2]
		
		Count0 = prop.psi.GetData().shape[0]
		Count1 = prop.psi.GetData().shape[1]
		Count2 = prop.psi.GetData().shape[2]

		tMap = TensorPotentialIndexMap
		for idxL, idxR, i, j, k in tMap(psiShape, potential):
			BigMatrix[idxL, idxR] += potential.PotentialData[i, j, k].real

	return BigMatrix




def SetupPotentialMatrixLL(prop, whichPotentials, eps=1e-14):
	"""
	Sets up a matrix in the LL format of pysparse from some of 
	the potentials of a Problem object
	"""

	def SetElement():
		MatrixLL[idxL, idxR] += potential.PotentialData[i, j, k].real

	print "Setting up potential matrix..."
	matrixSize = prop.psi.GetData().size

	#Find potential size for largest potential
	potList = prop.Propagator.BasePropagator.PotentialList
	potentialSize = max([potList[w].PotentialData.size for w in whichPotentials])

	#Set up linked list matrix
	MatrixLL = pysparse.spmatrix.ll_mat_sym(matrixSize, potentialSize)
	#MatrixLL = pysparse.spmatrix.ll_mat(matrixSize, matrixSize)

	countSize = 1e4

	for potNum in whichPotentials:
		potential = prop.Propagator.BasePropagator.PotentialList[potNum]
		print "    Processing potential %i: %s" % (potNum, potential.Name)
		curPotSize = potList[potNum].PotentialData.size

		count = 0
		outStr = ""
		for idxL, idxR, i, j, k in TensorPotentialIndexMap(prop.psi.GetData().shape, potential):
			#Skip upper triangle (we have a symmetric matrix)
			if idxL < idxR:
				continue

			#Skip current element if less than eps
			if abs(potential.PotentialData[i, j, k]) < eps:
				continue

			#Print progress info
			if mod(count, countSize) == 0:
				outStr = " " * 8
				outStr += "Progress: %i/%i" % (count/countSize, round(curPotSize/2./countSize))
				#outStr += "Progress: %i/%i" % (count/countSize, round(curPotSize/countSize))
				sys.stdout.write("\b"*len(outStr) + outStr)
				sys.stdout.flush()

			#Store element in linked-list matrix
			#MatrixLL[idxL, idxR] += potential.PotentialData[i, j, k].real
			SetElement()
			count += 1

		print


	return MatrixLL


def TensorPotentialIndexMap(psiShape, tensorPotential):
	"""
	Returns a generator for a map between indices in an m x m matrix and 
	"""

	basisPairs0 = tensorPotential.BasisPairs[0]
	basisPairs1 = tensorPotential.BasisPairs[1]
	basisPairs2 = tensorPotential.BasisPairs[2]

	basisCount0 = basisPairs0.shape[0]
	basisCount1 = basisPairs1.shape[0]
	basisCount2 = basisPairs2.shape[0]
	
	Count0 = psiShape[0]
	Count1 = psiShape[1]
	Count2 = psiShape[2]

	for i in xrange(basisCount0):
		xIndex0 = (basisPairs0[i,0] * Count1 * Count2)
		xIndex0p = (basisPairs0[i,1] * Count1 * Count2) 
		for j in xrange(basisCount1):
			xIndex1 = (basisPairs1[j,0] * Count2)
			xIndex1p = (basisPairs1[j,1] * Count2)
			for k in xrange(basisCount2):
				indexLeft = basisPairs2[k,0] + xIndex1 + xIndex0
				indexRight = basisPairs2[k,1] + xIndex1p + xIndex0p 
				yield indexLeft, indexRight, i, j, k



def TensorPotentialIndexMapOld(psiShape, tensorPotential):
	"""
	Returns a generator for a map between indices in an m x m matrix and 
	"""
	basisPairs0 = tensorPotential.BasisPairs[0]
	basisPairs1 = tensorPotential.BasisPairs[1]
	basisPairs2 = tensorPotential.BasisPairs[2]
	
	Count0 = psiShape[0]
	Count1 = psiShape[1]
	Count2 = psiShape[2]

	for i, (x0,x0p) in enumerate(basisPairs0):
		for j, (x1,x1p) in enumerate(basisPairs1):
			for k, (x2,x2p) in enumerate(basisPairs2):
				indexLeft = x2 + (x1 * Count2) + (x0 * Count1 * Count2) 
				indexRight = x2p + (x1p * Count2) + (x0p * Count1 * Count2) 
				yield indexLeft, indexRight, i, j, k


class BlockPreconditioner:
	def __init__(self, A, blockSize=1):
		self.Matrix = A
		self.shape = A.shape
		self.BlockSize = blockSize
		self.ProbSize = A.shape[0]
		self.NumCalls = 0
		
		#Check that blocksize divides matrix shape
		if mod(self.ProbSize, self.BlockSize) != 0:
			raise Exception("Preconditioner block size must divide matrix size!")

		#Set up preconditioner matrix
		self.Preconditioner = pysparse.spmatrix.ll_mat(self.shape[0], self.shape[1])
		curBlock = zeros((self.BlockSize, self.BlockSize))
		for k in range(self.ProbSize / self.BlockSize):
			curBlock[:] = 0.0
			curSlice = [slice(k*(blockSize), (k+1)*blockSize)] * 2
			
			for i in range(self.BlockSize):
				for j in range(self.BlockSize):
					I = k * blockSize + i
					J = k * blockSize + j
					curBlock[i,j] = self.Matrix[I, J]

			invBlock = linalg.inv(curBlock)

			for i in range(self.BlockSize):
				for j in range(self.BlockSize):
					I = k * blockSize + i
					J = k * blockSize + j
					self.Preconditioner[I,J] = invBlock[i,j]

	def precon(self, x, y):
		self.Preconditioner.matvec(x, y)
		self.NumCalls += 1



