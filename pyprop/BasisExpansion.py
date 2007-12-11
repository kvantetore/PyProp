
def OuterProduct(curFuncs, outData):
	rank = len(curFuncs)
	outerProductFunction = eval("core.OuterProduct_%s" % rank)
	outerProductFunction(curFuncs, outData)

def InnerProduct(a, b):
	return a.InnerProduct(b)

def IndexIterator(shape, rank=0):
	for localIndex in xrange(shape[rank]):
		if rank < len(shape) - 1:
			for j in IndexIterator(shape, rank+1):
				yield [localIndex] + j
		else:
			yield [localIndex]
		

class BasisFunctionGenerator(object):
	def __init__(self, prop):
		rank = prop.psi.GetRank()
		
		repr = prop.psi.GetRepresentation()
		grid = [repr.GetGlobalGrid(i) for i in xrange(rank)]
		self.BasisShape = array([len(grid[i]) for i in xrange(rank)])

		self.basisFun = []
		for j in xrange(rank):
			self.basisFun.append([array(prop.Propagator.GetBasisFunction(j, i), dtype=complex) for i in xrange(self.BasisShape[j])])

	def GetBasisSize(self):
		return self.BasisShape.prod()

	def Iterate(self, psi):
		basisCount = 0
		rank = len(self.BasisShape)
		for basisIndex in IndexIterator(self.BasisShape):
			curFuncs = [ self.basisFun[i][basisIndex[i]] for i in xrange(rank) ]
			OuterProduct(curFuncs, psi.GetData()[:])
			psi.Normalize()

			yield basisCount
			basisCount += 1


def GetBasisExpansionMatrix(prop):
	tempPsi = prop.GetTempPsi()
	psi = tempPsi.Copy()

	basisFuncs = BasisFunctionGenerator(prop)

	matrixSize = basisFuncs.GetBasisSize()
	matrix = zeros((matrixSize, matrixSize), dtype=complex)

	for j in basisFuncs.Iterate(prop.psi):
		#calculate |tempPsi> = H | prop.psi >
		tempPsi.GetData()[:] = 0 
		prop.MultiplyHamiltonian(tempPsi)

		for i in basisFuncs.Iterate(psi):
			#calculate < psi | H | prop.psi > = < psi | tempPsi >
			matrix[i, j] = InnerProduct(psi, tempPsi)

	return matrix

