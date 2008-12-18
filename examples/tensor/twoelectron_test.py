
def SetupBigMatrix(prop, whichPotential):
	matrixSize = prop.psi.GetData().size
	potential = prop.Propagator.BasePropagator.PotentialList[whichPotential]

	BigMatrix = zeros((matrixSize, matrixSize), dtype="complex")

	basisPairs1 = potential.BasisPairs[0]
	basisPairs2 = potential.BasisPairs[1]
	basisPairs3 = potential.BasisPairs[2]
	
	GetCoupledIndex = prop.psi.GetRepresentation().GetRepresentation(2).Range.GetCoupledIndex
	r1Count = prop.psi.GetData().shape[0]
	r2Count = prop.psi.GetData().shape[1]
	lCount = prop.psi.GetData().shape[2]

	for i, (A, Ap) in enumerate(basisPairs3):
		for j, (r1,r1p) in enumerate(basisPairs1):
			for k, (r2,r2p) in enumerate(basisPairs2):
				#indexLeft = (A*r1Count*r2Count) + (r1 * r2Count) + r2
				#indexRight = (Ap*r1Count*r2Count) + (r1p * r2Count) + r2p
				indexLeft = (r1 * r2Count*lCount) + (r2*lCount) + A
				indexRight = (r1p * r2Count*lCount) + (r2p*lCount) + Ap
				BigMatrix[indexLeft, indexRight] = potential.PotentialData[j, k, i]


	return BigMatrix

