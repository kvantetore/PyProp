
def RunSphericalHarmonicDistribution(wavefunctionFile):
	rightPsi = pyprop.CreateWavefunctionFromFile(wavefunctionFile)
	leftPsi = rightPsi.Copy()

	#data = <left | right> = sum conj(left) * S * right
	repr = rightPsi.GetRepresentation()
	angRepr = repr.GetRepresentation(0)
	repr.MultiplyIntegrationWeights(rightPsi)
	data = real(conj(leftPsi.GetData()) * rightPsi.GetData())

	data = sum(sum(data, axis=2), axis=1)

	for i in range(pyprop.ProcCount):
		if i == pyprop.ProcId:
			print "Proc %i" % pyprop.ProcId
			for dataIndex, angIndex in enumerate(repr.GetLocalGrid(0)):
				coupledIndex = angRepr.Range.GetCoupledIndex(int(angIndex))
				print "%s = %s" % (coupledIndex, data[dataIndex])
		pyprop.pypar.barrier()

