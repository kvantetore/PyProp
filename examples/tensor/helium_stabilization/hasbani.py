execfile("example.py")


def SubmitStabilizationRun(workingDir):
	"""
	Calculate total ionization for a range of intensities to determine stabilization
	"""
	outputDir = "stabilization_freq_5/"
	frequency = 5.0
	amplitudeList = arange(2.0, 31.0)
	
	for I in amplitudeList:
		name = outputDir + "stabilization_I_%i.h5" % I
		Submit(executable="run_stabilization.py", \
			runHours=1, \
			jobname="stabilization", \
			numProcs=65, \
			config="config.ini", \
			amplitude=I/frequency, \
			outputCount=300, \
			workingDir=workingDir, \
			outputFilename=name, \
			findGroundstate = False, \
			writeScript=False)


def FindIonizationProbability(datafile, boundstateFiles, ionizationThreshhold=-2.0):
	"""
	Find total single and double ionization for the Hasbani example by
	projecting on states with energy < 2.0 a.u.
	"""

	conf = pyprop.Config(pyprop.serialization.GetConfigFromHDF5(datafile))
	lmax = conf.AngularRepresentation.index_iterator.lmax
	Lmax = conf.AngularRepresentation.index_iterator.L[-1]

	conf.Propagation.grid_potential_list = []
	conf.Propagation.preconditioner = None

	#h5file = tables.openFile(datafile)
	#try:
	#	ionizationProbability = h5file.root.Norm[0]
	#finally:
	#	h5file.close()
	ionizationProbability = 1.0
		
	#Set up problem
	#conf.AngularRepresentation.index_iterator = pyprop.DefaultCoupledIndexIterator(lmax=lmax, L=L)
	prop = pyprop.Problem(conf)
	tmpPsi = prop.psi.Copy()
	totalIdxIterator = pyprop.DefaultCoupledIndexIterator(lmax=lmax, L=range(Lmax))

	#Load wavefunction
	h5file = tables.openFile(datafile, "r")
	try:
		prop.psi.GetData()[:] = h5file.root.wavefunction[:]
	finally:
		h5file.close()
	for L in range(Lmax + 1):
		#Project on all bound states for current L
		#curIdx = [i for i in pyprop.DefaultCoupledIndexIterator(lmax=lmax, L=L)]
		h5file = tables.openFile(boundstateFiles.pop(0), "r")
		numEigs = size(h5file.root.Eig.Eigenvalues)
		for i in range(numEigs):
			tmpPsi.Clear()
			for j,cur in enumerate(totalIdxIterator):
				#if cur.L == L and h5file.root.Eig.Eigenvalues[i] < -1.99:
				if cur.L == L and h5file.root.Eig.Eigenvalues[i] < ionizationThreshhold:
					#print cur.l1, cur.l2, L, cur.L, h5file.root.Eig.Eigenvalues[i]
					tmpPsi.GetData()[j,:,:] += array(h5file.getNode("/Eig/Eigenvector%03i" % i))[cur.l1, :, :]
			ionizationProbability -= abs(prop.psi.InnerProduct(tmpPsi))**2

		h5file.close()

	return ionizationProbability

