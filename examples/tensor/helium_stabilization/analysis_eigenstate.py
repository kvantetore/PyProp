
#-----------------------------------------------------------------------------
#             tools for projecting out bound states
#-----------------------------------------------------------------------------

def GetEigenstateFileConfig(filename):
	f = tables.openFile(filename, "r")
	try:
		#load config from file
		conf = pyprop.Config(f.root.Eig._v_attrs.configObject)

		return conf
	except:
		f.close()


INFO_L = "L"
INFO_lmax = "lmax"
INFO_RadialGrid = "radial_grid"

def GetConfigInfo(conf, infoId):
	if infoId == INFO_L:
		indexIterator = conf.AngularRepresentation.index_iterator
		Llist = array([idx.L for idx in indexIterator])
		L = Llist[0]
		if any(Llist != L):
			raise Exception("L is not unique for the given config (%s)" % Llist)
		return L

	elif infoId == INFO_lmax:
		indexIterator = conf.AngularRepresentation.index_iterator
		lmax = max(array([idx.l1 for idx in indexIterator]))
		return lmax

	elif infoId == INFO_RadialGrid:
		return GetRadialGrid(config=conf)

	else:
		raise Exception("Unknown infoId %s" % infoId)


def GetEigenstateFileInfo(filename, infoId):
	conf = GetEigenstateFileConfig(filename)
	return GetConfigInfo(conf, infoId)


def RemoveBoundStateProjection(psi, boundstateFiles):
	assert(pyprop.IsSingleProc())

	projectionList = []
	eigenvaluesList = []

	for curFilename in boundstateFiles:
		L = GetEigenstateFileInfo(curFileName, INFO_L)

		#Get the local indices corresponding to the local L
		LFilter = lambda idx: idx.L == L
		indexL = GetLocalCoupledSphericalHarmonicIndices(psi, LFilter)
	
		eigPsi = pyprop.CreateWavefunctionFromFile(eigenstateFile, GetEigenvectorDatasetPath(0))
		projPsi = eigPsi.Copy()
		projPsi.GetData()[:] = psi.GetData()[indexL, :, :]

		f = tables.openFile(curFilename, "r")
		try:
			eigenvalues = f.root.Eig.Eigenvalues[:]
		finally:
			f.close()

		curProjList = []

		for i, E in enumerate(eigenvalues):
			#load eigenstate
			pyprop.serialization.LoadWavefunctionHDF(eigenstateFile, GetEigenvectorDatasetPath(i), eigPsi)
			#calculate projection
			proj = eigPsi.InnerProduct(projPsi)
			curProjList.append(proj)
			#remove projection
			psi.GetData()[indexL, :, :] -= proj * eigPsi.GetData()

		projectionList.append(curProjList)
		eigenvaluesList.append(eigenvalues)

	return projectionList, eigenvaluesList


def RunRemoveBoundStateProjection(wavefunctionFile, boundstateFiles):
	psi = pyprop.CreateWavefunctionFromFile(wavefunctionFile)
	psi.Normalize()
	projList, evList = RemoveBoundStateProjection(psi, boundstateFiles)
	print "Probability not in projected    = %s" % real(psi.InnerProduct(psi))
	print "Probability in projected states = %s" % (sum([sum(abs(array(p2))**2) for p1 in projList]))


#-----------------------------------------------------------------------------
#             tools for calculating and plotting two-particle dP/dE
#-----------------------------------------------------------------------------

def GetEigenstateProjection(psi, eigenstateFile, eigenstateL):
	"""
	Project psi on all eigenstates in the file eigenstateFile, and 
	return the energies as well as projection coefficients
	"""
	assert(pyprop.IsSingleProc())

	f = tables.openFile(eigenstateFile, "r")
	try:
		E = f.root.Eig.Eigenvalues[:]
	finally:
		f.close()

	eigPsi = pyprop.CreateWavefunctionFromFile(eigenstateFile, GetEigenvectorDatasetPath(0))

	projPsi = eigPsi.Copy()

	#Get the local indices corresponding to the local L
	LFilter = lambda idx: idx.L == eigenstateL
	angularIndices = GetLocalCoupledSphericalHarmonicIndices(psi, LFilter)
	print "Psi has %i L=%i" % (len(indexL), eigenstateL)
	print "ProjPsi has %i angular states" % (projPsi.GetData().shape[angularRank])
	projPsi.GetData()[:] = psi.GetData()[indexL, :, :]
	
	proj = []
	for i, curE in enumerate(E):
		infoStr =  "Progress: %3i%%" % ((i * 100)/ len(E))
		sys.stdout.write("\b"*15 + infoStr)
		sys.stdout.flush()

		pyprop.serialization.LoadWavefunctionHDF(eigenstateFile, GetEigenvectorDatasetPath(i), eigPsi)
		proj.append( eigPsi.InnerProduct(projPsi) )

	return E, array(proj)


def RunEigenstateProjection(wavefunctionFile, eigenstateFile, outputFile, L):
	psi = pyprop.CreateWavefunctionFromFile(wavefunctionFile)
	E, proj = GetEigenstateProjection(psi, eigenstateFile, L)

	f = tables.openFile(outputFile, "w")
	try:
		group = f.createGroup(f.root, "L%02i" % L)
		f.createArray(group, "Eigenvalues", E)
		f.createArray(group, "EigenstateProjection", proj)

	finally:
		f.close()


def LoadEigenstateProjection(projectionFile, L):
	f = tables.openFile(projectionFile, "r")
	try:
		E = f.getNode(f.root, "/L%02i/Eigenvalues" % L)[:]
		proj = f.getNode(f.root, "/L%02i/EigenstateProjection" % L)[:]

	finally:
		f.close()

	return E, proj	
	

def CalculateEnergyDistribution(E, projection, bins=50):
	minE = min(E)
	maxE = max(E)

	binE, dE = linspace(minE, maxE, bins, endpoint=True, retstep=True)
	binProj = zeros(len(binE), dtype=double)

	curBin = 0
	for E, proj in zip(E, projection):
		while E > binE[curBin+1]:
			curBin += 1
		binProj[curBin] += abs(proj)**2

	return binE, binProj/dE
		
def PlotEnergyDistribution(bins=50):
	for L in [0,1,2,3,4,5]:
		E, projection = LoadEigenstateProjection("projection_L%02i.h5" % L, L)
		binE, binProj = CalculateEnergyDistribution(E, projection, bins)
		semilogy(binE, binProj, label="L = %i" % L)

	legend()
