
#-----------------------------------------------------------------------------
#             tools for projecting out bound states
#-----------------------------------------------------------------------------

def GetEigenstateFileConfig(filename):
	f = tables.openFile(filename, "r")
	try:
		#load config from file
		conf = pyprop.Config(f.root.Eig._v_attrs.configObject)

		return conf
	finally:
		f.close()


INFO_L = "L"
INFO_lmax = "lmax"
INFO_RadialGrid = "radial_grid"
INFO_Eigenvalues = "eigenvalues"

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
	if infoId == "eigenvalues":
		#Get energies
		f = tables.openFile(filename, "r")
		try:
			eigenvalues = f.root.Eig.Eigenvalues[:]
		finally:
			f.close()
		return eigenvalues

	else:
		conf = GetEigenstateFileConfig(filename)
		return GetConfigInfo(conf, infoId)


def GetBoundStateFiles(**args):
	locations = ["output/boundstates", "output/boundstates_mine"]
	conf = SetupConfig(**args)
	lmax = max([l1 for l1, l2, L, M in conf.AngularRepresentation.index_iterator])
	Llist = unique([L for l1, l2, L, M in conf.AngularRepresentation.index_iterator])
	getPostfix = lambda L: "_".join(GetRadialGridPostfix(config=conf, lmax=lmax, L=L) + GetAngularGridPostfix(config=conf, lmax=lmax, L=L)) 
	#getFilename = lambda L: "output/boundstates/boundstates_%s.h5" % (getPostfix(L))
	getFilename = lambda L: filter(os.path.exists, ["%s/boundstates_%s.h5" % (loc, getPostfix(L)) for loc in locations])[0]

	return map(getFilename, Llist)


def GetBoundStates(ionizationThreshold=-2.0, **args):
	boundFiles = GetBoundStateFiles(**args)
	print boundFiles

	def loadStates(filename):
		L = GetEigenstateFileInfo(filename, INFO_L)
		eigPsi = pyprop.CreateWavefunctionFromFile(filename, GetEigenvectorDatasetPath(0))

		curPsiList = []
		curEnergyList = []

		eigenvalues = GetEigenstateFileInfo(filename, INFO_Eigenvalues)

		boundIdx = filter(lambda i: eigenvalues[i]<ionizationThreshold, r_[:len(eigenvalues)])
		for i in boundIdx:
			psi = eigPsi.Copy()
			pyprop.serialization.LoadWavefunctionHDF(filename, GetEigenvectorDatasetPath(i), psi)
			curPsiList.append(psi)
			curEnergyList.append(eigenvalues[i])

		return curEnergyList, (curPsiList, L)

	energies, boundStates = zip(*map(loadStates, boundFiles))

	return energies, boundStates


def RemoveBoundStateProjection(psi, boundStates):
	assert(pyprop.IsSingleProc())

	projectionList = []

	for curPsiList, L in boundStates:
		if len(curPsiList) > 0:
			#Get the local indices corresponding to the local L
			LFilter = lambda idx: idx.L == L
			indexL = GetLocalCoupledSphericalHarmonicIndices(psi, LFilter)
		
			#Copy the part of psi corresponding to the current L to a 
			#single-L wavefunction to do projection.
			projPsi = curPsiList[0].Copy()
			projPsi.GetData()[:] = psi.GetData()[indexL, :, :]
		
			curProjList = []
			for eigPsi in curPsiList:
				#calculate projection
				proj = projPsi.InnerProduct(eigPsi)
				curProjList.append(proj)
				#remove projection
				psi.GetData()[indexL, :, :] -= proj * eigPsi.GetData()
		
			projectionList.append(curProjList)

	return projectionList


def RunRemoveBoundStateProjection(wavefunctionFile, ionizationThreshold=-2.0):
	#Get Boundstate Files
	conf = pyprop.LoadConfigFromFile(wavefunctionFile)
	boundEnergies, boundStates = GetBoundStates(ionizationThreshold, config=conf)

	#Load Wavefunction
	psi = pyprop.CreateWavefunctionFromFile(wavefunctionFile)
	
	absorbedProbability = real(psi.InnerProduct(psi))
	
	projList = RemoveBoundStateProjection(psi, boundStates)
	print "Probability absorbed            = %s" % (absorbedProbability)
	print "Probability not in projected    = %s" % real(psi.InnerProduct(psi))
	print "Probability in projected states = %s" % (sum([sum(abs(array(p1))**2) for p1 in projList]))

	return psi


#-----------------------------------------------------------------------------
#             tools for testing the eigenstates
#-----------------------------------------------------------------------------

def RunTestBoundStates(**args):
	"""
	Tests the eigenstates corresponding to the problem
	set up by args
	"""
	conf = SetupConfig(**args)
	prop = SetupProblem(config=conf)
	
	boundstateFiles = GetBoundStateFiles(config=conf)
	if len(boundstateFiles) != 1:
		raise Exception("GetBoundStateFiles returned more than one file (%s) for arguments '%s'" % (boundstateFiles, args))

	filename = boundstateFiles[0]
	fileConf = GetEigenstateFileConfig(filename)
	eigenvalues = GetEigenstateFileInfo(filename, INFO_Eigenvalues)

	if not CheckCompatibleRadialGrid(conf, fileConf):
		raise Exception("Incompatible radial grids in eigenstatefile %s and args %s" % (filename, args))

	if not CheckCompatibleAngularGrid(conf, fileConf):
		raise Exception("Incompatible angular grids in eigenstatefile %s and args %s" % (filename, args))

	for i, ev in enumerate(eigenvalues):
		pyprop.serialization.LoadWavefunctionHDF(filename, GetEigenvectorDatasetPath(i), prop.psi)

		error = TestEigenstate(prop, prop.psi)
		PrintOut("Eigenvalue %s, Eigenvector Error = %.15f" % (ev, error))
	
		

def TestEigenstate(prop, psi):
	"""
	Tests to which extent psi is an eigenstate of the 
	hamiltonian set up by prop. The error calculated by

	error = || (E - H) |psi> ||

	is returned, where E is the eigenvalue approximated 
	by the rayleigh coefficienct E = <psi | H | psi>
	
	On return, |psi> contains the error vector
	|psi>  <-  (E - H) |psi> 

	"""
	
	psi.Normalize()

	#Calculate H|psi>
	tempPsi = prop.GetTempPsi()
	tempPsi.Clear()
	prop.MultiplyHamiltonian(psi, tempPsi)

	#calculate eigenvalue estimate
	E = tempPsi.InnerProduct(psi)

	#calculate error vector
	psi.GetData()[:] *= E
	psi.GetData()[:] -= tempPsi.GetData()

	error = psi.GetNorm()

	return error

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
