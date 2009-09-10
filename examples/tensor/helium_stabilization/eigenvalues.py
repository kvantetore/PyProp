
#------------------------------------------------------------------------------------
#       Eigenvalue Functions using inverse iterations + pIRAM
#------------------------------------------------------------------------------------

def LoadBoundstate(psi, conf, L, eigenvalueIndex):
	"""
	Loads the boundstate of a given L-shell, corresponding
	to a given configuration.

	eigenvalueIndex is the index into the sorted eigenvalue list,
	such that eigenvalueIndex==0 corresponds to the lowest eigenstate
	for that L-shell.
	"""

	#Get boundstates filename	
	filename = GetBoundstatesFilename2(config=conf, L=L)

	#Get sorted index in the eigenvalue list
	f = tables.openFile(filename, "r")
	try:
		eigenvalues = f.root.Eig.Eigenvalues[:]
		idx = argsort(eigenvalues)[eigenvalueIndex]
	finally:
		f.close()

	#Load state
	return LoadBoundstateIndex(psi, filename, idx)
	

def LoadBoundstateIndex(psi, filename, eigenvectorIndex):
	"""
	Loads an eigenstate from a file <filename> into a wavefunction <psi>. 
	This function works in parallel.

	Because eigenstates can be found independently for each "L", 
	the number and order of angular indices are different in the
	stored eigenstates, from that of the wavefunction to load the 
	eigenstate into.

	The eigenstate config file is used to create an array of CoupledIndex elements
	to which the local CoupledIndex elements from the wavefunction are compared

	"""
	angularRank = 0

	PrintOut("Loading Bound State #%i from %s" % (eigenvectorIndex, filename))
	
	f = tables.openFile(filename, "r")
	try:
		#Get eigenvector dataset
		eigVec = f.getNode(GetEigenvectorDatasetPath(eigenvectorIndex))

		#Get corresponding eigenvalue
		E = f.root.Eig.Eigenvalues[eigenvectorIndex]

		#find angular indices of the saved eigenstate
		conf = pyprop.LoadConfigFromFile(filename, "/Eig")
		eigCoupledIndices = array([i for i in conf.AngularRepresentation.index_iterator])

		#find local angular indices of the wavefunction to load into
		repr = psi.GetRepresentation()
		angRepr = repr.GetRepresentation(angularRank)
		localIndices = array(repr.GetLocalGrid(angularRank), dtype=int)

		#Load the local angular indices into the wavefunction
		psi.Clear()
		for curLocalIdx, curGlobalIdx in enumerate(localIndices):
			#Get the CoupledIndex to this angular index
			curCoupledIndex =  angRepr.Range.GetCoupledIndex(curGlobalIdx)

			#Find out the index in the eigenvector corresponding
			eigIdx = find(eigCoupledIndices == curCoupledIndex)
			if len(eigIdx) == 0:
				continue
			elif len(eigIdx) > 1:
				raise Exception("What? %s" % eigIdx)
			else:
				eigIdx = eigIdx[0]
				

			#Slice the wavefunction at the correct angular index
			psiSlice = [s_[:]]*psi.GetRank()
			psiSlice[angularRank] = curLocalIdx

			#Slice the eigenstate at the correct angular index
			eigSlice = [s_[:]]*psi.GetRank()
			eigSlice[angularRank] = eigIdx

			psi.GetData()[curLocalIdx,:,:] = eigVec[eigIdx,:,:]

	finally:
		f.close()

	return E


def GetBoundstatesFilename(**args):
	radialPostfix = "_".join(GetRadialGridPostfix(**args))
	angularPostfix = "_".join(GetAngularGridPostfix(**args))
	locations = ["output/boundstates", "tore/output/boundstates"]
	boundstatesFilename = filter(os.path.exists, ["%s/boundstates_%s_%s.h5" % (loc, radialPostfix, angularPostfix) for loc in locations])
	
	return boundstatesFilename[0]


def GetBoundstatesFilename2(**args):
	customPostfix = ""
	if "customPostfix" in args:
		customPostFix = "_%s" % args["customPostfix"]
	radialPostfix = "_".join(GetRadialGridPostfix(**args))
	angularPostfix = "_".join(GetAngularGridPostfix(**args))
	location = "output/boundstates"
	boundstatesFilename = "%s/boundstates_%s_%s%s.h5" % (location, radialPostfix, angularPostfix, customPostfix)
	
	return boundstatesFilename


def RunFindBoundstates(**args):
	#use default eigenvaluecount and eigenvaluebasissize
	#if not specified
	eigenvalueCount = args.get("eigenvalueCount", 30)
	eigenvalueBasisSize = args.get("eigenvalueBasisSize", 70)
	args["eigenvalueCount"] = eigenvalueCount
	args["eigenvalueBasisSize"] = eigenvalueBasisSize
	shift = args.get("shift", -2.9)
	args["shift"] = shift

	outFolder = "output/boundstates"
	if pyprop.ProcId == 0:
		if not os.path.exists(outFolder):
			os.makedirs(outFolder)
	outFileName = GetBoundstatesFilename2(**args)
	FindEigenvaluesInverseIterations(outFileName=outFileName, **args)
			

def FindEigenvaluesInverseIterations(config="config_eigenvalues.ini", \
	outFileName="out/eig_inverseit.h5", **args):
	"""
	Finds a subset of the eigenvalues/eigenvectors using inverse iterations
	and pIRAM and saves the result in outFileName
	"""
	shift = args.get("shift", -2.9)
	args.pop("shift")
	#args["shift"] = shift
	useAnasazi = args["useAnasazi"]

	if useAnasazi:
		solver, shiftInvertSolver, eigenvalues = FindEigenvaluesNearShiftAnasazi(shift, config=config, **args)
	else:
		solver, shiftInvertSolver, eigenvalues = FindEigenvaluesNearShift(shift, config=config, **args)
	SaveEigenvalueSolverShiftInvert(outFileName, solver, shiftInvertSolver, shift)

	#Print Statistics
	if pyprop.ProcId == 0:
		shiftInvertSolver.PrintStatistics()
	PrintOut("SHIFT = %s" % shift)
	PrintOut("GMRESERRORS = %s" % shiftInvertSolver.Solver.GetErrorEstimateList())
	PrintOut("EIGENVALUES = %s" % eigenvalues)

	return solver, shiftInvertSolver


def FindEigenvaluesNearShift(shift, **args):
	"""
	Calculates eigenvalues and eigenvectors for **args around shift
	by using inverse iterations and pIRAM
	"""

	#Setup Problem
	prop = SetupProblem(silent=True, eigenvalueShift=shift, disablePreconditioner=True, **args)

	#Setup shift invert solver in order to perform inverse iterations
	shiftInvertSolver = pyprop.GMRESShiftInvertSolver(prop)
	prop.Config.Arpack.inverse_iterations = True
	prop.Config.Arpack.matrix_vector_func = shiftInvertSolver.InverseIterations

	#Setup eiganvalue solver
	solver = pyprop.PiramSolver(prop)
	solver.Solve()

	#Get the converged eigenvalues
	ev = solver.Solver.GetEigenvalues().copy()
	estimates = solver.Solver.GetConvergenceEstimates().copy()
	idx = where(estimates < 0)[0]
	eigenvalues = ev[idx]

	#convert from shift inverted eigenvalues to "actual" eigenvalues
	eigenvalues = 1.0 / eigenvalues + shift
	PrintOut("eigenvalues = %s" % eigenvalues)

	return solver, shiftInvertSolver, eigenvalues


def GetEigenvectorDatasetName(eigenvectorIndex):
	return "Eigenvector_%i" % eigenvectorIndex 


def GetEigenvectorDatasetPath(eigenvectorIndex):
	return "/Eig/%s" % GetEigenvectorDatasetName(eigenvectorIndex) 


def SaveEigenvalueSolverShiftInvert(filename, solver, shiftInvertSolver, shift):
	"""
	Saves the output of FindEigenvaluesNearShift, including error estimates 
	to a hdf file.
	"""

	#Get eigenvalue error estimates
	errorEstimatesPIRAM = solver.Solver.GetErrorEstimates()
	convergenceEstimatesEig = solver.Solver.GetConvergenceEstimates()
	errorEstimatesGMRES = shiftInvertSolver.Solver.GetErrorEstimateList()


	#Get eigenvalues
	prop = solver.BaseProblem
	E = 1.0 / array(solver.GetEigenvalues()) + shift

	#remove file if it exists
	try:
		if os.path.exists(filename):
			if pyprop.ProcId == 0:
				os.remove(filename)
	except:
		PrintOut("Could not remove %s (%s)" % (filename, sys.exc_info()[1]))

	#Store eigenvalues and eigenvectors
	PrintOut("Now storing eigenvectors...")
	for i in range(len(E)):
		solver.SetEigenvector(prop.psi, i)
		prop.SaveWavefunctionHDF(filename, GetEigenvectorDatasetPath(i))

	if pyprop.ProcId == 0:
		RemoveExistingDataset(filename, "/Eig/Eigenvalues")
		RemoveExistingDataset(filename, "/Eig/ErrorEstimateListGMRES")
		RemoveExistingDataset(filename, "/Eig/ErrorEstimateListPIRAM")
		RemoveExistingDataset(filename, "/Eig/ConvergenceEstimateEig")
		h5file = tables.openFile(filename, "r+")
		try:
			#myGroup = h5file.createGroup("/", "Eig")
			myGroup = h5file.getNode("/Eig")
			h5file.createArray(myGroup, "Eigenvalues", E)
			h5file.createArray(myGroup, "ErrorEstimateListGMRES", errorEstimatesGMRES)
			h5file.createArray(myGroup, "ErrorEstimateListPIRAM", errorEstimatesPIRAM)
			h5file.createArray(myGroup, "ConvergenceEstimateEig", convergenceEstimatesEig)

			#Store config
			myGroup._v_attrs.configObject = prop.Config.cfgObj
			
			#PIRAM stats
			myGroup._v_attrs.opCount = solver.Solver.GetOperatorCount()
			myGroup._v_attrs.restartCount = solver.Solver.GetRestartCount()
			myGroup._v_attrs.orthCount = solver.Solver.GetOrthogonalizationCount()
		finally:
			h5file.close()

#------------------------------------------------------------------------------------
#                       Eigenvalue Functions not using pyprop
#------------------------------------------------------------------------------------

def FindEigenvaluesJD(howMany, shift, tol = 1e-10, maxIter = 200, dataSetPath="/", \
	configFileName="config_eigenvalues.ini", L=0, lmax=4, outFileName = "eig_jd.h5", \
	preconType = None):
	"""
	Find some eigenvalues for a given L-subspace using Jacobi-Davidson method
	"""

	#Set up problem
	idxIt = pyprop.DefaultCoupledIndexIterator(lmax = lmax, L = L)
	prop = SetupProblem(config = configFileName, index_iterator = idxIt)

	#Set up hamilton matrix
	H_ll = SetupPotentialMatrixLL(prop,[0,1])
	H = H_ll.to_sss()

	#Set up overlap matrix
	S_ll = SetupPotentialMatrixLL(prop,[2])
	S = S_ll.to_sss()

	#Set up preconditioner
	Precon = None
	if preconType:
		Precon = preconType(H_ll)

	#Call Jacobi-Davison rountine
	numConv, E, V, numIter, numIterInner = \
		pysparse.jdsym.jdsym(H, S, Precon, howMany, shift, tol, maxIter, pysparse.itsolvers.qmrs)

	#Store eigenvalues and eigenvectors
	h5file = tables.openFile(outFileName, "w")
	try:
		myGroup = h5file.createGroup("/", "Eig")
		h5file.createArray(myGroup, "Eigenvectors", V)
		h5file.createArray(myGroup, "Eigenvalues", E)
		myGroup._v_attrs.NumberOfIterations = numIter
		myGroup._v_attrs.NumberOfInnerIterations = numIterInner
		myGroup._v_attrs.NumberOfConvergedEigs = numConv
		myGroup._v_attrs.configObject = prop.Config.cfgObj
	finally:
		h5file.close()


	return numConv, E, V, numIter, numIterInner


def FindEigenvaluesDirectDiagonalization(L=0, lmax=3, storeResult=False, checkSymmetry=False, \
	outFileName = "eig_direct.h5", config="config_eigenvalues.ini"):
	"""
	Get energies and eigenstates by direct diagonalization of L-subspace matrix
	"""
	#Coupled spherical index iterator based on given lmax and L
	index_iterator = pyprop.DefaultCoupledIndexIterator(lmax=lmax, L=L)
	
	#Set up problem
	prop = SetupProblem(config=config, index_iterator=index_iterator)

	#Set up hamilton and overlap matrices
	HamiltonMatrix = SetupBigMatrixReal(prop, [0,1])
	OverlapMatrix = SetupBigMatrixReal(prop, [2])

	#Calculate generalized eigenvalues and eigenvectors
	print "Calculating generalized eigenvalues and eigenvectors..."
	sys.stdout.flush()
	E, V = scipy.linalg.eig(HamiltonMatrix, b=OverlapMatrix)

	#Sort eigenvalues and eigenvectors
	sortIdx = argsort(E)
	E = E[sortIdx].real
	V = V[:,sortIdx]

	#Check symmetry of eigenstates
	print "Checking symmetry of eigenstates..."
	psiShape = prop.psi.GetData().shape
	lmax = psiShape[0]
	eps = 1e-9
	symmetryList = []
	if checkSymmetry:
		for stateIdx in range(len(E)):
			
			#Check symmetry of l-component with largest norm
			lSpaceNorm = \
				[scipy.linalg.norm(V[:,stateIdx].reshape(psiShape)[i,:]) for i in range(lmax)]
			largestComponentIdx = argmax(lSpaceNorm)
			v_radial = V[:,stateIdx].reshape(psiShape)[largestComponentIdx,:,:]

			#Check if even/odd or not symmetric
			if scipy.linalg.norm(v_radial - transpose(v_radial)) < eps:
				symmetryList.append("Even")
			elif scipy.linalg.norm(v_radial + transpose(v_radial)) < eps:
				symmetryList.append("Odd")
			else:
				print "E = %s, idx = %s" % (E[stateIdx], stateIdx)
				symmetryList.append("NoSym")

	#Store result
	if storeResult:
		h5file = tables.openFile(outFileName, "w")
		try:
			myGroup = h5file.createGroup("/", "Eig")
			h5file.createArray(myGroup, "Eigenvectors", V)
			h5file.createArray(myGroup, "Eigenvalues", E)
			if symmetryList:
				h5file.createArray(myGroup, "SymmetryList", symmetryList)
			myGroup._v_attrs.configObject = prop.Config.cfgObj
		finally:
			h5file.close()

	#...or return result
	else:
		return prop, HamiltonMatrix, OverlapMatrix, E, V





#------------------------------------------------------------------------------------
#       Eigenvalue Functions using Trilinos::Anasazi
#------------------------------------------------------------------------------------

execfile("AnasaziSolver.py")

def FindEigenvaluesNearShiftAnasazi(shift, **args):
	"""
	Calculates eigenvalues and eigenvectors for **args around shift
	by using inverse iterations and Anasazi
	"""
	
	#Setup Problem
	prop = SetupProblem(silent=True, eigenvalueShift=shift, disablePreconditioner=False, **args)

	###Setup shift invert solver in order to perform inverse iterations
	prop.Config.GMRES.shift = shift
	shiftInvertSolver = pyprop.GMRESShiftInvertSolver(prop)
	prop.Config.Arpack.inverse_iterations = True
	prop.Config.Arpack.matrix_vector_func = shiftInvertSolver.InverseIterations

	solver = AnasaziSolver(prop)
	#solver.ApplyMatrix = lambda src, dst, a, b: solver.Preconditioner.Preconditioner.Multiply(src, dst)
	solver.Solve()
	
	PrintOut("GMRES Error Estimates: %s " % shiftInvertSolver.Solver.GetErrorEstimateList())
	shiftInvertSolver.PrintStatistics()

	#Get the converged eigenvalues
	eigenvalues = solver.Solver.GetEigenvalues().copy()

	#convert from shift inverted eigenvalues to "actual" eigenvalues
	eigenvalues = 1.0 / eigenvalues + shift
	PrintOut("eigenvalues = %s" % eigenvalues)

	return solver, shiftInvertSolver, eigenvalues

