
def FindGroundstate(**args):
	"""
	Uses imaginary time to find the groundstate
	to the problem specified in **args
	"""
	args["imtime"] = True

	prop = SetupProblem(**args)

	for t in prop.Advance(10):
		E = prop.GetEnergy()
		print "t = %3.2f, E = %2.8f" % (t, E)

	return prop


def FindEigenstates(**args):
	"""
	Uses pIRAM to find the the lowest eigenvectors of the
	problem specified in **args
	"""
	prop = SetupProblem(**args)

	#use custom initial residual if provided
	initialResidual = args.get("initialResidual")
	if initialResidual != None:
		prop.psi.GetData()[:] = initialResidual.GetData()

	#find eigenstates
	solver = pyprop.ArpackSolver(prop)
	solver.Solve()
	return solver


def TestAntisymmetric(psi, tempPsi=None):
	psi.Normalize()
	if tempPsi == None:
		tempPsi = psi.Copy()

	#Test if function is antisymmetric if x[0] -> x[1], and x[1] -> x[0]
	tempPsi.GetData()[:,:] = psi.GetData()[:,:] + psi.GetData()[:,:].transpose()
	normAntisymmetric = tempPsi.GetNorm()

	tempPsi.GetData()[:,:] = psi.GetData()[:,:] - psi.GetData()[:,:].transpose()
	normSymmetric = tempPsi.GetNorm()

	isAntiSymmetric = False
	if normAntisymmetric < 1e-5 and abs(normSymmetric - 2) < 1e-5:
		isAntiSymmetric = True
	elif abs(normAntisymmetric-2) < 1e-5 and normSymmetric < 1e-5:
		isAntiSymmetric = False
	else:
		print "Warning: psi is neither symmetric nor antisymmetric (%f, %f)" % (normAntisymmetric, normSymmetric)
		isAntiSymmetric = normAntisymmetric < normSymmetric

	return isAntiSymmetric
	

def FindElectronicEnergyCurves(**args):
	args["silent"] = True

	nuclearSeparationList = \
		array([0.4000, 0.4500, 0.5000, 0.5500,\
 		0.6500, 0.7000, 0.7500, 0.8000, 0.9000, \
		1.0000, 1.1000, 1.2000, 1.3000, 1.3500, \
		1.3900, 1.4000, 1.4010, 1.4011, 1.4100, \
		1.4500, 1.5000, 1.6000, 1.7000, 1.8000, \
		1.9000, 2.0000, 2.1000, 2.2000, 2.3000, \
		2.4000, 2.5000, 2.6000, 2.7000, 2.8000, \
		2.9000, 3.0000, 3.1000, 3.2000, 3.3000, \
		3.4000, 3.5000, 3.6000, 3.7000, 3.8000, \
		3.9000, 4.0000, 4.1000, 4.2000, 4.3000, \
		4.4000, 4.5000, 4.6000, 4.7000, 4.8000, \
		4.9000, 5.0000, 5.1000, 5.2000, 5.3000, \
		5.4000, 5.5000, 5.6000, 5.7000, 5.8000, \
		5.9000, 6.0000, 6.1000, 6.2000, 6.3000, \
		6.4000, 6.5000, 6.6000, 6.7000, 6.8000, \
		6.9000, 7.0000, 7.2000, 7.4000, 7.6000, \
		7.8000, 8.0000, 8.2500, 8.5000, 9.0000, \
		9.5000, 10.0000], dtype=double)

	nuclearSeparationList = nuclearSeparationList[::5]

	energyList = []

	#Set up a temporary psi
	prop = SetupProblem(**args)
	tempPsi = prop.psi
	initialResidual = tempPsi.Copy()

	#Solve eigenstate for each nuclear separation
	estimatedTime = False
	for i, separation in enumerate(nuclearSeparationList):
		if estimatedTime:
			print "Estimated Time Remaining: %s" % formatTimedelta(estimatedTime)

		currentTime = - time.time()
		solver = FindEigenstates(nuclearSeparation=separation, initialResidual=initialResidual, **args)
		E = solver.GetEigenvalues()
		restartCount = solver.Solver.GetIterationCount()
		operationCount = solver.Solver.GetOperationCount()

		print "Found %i eigenvalues, used %i restart iterations and %i multiply hamiltonian operations" % (len(E), restartCount, operationCount)
		 
		curEnergyList = []
	
		#Find the energies which correspond to antisymmetric wavefunctions
		psi = solver.BaseProblem.psi
		initialResidual.GetData()[:] = 0
		for j, curE  in enumerate(E):
			solver.SetEigenvector(psi, j)
			initialResidual.GetData()[:] += psi.GetData()
			#Only the symmetric (in space) wavefunctions are bonding!
			if not TestAntisymmetric(psi, tempPsi=tempPsi):
				curEnergyList.append(curE)

		energyList.append(curEnergyList)

		currentTime += time.time()
		estimatedTime = currentTime * (len(nuclearSeparationList) - i - 1)
		
	#Convert it from a list of lists to an array
	try:
		energyList = array(energyList)
	except:
		print "Warning: energyList is not a rectangular array"

	return nuclearSeparationList, energyList


def formatTimedelta(seconds):
		h = seconds // 3600
		m = (seconds % 3600) // 60
		s = seconds % 60
		return "%02i:%02i:%02i" % (h,m,s)

