

def GetBoundStates(ionizationThreshold=-2.0, **args):
	filename = GetEigenvectorFilename(**args)

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

	return curEnergyList, curPsiList


def RemoveBoundStateProjection(psi, boundStates):
	"""
	Remove bound state projection projection from psi.

	boundStates should be a list of (bound) eigenstates 
	wavefunctions,	and psi should be a Wavefunction.
	"""
	for b in boundStates:
		curProj = psi.InnerProduct(b)
		psi.GetData()[:] -= curProj * b.GetData()[:]



def CalculateBoundstateProjections(psi, boundStates):
	"""
	Calculate bound state projection projection from psi.

	boundStates should be a list of (bound) eigenstates 
	wavefunctions,	and psi should be a Wavefunction.
	"""
	populations = []
	for b in boundStates:
		populations += [abs(psi.InnerProduct(b))**2]


	return populations


def CreateSuperpositionOfBoundstates(**args):
	conf = SetupConfig(**args)
	psi = pyprop.CreateWavefunction(conf)
	psi.Clear()
	E, V = GetBoundStates(config = conf)

	for b in V:
		psi.GetData()[:] += b.GetData()[:]

	#psi.Normalize()
	
	return psi

