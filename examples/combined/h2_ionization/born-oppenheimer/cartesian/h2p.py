

def GetH2pGroundstate(**args):
	assert("configFile" not in args)
	args["configFile"] = "config-h2p.ini"
	args["eigenvalueCount"] = 1
	
	prop = SetupProblem(**args)
	solver = pyprop.ArpackSolver(prop)
	solver.Solve()

	E = solver.GetEigenvalues()[0]
	solver.SetEigenvector(prop.psi, 0)

	return E, prop.psi


