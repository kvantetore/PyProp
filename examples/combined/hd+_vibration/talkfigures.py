

def MakeTalkFigure1():
	propIon = SetupProblem(molecule="d2+", species="Ion", radialScaling=2)
	propNeutral = SetupProblem(molecule="d2+", species="Neutral", radialScaling=2)

	r = propIon.psi.GetRepresentation().GetLocalGrid(0)
	dt = propIon.TimeStep

	pot1 = propNeutral.Propagator.PotentialList[0].GetPotential(dt)[:,0]
	pot2 = propIon.Propagator.PotentialList[0].GetPotential(dt)[:,0]
	pot3 = propIon.Propagator.PotentialList[0].GetPotential(dt)[:,1]

	return r, pot1, pot2, pot3, propIon.Propagator.PotentialList[0]
