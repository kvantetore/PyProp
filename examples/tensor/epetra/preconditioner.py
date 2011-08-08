from pyprop import EpetraPotential

class TwoElectronPreconditionerIfPack:
    """
    Preconditioner for GMRES for solving systems of the type

    (1)     (a H + b S) x = y  ->  x

    where H is the Hamiltonian, S is the overlap matrix
    and a and b are complex scaling factors (scalingH and scalingS).

    The preconditioner sets up an Epetra potential of all potentials
    in the 'potential_evaluation' list given in the config section
    for the preconditioner. This is the passed to IfPack for computation
    of approximate L and U factors, where the overlap between procs are
    given by the 'overlap_level' config keyword.
    """

    def __init__(self, psi):
        self.Rank = psi.GetRank()
        self.Psi = psi

    def ApplyConfigSection(self, conf):
        self.OverlapSection = conf.Config.GetSection("OverlapPotential")
        self.PotentialSections = [conf.Config.GetSection(s) for s in conf.potential_evaluation]
        self.Cutoff = conf.cutoff
        self.OverlapLevel = conf.overlap_level

    def SetHamiltonianScaling(self, scalingH):
        self.HamiltonianScaling = scalingH

    def SetOverlapScaling(self, scalingS):
        self.OverlapScaling = scalingS

    def GetHamiltonianScaling(self):
        return self.HamiltonianScaling

    def GetOverlapScaling(self):
        return self.OverlapScaling

    def Setup(self, prop):
        """
        Set up a tensor potential for overlap potential and all other potentials
        and add them together, assuming they have the same layout
        ending up with a potential containing S + scalingH * (P1 + P2 + ...)

        """
        pot = self.GeneratePotentials(prop)
        pyprop.PrintOut("Setting up IfPack ILU preconditioner for potentials: %s" % pot.Name)
        self.Solver = pyprop.CreateInstanceRank("core.IfpackRadialPreconditioner", 3)
        self.Solver.Setup(prop.psi.GetData(), pot.Potential, self.OverlapLevel)

    def GeneratePotentials(self, prop):
        """
        Compute potentials and add them together, forming an Epetra potential. 
        """
        potentialList = []

        for configSection in self.PotentialSections:
            #generate potential 
            tensorPot = prop.BasePropagator.GeneratePotential(configSection)
            localBasisPairs = [geom.GetBasisPairs() for geom in tensorPot.GeometryList]
            if tensorPot.Name == "OverlapPotential":
                tensorPot.PotentialData *= self.GetOverlapScaling()
            else:
                tensorPot.PotentialData *= self.GetHamiltonianScaling()

            #Setup a new Epetra potential
            epetraPot = EpetraPotential(self.Psi)
            epetraPot.ApplyConfigSection(configSection)

            #check if current potential can be consolidated with an existing one
            for existingPot in potentialList:
                if existingPot.CanConsolidate(epetraPot):
                    existingPot.AddTensorPotentialData(tensorPot.PotentialData, localBasisPairs, self.Cutoff)
                    existingPot.Name += "+" + tensorPot.Name
                    tensorPot = None
                    epetraPot = None
                    break
                    
            #add new potential to list
            if epetraPot != None:
                epetraPot.AddTensorPotentialData(tensorPot.PotentialData, localBasisPairs, self.Cutoff)
                potentialList.append(epetraPot)
                tensorPot = None

        #Global assemble all potentials
        for pot in potentialList:
            pot.GlobalAssemble()

        return potentialList[0]

    def Solve(self, psi):
        self.Solver.Solve(psi.GetData())
