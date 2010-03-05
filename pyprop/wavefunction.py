import core

def CreateWavefunctionInstance(representation, allocateData=True):
	#Create instance
	print "    Creating instance"
	rank = len(representation.GetFullShape())
	psi = CreateInstanceRank("core.Wavefunction", rank)
	
	#Set reresentation
	print "    Setting representation"
	psi.SetRepresentation(representation)
	
	#Allocate data
	if allocateData:
		print "    Allocating data"
		psi.AllocateData()
	
	return psi
	
	
