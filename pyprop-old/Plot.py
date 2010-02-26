#1D
def Plot1D(prop, **args):
	data = prop.psi.GetData()
	grid = prop.GetGrid()[0]
	if IsMaster():
		plot(grid, real(conj(data) * data), **args)
	
#2D
def Plot2DFull(prop, **args):
	if 'shading' not in args:
		args['shading'] = 'flat'
	data = GetFullWavefunctionData(prop)

	xvector = prop.psi.GetRepresentation().GetLocalGrid(0)
	yvector = prop.psi.GetRepresentation().GetLocalGrid(1)
	x, y = meshgrid(xvector,yvector)
	
	plotdata = real(conj(data) * data)
	plotdata = plotdata * 100 / max(plotdata.flatten())
	if IsMaster():
		pcolor(x, y, plotdata,  **args)
	

def Plot2DRank(prop, rank, **args):
	grid = prop.GetGrid()[rank]
	data = GetFullWavefunctionData(prop)
	plotdata = sum(abs(data)**2, 1 - rank) * diff(grid)[0] #assume equidistant grid
	if IsMaster():
		plot(grid, plotdata, **args)

def GetProjectedWavefunction(prop, rankList):
	"Integrates the absolute square of the wavefunction in all other "
	"ranks than those specified in ranklist, and returns the resulting "
	"len(rankList)-dimensional array"

	#data to integrate
	data = abs(prop.psi.GetData())**2

	#sort rankList in reverse order
	rankList.sort()
	rankList.reverse()

	#create a list of ranks to be summed over
	removeRanks  = range(prop.psi.GetRank())
	for rank in rankList:
		del removeRanks[rank]
	#in reversed order
	removeRanks.reverse()

	#sum over ranks
	for rank in removeRanks:
		data = sum(data, rank)

	return data
	
#4D
def Plot4DRank2D(prop, rank1, rank2, **args):
	if 'shading' not in args:
		args['shading'] = 'flat'

	if rank1 < rank2:
		rank1, rank2 = rank2, rank1
	
	data = abs(prop.psi.GetData())
	
	removeRanks = range(prop.psi.GetRank())
	del removeRanks[rank1]
	del removeRanks[rank2]
	removeRanks.reverse()
	for rank in removeRanks:
		data = sum(data, rank)
	
	data = 100 * data / max(data.flatten())
	pcolor(data, **args)

def Plot4DRank1D(prop, rank, **args):
	data = abs(prop.psi.GetData())
	
	removeRanks = range(prop.psi.GetRank())
	del removeRanks[rank]
	removeRanks.reverse()
	for rank in removeRanks:
		data = sum(data, rank)
	
	grid = prop.GetGrid()[rank]
	plot(grid, data, **args)

