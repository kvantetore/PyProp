from pylab import *
from numpy import *
import tables

def MakePlotSingleIonizationEnergyDistribution(filename="dpde_scan.h5"):
	f = tables.openFile(filename, "r")
	try:
		E = array(f.root.energy[:][0])
		E0 = array(f.root.intensity[:])
		dpde = array(f.root.dpde[:])
	finally:
		f.close()

	#filter zero-bins:
	dpdeList = [dpde[:,i] for i in range(dpde.shape[1])]
	fE, fdpde = zip(*filter(lambda d: sum(d[1])>0.0, zip(E, dpdeList)))
	fE = array(fE)
	fdpde = transpose(array(fdpde))

	pcolormesh(fE, E0, fdpde)
	xlim(0, 14)
	ylim(1, 35)
	title("dP/dE for $\omega=5$")
	xlabel("Electron Energy (a.u.)")
	ylabel("Field Strength ($E_0$)")

