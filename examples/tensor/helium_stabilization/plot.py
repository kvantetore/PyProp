from pylab import *
from numpy import *
import tables
from scipy import interpolate

execfile("load_cmap.py")


def LaTeXFigureSettings(fig_width=16, fig_height=12):
	cm_to_inch = 0.393700787
	fig_width *= cm_to_inch
	fig_height *= cm_to_inch
	fig_size =  [fig_width,fig_height]
	params = {'backend': 'ps',
			  'axes.labelsize': 16,
			  'text.fontsize': 16,
			  'legend.fontsize': 14,
			  'xtick.labelsize': 14,
			  'ytick.labelsize': 14,
			  'text.usetex': False,
			  'figure.figsize': fig_size,
			  'font.family' : 'serif',
			  'font.serif' : 'Times New Roman'}
	pylab.rcParams.update(params)


def MakePlotSingleIonizationEnergyDistribution(filename="dpde_scan.h5"):
	f = tables.openFile(filename, "r")
	try :
		E = array(f.root.energy[:][0])
		E0 = array(f.root.intensity[:])
		dpde = array(f.root.dpde[:])
	finally :
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


def MakePlotDoubleIonizationEnergyDistribution(filename):
	gradientFile = "gradient_uib.txt"

	f = tables.openFile(filename, "r")
	try :
		E = array(f.root.dpde_double._v_attrs.energy)
		dpde = array(f.root.dpde_double)
	finally :
		f.close()

	#2d interpolation
	newEnergies = linspace(E[0], E[-1], 500)
	interpolator = interpolate.RectBivariateSpline(E, E, dpde, kx=3, ky=3)
	dpdeNew = interpolator(newEnergies, newEnergies)

	#load color map
	cmap = LoadColormap(gradientFile, reverse=True)
	figure()
	#pcolormesh(newEnergies, newEnergies, dpdeNew, cmap=cmap, vmin=1e-9, vmax=1e-1)
	pcolormesh(newEnergies, newEnergies, dpdeNew, cmap=cmap)
	xlabel("Energy (a.u.)")
	ylabel("Energy (a.u.)")



def MakePlotSingleAndDoubleIonizationEnergyDistribution(filename):
	figRatio = 16.0 / 12.0
	offset = figRatio
	gradientFile = "gradient_uib.txt"

	f = tables.openFile(filename, "r")
	try :
		energyDouble = array(f.root.dpde_double._v_attrs.energy)
		dpdeDouble = array(f.root.dpde_double)

		energySingle = array(f.root.dpde_single._v_attrs.energy)
		dpdeSingle = array([f.root.dpde_single[:]])
		
	finally :
		f.close()

	#2d interpolation
	#newEnergies = linspace(energyDouble[0], energyDouble[-1], 500)
	newEnergies = energySingle[:]
	interpolator = interpolate.RectBivariateSpline(energyDouble, energyDouble, dpdeDouble, kx=3, ky=3)
	dpdeDoubleNew = interpolator(newEnergies, newEnergies)

	#load color map
	cmap = LoadColormap(gradientFile, reverse=True)
	figure()

	#plot single dpde
	ax1 = axes([0.1, 0.12 , 0.02, 0.8])
	ax1.set_xticks(())
	ax1.pcolormesh(r_[0:2], energySingle, transpose(dpdeSingle), cmap=cmap)
	ylabel("Energy (a.u.)")
	
	ax2 = axes([0.12,0.1,0.8,0.02])
	ax2.set_yticks(())
	ax2.pcolormesh(energySingle, r_[0:2], dpdeSingle, cmap=cmap)
	xlabel("Energy (a.u.)")

	ax3 = axes([0.13,0.13,0.79,0.79])
	ax3.pcolormesh(newEnergies, newEnergies, dpdeDoubleNew, cmap=cmap)
	ax3.set_xticks(())
	ax3.set_yticks(())
	#xlabel("Energy (a.u.)")
	#ylabel("Energy (a.u.)")

	show()

