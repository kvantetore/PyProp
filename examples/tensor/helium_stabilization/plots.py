from scipy import interpolate

def GetGradientFile():
	#return "gradient_chessboard.txt"
	return "gradient_uib.txt"

def TraverseGroups(f):
	for i, p in enumerate(f.root.scanParameter[:]):
		grp = f.getNode("/parameter_%i" % i)
		yield p, grp


def MakePlotsStabilizationFreq5ScanDoubleIonizationDPDE():
	for i in range(36):
		E0, E, dpde = FromFileCalculateDoubleIonizationDPDE("output/stabilization_freq5_scan.h5", i)
		clf()
		pcolormesh(E, E, dpde, vmax=0.2)
		title("$dP/dE_1 dE_2$ for $E_0 = %i$" % E0)
		xlabel("$E_1$ (a.u.)")
		ylabel("$E_2$ (a.u.)")
		savefig("figs/stabilization_freq5_scan/doubleionization_dpde_%i.png" % E0)


def MakePlotsStabilizationFreq5ScanDoubleIonizationDPDELog():
	for i in range(10,11):
		E0, E, dpde = FromFileCalculateDoubleIonizationDPDE("output/stabilization_freq5_scan.h5", i)
		clf()
		pcolormesh(E, E, log(dpde), vmin=-7)
		title("$dP/dE_1 dE_2$ for $E_0 = %i$" % E0)
		xlabel("$E_1$ (a.u.)")
		ylabel("$E_2$ (a.u.)")

def SaveIonizationDataFreq5Scan():
	f = tables.openFile("output/stabilization_freq5_scan.h5")
	try:
		E0, total = zip(*map(lambda p: (p[0], p[1]._v_attrs.TotalIonization), TraverseGroups(f)))
		E0, single = zip(*map(lambda p: (p[0], p[1]._v_attrs.SingleIonization), TraverseGroups(f)))
		E0, absorbed = zip(*map(lambda p: (p[0], p[1]._v_attrs.AbsorbedProbability), TraverseGroups(f)))
		E0, sym = zip(*map(lambda p: (p[0], p[1]._v_attrs.SymmetrizedProbability), TraverseGroups(f)))
		E0, anti = zip(*map(lambda p: (p[0], p[1]._v_attrs.AntiSymmetrizedProbability), TraverseGroups(f)))
	finally:
		f.close()

	f = tables.openFile("output/stabilization_freq5_scan_summary.h5", "w")
	f.createArray(f.root, "amplitude", E0)
	f.createArray(f.root, "total", total)
	f.createArray(f.root, "single", single)
	f.createArray(f.root, "absorbed", absorbed)
	f.createArray(f.root, "symetric", sym)
	f.createArray(f.root, "antisymmetric", anti)


def MakePlotDpDomegaPolar(filename, thCut, energy1, energy2=None, vmin=None, vmax=None):
	if energy2 == None:
		energy2 = energy1

	#setup bounds
	rectBound = (.20, .20, .65, .65)
	singleWidth = .025
	singleSpacing = 0.005
	
	#create figure
	fig = figure()
	rectMain = (rectBound[0]+singleWidth, rectBound[1]+singleWidth, rectBound[2]-singleWidth, rectBound[3]-singleWidth)
	axMain = fig.add_axes(rectMain, polar=True, axisbelow=True)
	
	#get data
	en, th, dp = GetDpDomegaDoubleFromFile(filename)
	energy_double, dpde_double = GetDpDeDoubleFromFile(filename)

	#slice at energies and theta1
	idx1 = argmin(abs(en-energy1))
	idx2 = argmin(abs(en-energy2))
	idx3 = argmin(abs(th - thCut))
	dpSlice = dp[idx3, : ,idx1, idx2]

	#interpolate
	th2 = linspace(0,pi,512)
	dp2 = scipy.interpolate.UnivariateSpline(th, dpSlice, s=0)(th2)

	#normalize
	dp2Norm = sum(dp2 * sin(th2)) * diff(th2)[0]
	dp2 /= dp2Norm

	maxValDpDomega = max(dp2)

	#plot
	axMain.plot(th2, dp2, color = UiB_Blue, linestyle = '-')
	axMain.plot(-th2, dp2, color = UiB_Blue, linestyle = '-')


	#fig.figurePatch.set_alpha(1.0)
	PaperUpdatePolarFigure(fig)
	axMain.set_position(GetOptimalAxesPosition(fig, axMain))
	axMain.set_xlim([0, maxValDpDomega])

	#Update r and theta grid lines/labels
	rGridLines = linspace(0, maxValDpDomega, 5)[1:-1]
	rgrids(rGridLines, ["" for p in rGridLines])
	thetagrids(range(0,360,45), [r"%s$^\circ$" % t for t in range(0,360,45)])

	#put arrow indicating direction of other electron
	axMain.plot([0,0], [0, maxValDpDomega/3.0], "k-", linewidth=0.8)
	axMain.plot([0], [maxValDpDomega/3.0], "k>", markersize=3.0)

	draw()

	return fig


def MakePlotStabilizationFreq5ScanIonization():
	f = tables.openFile("output/stabilization_freq5_scan.h5")
	try:
		E0, total = zip(*map(lambda p: (p[0], p[1]._v_attrs.TotalIonization), TraverseGroups(f)))
		E0, single = zip(*map(lambda p: (p[0], p[1]._v_attrs.SingleIonization), TraverseGroups(f)))
		#E0, absorbed = zip(*map(lambda p: (p[0], p[1]._v_attrs.AbsorbedProbability), TraverseGroups(f)))
		#E0, sym = zip(*map(lambda p: (p[0], p[1]._v_attrs.SymmetrizedProbability), TraverseGroups(f)))
		#E0, anti = zip(*map(lambda p: (p[0], p[1]._v_attrs.AntiSymmetrizedProbability), TraverseGroups(f)))
	finally:
		f.close()

	plot(E0, total, "b-", label="Total Ionization")
	plot(E0, single, "g--", label="Single Ionization")
	plot(E0, array(total) - array(single), "r:", label="Double Ionization")
	xlabel("$E_0$ (a.u.)")
	ylabel("Probability")
	legend(loc="lower right")

	savefig("figs/stabilization_freq5_scan/ionization_probability_scan.png")
	

def MakePlotStabilizationFreq5ScanSingleIonizationDPDE():
	E0, E, dpde = FromFileCalculateSingleIonizationDPDE("output/stabilization_freq5_scan.h5")
	dpde = array(dpde)
	pcolormesh(E, E0, dpde)

	savefig("figs/stabilization_freq5_scan/singleionization_dpde_scan.png")


#convergence
def GetDpDomegaDoubleFromFile(filename, phiEvalType='avg'):
	f = tables.openFile(filename)
	try:
		if phiEvalType == 'avg':
			e = f.root.dpdomega_double_avg._v_attrs.energy
			th = f.root.dpdomega_double_avg._v_attrs.theta
			dp = f.root.dpdomega_double_avg[:]
		elif phiEvalType == 'coplanar':
			e = f.root.dpdomega_double_coplanar._v_attrs.energy
			th = f.root.dpdomega_double_coplanar._v_attrs.theta
			dp = f.root.dpdomega_double_coplanar[:]
	finally:
		f.close()
	return e, th, dp

def GetDpDomegaSingleFromFile(filename):
	f = tables.openFile(filename)
	try:
		e = f.root.dpdomega_single._v_attrs.energy
		th = f.root.dpdomega_single._v_attrs.theta
		dp = f.root.dpdomega_single[:]
	finally:
		f.close()
	return e, th, dp

def GetDpDeDoubleFromFile(filename):
	f = tables.openFile(filename)
	try:
		e = f.root.dpde_double._v_attrs.energy
		dp = f.root.dpde_double[:]
	finally:
		f.close()
	return e, dp

def GetDpDeSingleFromFile(filename):
	f = tables.openFile(filename)
	try:
		e = f.root.dpde_single._v_attrs.energy
		dp = f.root.dpde_single[:]
	finally:
		f.close()
	return e, dp


def InterpolateDpDe2D(oldEnergies, dpde, newEnergies):
	interpolator = interpolate.RectBivariateSpline(oldEnergies, oldEnergies, dpde, kx=3, ky=3)
	dpdeNew = interpolator(newEnergies, newEnergies)
	return dpdeNew

def MakePlotConvergence():
	folder = "convergence/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0"
	e0 = 1
	dt = "1e-02"

	#reference
	filename = "convergence/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0/stabilization_I_%s_kb20_dt_1e-02_T_15.1_0022.h5" % (e0)
	energy, th, dp_ref = GetDpDomegaFromFile(filename)
	dpdomega_ref = sum(sum(dp_ref, axis=3), axis=2) * diff(energy)[0]**2
	energy2, dpde_double_ref = GetDpDeDoubleFromFile(filename)
	energy3, dpde_single_ref = GetDpDeSingleFromFile(filename)
	dp_ref_norm = sum(sum(sum(dp_ref, axis=3), axis=2) * outer(sin(th),sin(th))) * diff(energy)[0]**2 * diff(th)[0]**2 * (2*pi)**2
	dpde_double_ref_norm = sum(dpde_double_ref) * diff(energy2)[0]**2 
	dpde_single_ref_norm = sum(dpde_single_ref) * diff(energy3)[0]

	print "folder = %s" % folder

	f1 = figure()
	ax1 = f1.get_axes()

	energyIndex = argmin(abs(energy-(2*5-2.9)/2))
	energyIndex2 = argmin(abs(energy2-(2*5-2.9)/2))

	#for i, fileIndex in enumerate([15, 15+7, 29]):
	for i, fileIndex in enumerate([15, 22, 25, 29]):
		filename = os.path.join(folder,"stabilization_I_%s_kb20_dt_%s_T_15.1_%04i.h5" % (e0, dt, fileIndex) )

		#dpdomega
		energy, th, dp = GetDpDomegaFromFile(filename)
		dpdomega = sum(sum(dp, axis=3), axis=2) * diff(energy)[0]**2
		
		dp_err = sum(sum(sum(abs(dp - dp_ref), axis=3), axis=2) * outer(sin(th),sin(th))) * diff(energy)[0]**2 * diff(th)[0]**2 * (2*pi)**2
		dpdomega_err = sum(abs(dpdomega - dpdomega_ref) * outer(sin(th), sin(th))) * diff(th)[0]**2 * (2*pi)**2 / dp_ref_norm
	
		#dpde(double)
		energy2, dpde_double = GetDpDeDoubleFromFile(filename)
		dpde_double_err = sum(abs(dpde_double - dpde_double_ref)) * diff(energy2)[0]**2 / dpde_double_ref_norm

		#dpde(single)
		energy3, dpde_single = GetDpDeSingleFromFile(filename)
		dpde_single_err = sum(abs(dpde_single - dpde_single_ref)) * diff(energy3)[0] / dpde_single_ref_norm

		#subplot(3,1,i+1)
		#pcolormesh(th, th, dpdomega)
		#pcolormesh(energy2, energy2, dpde_double)
		#pcolormesh(th, th, dp[:,:,energyIndex,energyIndex])
		#plot(energy, dp[0,0,energyIndex,:])
		#axis("equal")
		plot(th, dp[0,:,energyIndex,energyIndex])
		#plot(energy2, dpde_double[:, energyIndex2], "r")
		#plot(energy2, dpde_double[energyIndex2,:], "g")

		print "file index = %s" % fileIndex
		print "  dp_err = %s" % (dp_err,)
		print "  dpdomega_err = %s" % (dpdomega_err,)
		print "  dpde_double_err = %s" % (dpde_double_err,)
		print "  dpde_single_err = %s" % (dpde_single_err,)

	#referenceplot
	#subplot(2,1,i+2)
	#pcolormesh(th, th, dpdomega_ref)
	#pcolormesh(th, th, dp_ref[:,:,energyIndex, energyIndex])
	#pcolormesh(energy2, energy2, dpde_double_ref)
	#axis("equal")
	#plot(energy2, dpde_double_ref[:, energyIndex2], "k--")
	#plot(energy2, dpde_double_ref[energyIndex2,:], "y--")

	#figure()
	#pcolormesh(th, th, dpdomega_ref - dpdomega)



def MakeDpDePlotScan():
	e0list = [1,5,10,15,30]
	for e0 in e0list:
		f = tables.openFile("output/stabilization_freq5/stabilization_I_%i_kb20_dt_1e-02.h5" % e0)
		try:
			energy_double = f.root.dpde_double._v_attrs.energy
			dpde_double = f.root.dpde_double[:]
			energy_single = f.root.dpde_single._v_attrs.energy
			dpde_single = f.root.dpde_single[:]
			MakeDpDePlot(energy_double, dpde_double, energy_single, dpde_single)

			ax = gcf().axes[0]
			ax.set_title("$E_0 = %i$, $\omega = 5$" % e0)
			draw()

		finally:
			f.close()
	

def CalculatePartialIonizationProbability(energy_double, dpde_double, minE, maxE, doPlot=False):
	"""
	Calculates the double ionization probability for getting electrons with
	combined energy between minE and maxE (minE < E1 + E2 < maxE)

	dpde_double is the calculated differential probability at equidistant energies
	energy_double is the list of the corresponding energies in dpde_double
	"""

	E1, E2 = meshgrid(energy_double, energy_double)
	#flatten to make indices from find work
	dp = dpde_double.copy().reshape(dpde_double.size)
	#zero out energies outside the requested range
	dp[find(E1 + E2 < minE)] = 0
	dp[find(E1 + E2 >= maxE)] = 0
	#reshape back to original shape
	dp = dp.reshape(dpde_double.shape)

	#plot if requested
	if doPlot:
		pcolormesh(energy_double, energy_double, dp)
		title("minE = %f, maxE = %f" % (minE, maxE))

	#integrate to get partial ionization probability
	#partialProb = sum(triu(dp)) * diff(energy_double)[0]**2
	partialProb = sum(dp) * diff(energy_double)[0]**2

	return partialProb

def CalculatePartialIonizationProbabilityScan(fileList, limits):
	partialList = []
	e0List = []
	for filename in fileList:
		f = tables.openFile(filename)
		try:
			energy_double = f.root.dpde_double_L_1._v_attrs.energy
			dpde_double = f.root.dpde_double_L_1[:]
			
			a0 = f.root.wavefunction._v_attrs.configObject.get("PulseParameters", "amplitude")
			frequency = f.root.wavefunction._v_attrs.configObject.get("PulseParameters", "frequency")
			e0 = float(a0) * float(frequency)
			e0List.append(e0)

			p1 = CalculatePartialIonizationProbability(energy_double, dpde_double, limits[0], limits[1])
			p2 = CalculatePartialIonizationProbability(energy_double, dpde_double, limits[1], limits[2])
			p3 = CalculatePartialIonizationProbability(energy_double, dpde_double, limits[2], limits[3])
			p4 = CalculatePartialIonizationProbability(energy_double, dpde_double, limits[3], limits[4])
			partialList.append([p1, p2, p3, p4])
		finally:
			f.close()
	
	return array(e0List), array(partialList)

def GetIonizationProbabilityScan(fileList, symmetry="sym"):
	ionList = []
	for fileName in fileList:
		f = tables.openFile(fileName)
		try:
			if symmetry == "sym":
				doubleIon = f.root.wavefunction._v_attrs.DoubleIonization
				singleIon = f.root.wavefunction._v_attrs.SingleIonization
			elif symmetry == "antisym":
				doubleIon = f.root.wavefunction._v_attrs.AntiDoubleIonization
				singleIon = f.root.wavefunction._v_attrs.AntiSingleIonization
			else:
				raise Exception("Unknown symmetry (%s), should be either 'sym' or 'antisym'!" % symmetry)
				
			a0 = f.root.wavefunction._v_attrs.configObject.get("PulseParameters", "amplitude")
			frequency = f.root.wavefunction._v_attrs.configObject.get("PulseParameters", "frequency")
			e0 = float(a0) * float(frequency)
			
			ionList.append([e0, singleIon, doubleIon])
		finally:
			f.close()

	e0list, singleList, doubleList = zip(*ionList)
	return e0list, array(singleList), array(doubleList)



def CalculatePartialAngularIonizationProbability(energy_double, theta, dpdomegade_double, minE, maxE):
	"""
	Calculates the angular double ionization probability for getting electrons with
	combined energy between minE and maxE (minE < E1 + E2 < maxE)
	"""

	E1, E2 = meshgrid(energy_double, energy_double)
	numE = len(energy_double)
	numTheta = len(theta)
	dpSize = len(energy_double)**2
	dpdomega_double = numpy.zeros((numTheta, numTheta))

	for i in range(numTheta):
		for j in range(numTheta):
			#flatten to make indices from find work
			dp = dpdomegade_double[i,j,:,:].copy().reshape(dpSize)
			
			#zero out energies outside the requested range (and E1 < E2)
			dp[find(E1 + E2 < minE)] = 0
			dp[find(E1 + E2 >= maxE)] = 0
			dp[find(E1 < E2 )] = 0
			
			dpdomega_double[i,j] = sum(dp)

	#Scale by energy spacing
	dpdomega_double[:] *= diff(energy_double)[0]**2

	return dpdomega_double


def CalculateRadialPartialAngularIonizationProbability(energy_double, theta, dpdomegade_double, energyPos1, energyPos2, energyRadius):
	"""
	Calculates the angular double ionization probability for getting electrons with
	combined energy between minE and maxE inside the circle 
	"""

	E1, E2 = meshgrid(energy_double-energyPos1, energy_double-energyPos2)
	numE = len(energy_double)
	numTheta = len(theta)
	dpSize = len(energy_double)**2
	dpdomega_double = numpy.zeros((numTheta, numTheta))

	for i in range(numTheta):
		for j in range(numTheta):
			#flatten to make indices from find work
			dp = dpdomegade_double[i,j,:,:].copy().reshape(dpSize)
			
			#zero out energies outside the requested energy circle
			dp[find(sqrt(E1**2 + E2**2) > energyRadius)] = 0
			
			dpdomega_double[i,j] = sum(dp)

	#Scale by energy spacing
	dpdomega_double[:] *= diff(energy_double)[0]**2

	return dpdomega_double
