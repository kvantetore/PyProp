
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
def GetDpDomegaFromFile(filename):
	f = tables.openFile(filename)
	try:
		e = f.root.dpdomega._v_attrs.energy
		th = f.root.dpdomega._v_attrs.theta
		dp = f.root.dpdomega[:]
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


def MakeDpDePlot(energy_double, dpde_double, energy_single, dpde_single, vmax=None):

	rectBound = (.12, .10, .80, .80)
	singleWidth = .025
	singleSpacing = 0.005
	energyLim = (0,14.5)

	fig = figure()
	
	if vmax == None:
		vmax = numpy.max(dpde_double)

	rectMain = (rectBound[0]+singleWidth, rectBound[1]+singleWidth, rectBound[2]-singleWidth, rectBound[3]-singleWidth)
	axMain = fig.add_axes(rectMain)
	axMain.pcolormesh(energy_double, energy_double, dpde_double, vmin=0, vmax=vmax)
	axMain.set_xticks([])
	axMain.set_yticks([])
	axMain.set_xlim(energyLim)
	axMain.set_ylim(energyLim)

	rectLeft = (rectBound[0], rectBound[1]+singleWidth, singleWidth-singleSpacing, rectBound[3]-singleWidth)
	axLeft = fig.add_axes(rectLeft)
	axLeft.pcolormesh(array([0,1]), energy_single, array([dpde_single, dpde_single]).transpose(), vmin=0, vmax=vmax)
	axLeft.set_xticks([])
	axLeft.set_ylim(energyLim)

	rectBottom = (rectBound[0]+singleWidth, rectBound[1], rectBound[2]-singleWidth, singleWidth-singleSpacing)
	axBottom = fig.add_axes(rectBottom)
	axBottom.pcolormesh(energy_single, array([0,1]), array([dpde_single, dpde_single]), vmin=0, vmax=vmax)
	axBottom.set_yticks([])
	axBottom.set_xlim(energyLim)

	draw()

	
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
	partialProb = sum(dp) * diff(energy_double)[0]**2

	return partialProb

def CalculatePartialIonizationProbabilityScan():
	#e0list = [1,5,10,15,30]
	e0list = r_[1:35]

	limits = [0, 4, 9.5, 14.5, 19.5]

	partialList = []
	for e0 in e0list:
		f = tables.openFile("output/stabilization_freq5/stabilization_I_%i_kb20_dt_1e-02.h5" % e0)
		try:
			energy_double = f.root.dpde_double._v_attrs.energy
			dpde_double = f.root.dpde_double[:]
			
			p1 = CalculatePartialIonizationProbability(energy_double, dpde_double, limits[0], limits[1])
			p2 = CalculatePartialIonizationProbability(energy_double, dpde_double, limits[1], limits[2])
			p3 = CalculatePartialIonizationProbability(energy_double, dpde_double, limits[2], limits[3])
			p4 = CalculatePartialIonizationProbability(energy_double, dpde_double, limits[3], limits[4])
			partialList.append([p1, p2, p3, p4])
		finally:
			f.close()
	
	return limits, e0list, array(partialList)

def GetIonizationProbabilityScan():
	e0list = r_[1:35]

	ionList = []
	for e0 in e0list:
		f = tables.openFile("output/stabilization_freq5/stabilization_I_%i_kb20_dt_1e-02.h5" % e0)
		try:
			doubleIon = f.root.wavefunction._v_attrs.DoubleIonization
			singleIon = f.root.wavefunction._v_attrs.SingleIonization
			
			ionList.append([singleIon, doubleIon])
		finally:
			f.close()

	singleList, doubleList = zip(*ionList)
	return e0list, array(singleList), array(doubleList)

