execfile("plot_paper_common.py")
execfile("example.py")
execfile("plots.py")
execfile("plot.py")

#Import helium_stabilization_sae example module to get model results
sys.path.append("../")
import helium_stabilization_sae

energyAvgRadius = 1.0

#PaperSetOutput("icpeac")
PaperSetOutput("paper")

InputFolder_1s1s = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s1s"
InputFolder_1s1s_cycle10 = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax90_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s1s_1/"
InputFolder_1s1s_cycle20 = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax160_xsize30_order5_xpartition10_gamma1.0_angular_lmax5_L0-6_M0_1s1s_1_20cycle/"
InputFolder_1s2s = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s2s"
InputFolder_1s2s_1 = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s2s_1"
InputFolder_1s2p_1 = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s2p_1"
InputFolder_1s2p_1 = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s2p_1"
InputFolder_1s2p_3 = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s2p_3"

def GetInputFilename(experiment, e0):
	if isinstance(e0, str):
		return e0

	if experiment == "1s1s":
		return os.path.join(InputFolder_1s1s, "stabilization_I_%i_kb20_dt_1e-02.h5" % e0)
	if experiment == "1s1s_cycle10":
		return os.path.join(InputFolder_1s1s_cycle10, "stabilization_I_%i_kb20_dt_1e-02_T_18.8_phase_zero.h5" % e0)
	if experiment == "1s1s_cycle20":
		return os.path.join(InputFolder_1s1s_cycle20, "stabilization_I_%s_kb20_dt_1e-02_T_37.7_phase_zero.h5" % e0)
	elif experiment == "1s2s":
		return os.path.join(InputFolder_1s2s, "stabilization_I_%i_kb20_dt_1e-02_T_11.3.h5" % e0)
	elif experiment == "1s2s_1":
		return os.path.join(InputFolder_1s2s_1, "stabilization_I_%i_kb20_dt_1e-02_T_11.3_phase_zero.h5" % e0)
	elif experiment == "1s2s_1_cycle10":
		return os.path.join(InputFolder_1s2s_1, "stabilization_I_%i_kb20_dt_1e-02_T_18.8_phase_zero.h5" % e0)
	elif experiment == "1s2p_1":
		return os.path.join(InputFolder_1s2p_1, "stabilization_I_%i_kb20_dt_1e-02_T_11.3_phase_zero.h5" % e0)
	elif experiment == "1s2p_3":
		return os.path.join(InputFolder_1s2p_3, "stabilization_I_%i_kb20_dt_1e-02_T_11.3_phase_zero.h5" % e0)
	else:
		raise Exception("Unknown experiment %s" % experiment)


def GetSymmetry(experiment):
	if experiment == "1s1s" or experiment == "1s1s_cycle10" or experiment == "1s1s_cycle20":
		return "sym"
	elif experiment == "1s2s":
		return "antisym"
	elif experiment == "1s2s_1":
		return "sym"
	elif experiment == "1s2s_1_cycle10":
		return GetSymmetry("1s2s_1")
	elif experiment == "1s2p_3":
		return "antisym"
	elif experiment == "1s2p_1":
		return "sym"
	else:
		raise Exception("Unknown experiment %s" % experiment)


def GetAllExperimentFiles(experiment):
	e0list = r_[1:100]
	fileList = ["%s" % fname for fname in [GetInputFilename(experiment, e0) for e0 in e0list] if os.path.exists(fname)]
	return fileList


def PaperGetBaseEnergy(experiment):
	if experiment == "1s1s" or experiment == "1s1s_cycle10" or experiment == "1s1s_cycle20":
		return 2.903
	elif experiment == "1s2s":
		return 2.17
	elif experiment == "1s2s_1":
		return 2.14
	elif experiment == "1s2s_1_cycle10":
		return PaperGetBaseEnergy("1s2s_1")
	elif experiment == "1s2p_1":
		return 2.13
	elif experiment == "1s2p_3":
		return 2.13
	else:
		raise Exception("Unknown experiment %s" % experiment)


def PaperGetDpDOmegaEnergies(experiment):
	E = array([5*i - PaperGetBaseEnergy(experiment) for i in [1,2,3]])
	if experiment == "1s1s" or experiment == "1s1s_cycle10" or experiment == "1s1s_cycle20":
		return zip(["1photon", "2photon", "3photon", "3photon_inv"], [E[0]/2, E[1]/2, E[1]/2, E[1]/2+5], [E[0]/2, E[1]/2, E[1]/2+5, E[1]/2.0])
	elif experiment == "1s2s":
		return zip(["2photon", "2photon_inv"], [4.55, 3.28], [3.28, 4.55])
	elif experiment == "1s2s_1":
		#return zip(["1photon", "2photon", "3photon"], E/2, [E[0]/2, E[1]/2, E[1]/2])
		return zip(["1photon", "2photon", "2photon_inv", "3photon", "3photon_inv"], [E[0]/2, 3.44, 4.56, 3.44, 9.58], [E[0]/2, 4.56, 3.44, 9.58, 3.44])
		#return [("2photon", 3.93, 3.93)]
	elif experiment == "1s2s_1_cycle10":
		return zip(["1photon", "2photon", "3photon"], [E[0]/2, 3.44, 3.44], [E[0]/2, 4.56, 9.58])
	elif experiment == "1s2p":
		return zip(["2photon", "2photon_inv"], [4.45, 3.42], [3.42, 4.45])
	else:
		raise Exception("Unknown experiment %s" % experiment)


def PaperGetSAEIonization(experiment):
	filename = "output/sae_results/%s.h5" % experiment
	f = tables.openFile(filename, "r")
	try:
		e0 = f.root.Amplitude[:]
		single = f.root.SingleIonization[:]
		double = f.root.DoubleIonization[:]

	finally:
		f.close()

	return e0, single, double


def PaperGetIonization(experiment, e0):
	filename = GetInputFilename(experiment, e0)
	expSym = GetSymmetry(experiment)
	f = tables.openFile(filename, "r")
	try:
		if expSym == "sym":
			single = f.root.wavefunction._v_attrs.SingleIonization
			double = f.root.wavefunction._v_attrs.DoubleIonization
		else:
			single = f.root.wavefunction._v_attrs.AntiSingleIonization
			double = f.root.wavefunction._v_attrs.AntiDoubleIonization

	finally:
		f.close()

	return single, double
		

def PaperMakePlotDpDomega(experiment, e0, energy1, energy2=None, vmin=None, vmax=None):
	if energy2 == None:
		energy2 = energy1

	#get data
	filename = GetInputFilename(experiment, e0)
	en, th, dp = GetDpDomegaDoubleFromFile(filename)
	
	#slice at energies
	idx1 = argmin(abs(en-energy1))
	idx2 = argmin(abs(en-energy2))
	dpSlice = dp[:,:,idx1, idx2]

	#interpolate
	th2 = linspace(0,pi,512)
	dp2 = scipy.interpolate.RectBivariateSpline(th, th, dpSlice)(th2, th2)

	#load color map
	gradientFile = GetGradientFile()
	cmap = LoadColormap(gradientFile, reverse=True)

	#plot
	fig = figure()

	ax = gca()
	ax.pcolorfast(th2/pi, th2/pi, dp2, cmap=cmap, vmin=vmin, vmax=vmax)
	ax.set_xlabel("$\\theta_2$")
	ax.set_ylabel("$\\theta_1$")
	ax.set_xticks([0,0.25,0.5,0.75,1])
	ax.set_xticklabels(("$0$", "$\pi/4$", "$\pi/2$", "$\pi\  3/4$", "$\pi$"))
	ax.set_yticks([0,0.25,0.5,0.75,1])
	ax.set_yticklabels(("$0$", "$\pi/4$", "$\pi/2$", "$\pi\  3/4$", "$\pi$"))
	ax.set_position(Bbox(array([[0.25, 0.25], [0.95, 0.95]])))

	#fig.figurePatch.set_alpha(1.0)
	PaperUpdateFigure(fig)
	ax.set_position(GetOptimalAxesPosition(fig, ax))
	draw()

	return fig


def PaperMakePlotDpDomegaPolarTimeScan(fileList, thCut, energy1, energy2=None):
	if energy2 == None:
		energy2 = energy1

	#setup bounds
	rectBound = (.20, .20, .65, .65)
	singleWidth = .025
	singleSpacing = 0.005
	
	lineStyleList = ["-", "--", ":"]
	colorList = [UiB_Blue, UiB_Green, UiB_Red]

	#create figure
	fig = figure()
	rectMain = (rectBound[0]+singleWidth, rectBound[1]+singleWidth, rectBound[2]-singleWidth, rectBound[3]-singleWidth)
	axMain = fig.add_axes(rectMain, polar=True, axisbelow=True)
	
	maxValDpDomega = 0
	for	filename in fileList:
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
	
		#get current line style and color
		curColor = colorList.pop(0)
		curStyle = lineStyleList.pop(0)

		#normalize
		#dp2Norm = sum(dp2 * sin(th2)) * diff(th2)[0]
		#dp2 /= dp2Norm

		#update max val
		maxValDpDomega = max(max(dp2), maxValDpDomega)

		#plot
		axMain.plot(th2, dp2, color = curColor, linestyle = curStyle)
		axMain.plot(-th2, dp2, color = curColor, linestyle = curStyle)

	
	fig.figurePatch.set_alpha(1.0)
	PaperUpdatePolarFigure(fig)
	#axMain.set_position(GetOptimalAxesPosition(fig, axMain))
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


def PaperMakePlotDpDomegaPolar(experiment, thCut, energy1, energy2=None, vmin=None, vmax=None, energyAverage=False, phiEvalType="avg"):
	if energy2 == None:
		energy2 = energy1

	#setup bounds
	rectBound = (.20, .20, .65, .65)
	singleWidth = .025
	singleSpacing = 0.005
	
	lineStyleList = ["-", "--", ":"]
	colorList = [UiB_Blue, UiB_Green, UiB_Red]
	e0List = [1, 10, 20]
	#e0List = range(1,21)
	#e0List = [15]

	#get photon energy interval
	limits = [0] + [ 5 * (i+0.5) - PaperGetBaseEnergy(experiment) for i in range(1,5)]
	totalEnergy = energy1 + energy2
	enLimLow = filter(lambda en: en<totalEnergy, limits)[-1]
	enLimHigh = filter(lambda en: en>totalEnergy, limits)[0]

	#create figure
	fig = figure()
	rectMain = (rectBound[0]+singleWidth, rectBound[1]+singleWidth, rectBound[2]-singleWidth, rectBound[3]-singleWidth)
	axMain = fig.add_axes(rectMain, polar=True, axisbelow=True)
	
	maxValDpDomega = 0
	for	e0 in e0List:
		#get data
		filename = GetInputFilename(experiment, e0)
		en, th, dp = GetDpDomegaDoubleFromFile(filename, phiEvalType)
		energy_double, dpde_double = GetDpDeDoubleFromFile(filename)

		#slice at energies and theta1
		idx1 = argmin(abs(en-energy1))
		idx2 = argmin(abs(en-energy2))
		idx3 = argmin(abs(th - thCut))
		#dpSlice = sum(sum(dp, axis=3), axis=2)[idx3, :]
		if energyAverage:
			print "Averaging over energy"
			#dpdomega_double = CalculatePartialAngularIonizationProbability(en, th, dp, enLimLow, enLimHigh)
			dpdomega_double = CalculateRadialPartialAngularIonizationProbability(en, th, dp, energy1, energy2, energyAvgRadius)
			dpSlice = dpdomega_double[:, idx3]
		else:
			dpSlice = dp[idx3, : ,idx1, idx2]

		#interpolate
		th2 = linspace(0,pi,512)
		dp2 = scipy.interpolate.UnivariateSpline(th, dpSlice, s=0)(th2)
	
		#get n-photon ionization probability, skip this e0 if too small
		curPhotonProb = CalculatePartialIonizationProbability(energy_double, dpde_double, enLimLow, enLimHigh)
		#get current line style and color
		curColor = colorList.pop(0)
		curStyle = lineStyleList.pop(0)
		#curColor = UiB_Blue
		#if e0 == 10:
		#	curColor = UiB_Green
		#if e0 == 20:
		#	curColor = UiB_Red
		#curStyle = "-"

		#print "Current n-photon prob: %s (lim = %s,%s)" % (curPhotonProb, enLimLow, enLimHigh)
		#if curPhotonProb < 1e-5:
		#	print "Skipping e0 = %s, since n-photon ionization probability is only %s" % (e0,curPhotonProb)
		#	continue

		#normalize
		dp2Norm = sum(dp2 * sin(th2)) * diff(th2)[0]
		dp2 /= dp2Norm

		#update max val
		maxValDpDomega = max(max(dp2), maxValDpDomega)

		#plot
		axMain.plot(th2, dp2, color = curColor, linestyle = curStyle)
		axMain.plot(-th2, dp2, color = curColor, linestyle = curStyle)

	#fig.figurePatch.set_alpha(1.0)
	PaperUpdatePolarFigure(fig)
	axMain.set_position(GetOptimalAxesPosition(fig, axMain))
	axMain.set_xlim([0, maxValDpDomega])

	#Update r and theta grid lines/labels
	#axMain.axes.axesFrame.set_visible(False)
	rGridLines = linspace(0, maxValDpDomega, 5)[1:-1]
	#rGridLines = [maxValDpDomega / 1.618]
	rgrids(rGridLines, ["" for p in rGridLines])
	thetagrids(range(0,360,45), [r"%s$^\circ$" % t for t in range(0,360,45)])

	#put arrow indicating direction of other electron
	axMain.plot([0,0], [0, maxValDpDomega/3.0], "k-", linewidth=0.8)
	axMain.plot([0], [maxValDpDomega/3.0], "k>", markersize=3.0)

	draw()

	return fig


def _PaperMakePlotDpDe(filename, baseEnergy, e0, photonEnergies, normFactor = 1.0, vmax=None):
	#setup bounds
	rectBound = (.15, .15, .80, .80)
	singleWidth = .025
	singleSpacing = 0.005
	energyLim = (0,14.5)

	#get data
	energy_double, dpde_double = GetDpDeDoubleFromFile(filename)
	energy_single, dpde_single = GetDpDeSingleFromFile(filename)

	#interpolate to get smoother plot
	energy_double2 = linspace(energy_double[0], energy_double[-1], 512)
	interp = scipy.interpolate.RectBivariateSpline(energy_double, energy_double, dpde_double)
	dpde_double2 = interp(energy_double2, energy_double2)

	#normalize
	dpde_double2 /= normFactor
	dpde_single /= normFactor
	
	fig = figure()
	
	if vmax == None:
		vmax = numpy.max(dpde_double2)
	#print "vmax = %s" % (vmax,)

	#load color map
	gradientFile = GetGradientFile()
	cmap = LoadColormap(gradientFile, reverse=True)

	#plot double dpde
	rectMain = (rectBound[0]+singleWidth, rectBound[1]+singleWidth, rectBound[2]-singleWidth, rectBound[3]-singleWidth)
	axMain = fig.add_axes(rectMain)
	quadMeshMain = axMain.pcolorfast(energy_double2, energy_double2, dpde_double2, vmin=0, vmax=vmax, cmap=cmap)
	axMain.set_xticks([])
	axMain.set_yticks([])
	axMain.set_xlim(energyLim)
	axMain.set_ylim(energyLim)

	#add photon lines
	for curE in [(5*i - baseEnergy) for i in range(1,4)]:
		lin = Line2D([0,curE,], [curE, 0], color=UiB_Black)
		axMain.add_artist(lin)

	#add circles around photon peaks
	th = linspace(0, 2*pi, 100)
	for photonName, E1, E2 in photonEnergies:
		ear = energyAvgRadius
		circ = Line2D(ear*sin(th) + E1, ear*cos(th) + E2, color=UiB_Black, linestyle="-", linewidth=.1)
		axMain.add_artist(circ)

	#plot left single dpde
	rectLeft = (rectBound[0], rectBound[1]+singleWidth, singleWidth-singleSpacing, rectBound[3]-singleWidth)
	axLeft = fig.add_axes(rectLeft)
	axLeft.pcolorfast(array([0,1]), energy_single, array([dpde_single, dpde_single]).transpose(), vmin=0, vmax=vmax, cmap=cmap)
	axLeft.set_xticks([])
	axLeft.set_ylim(energyLim)
	setp(axLeft.get_yticklines(), "markersize", 2)
	ylabel("Energy (a.u.)", fontsize=10)


	#plot bottom single dpde
	rectBottom = (rectBound[0]+singleWidth, rectBound[1], rectBound[2]-singleWidth, singleWidth-singleSpacing)
	axBottom = fig.add_axes(rectBottom)
	axBottom.pcolorfast(energy_single, array([0,1]), array([dpde_single, dpde_single]), vmin=0, vmax=vmax, cmap=cmap)
	axBottom.set_yticks([])
	axBottom.set_xlim(energyLim)
	setp(axBottom.get_xticklines(), "markersize", 2)
	xlabel("Energy (a.u.)")

	#add colorbar
	cbar = colorbar(quadMeshMain, ax=axMain, format="%3i")
	axColorbar = cbar.ax

	#update figure
	fig.figurePatch.set_alpha(1.0)
	PaperUpdateFigure(fig)

	#update colorbar fonts
	cbarFont = matplotlib.font_manager.FontProperties(size=5.0)
	for ytl in axColorbar.get_yticklabels():
		ytl.set_fontproperties(cbarFont)

	#reposition containers
	RepositionDpDe(fig, axMain, axLeft, axBottom, axColorbar)
	
	draw()

	return fig, axMain, axLeft, axBottom


def _PaperMakePlotDpDeAlt(filename, baseEnergy, e0, photonEnergies, normFactor = 1.0, vmax=None):
	#setup bounds
	rectBound = (.15, .15, .80, .80)
	singleWidth = .025
	singleSpacing = 0.005
	energyLim = (0,14.5)

	#get data
	energy_double, dpde_double = GetDpDeDoubleFromFile(filename)
	energy_single, dpde_single = GetDpDeSingleFromFile(filename)

	#interpolate to get smoother plot
	energy_double2 = linspace(energy_double[0], energy_double[-1], 512)
	interp = scipy.interpolate.RectBivariateSpline(energy_double, energy_double, dpde_double, s=0.0)
	dpde_double2 = interp(energy_double2, energy_double2)

	#normalize
	dpde_double2 /= normFactor
	dpde_single /= normFactor
	
	fig = figure()
	
	if vmax == None:
		vmax = numpy.max(dpde_double2)
	print "vmax = %s" % (vmax,)

	#load color map
	gradientFile = GetGradientFile()
	cmap = LoadColormap(gradientFile, reverse=True)

	#plot double dpde
	rectMain = (rectBound[0], rectBound[1], rectBound[2], rectBound[3])
	axMain = fig.add_axes(rectMain)
	quadMeshMain = axMain.pcolorfast(energy_double2, energy_double2, dpde_double2, vmin=0, vmax=vmax, cmap=cmap)
	#axMain.set_xticks([])
	#axMain.set_yticks([])
	axMain.set_xlim(energyLim)
	axMain.set_ylim(energyLim)
	xlabel("Energy (a.u.)")
	ylabel("Energy (a.u.)", fontsize=10)

	#add photon lines
	for curE in [(5*i - baseEnergy) for i in range(1,4)]:
		lin = Line2D([0,curE,], [curE, 0], color=UiB_Black)
		axMain.add_artist(lin)

	#add circles around photon peaks
	#th = linspace(0, 2*pi, 100)
	#for photonName, E1, E2 in photonEnergies:
	#	ear = energyAvgRadius
	#	circ = Line2D(ear*sin(th) + E1, ear*cos(th) + E2, color=UiB_Black, linestyle="-", linewidth=.1)
	#	axMain.add_artist(circ)

	#add colorbar
	cbar = colorbar(quadMeshMain, ax=axMain, format="%3i")
	axColorbar = cbar.ax

	#update figure
	fig.figurePatch.set_alpha(1.0)
	PaperUpdateFigure(fig)

	#update colorbar fonts
	cbarFont = matplotlib.font_manager.FontProperties(size=5.0)
	for ytl in axColorbar.get_yticklabels():
		ytl.set_fontproperties(cbarFont)

	#reposition containers
	RepositionDpDe(fig, axMain, axLeft, axBottom, axColorbar)
	
	draw()

	return fig, axMain


def PaperMakePlotDpDe(experiment, e0, normFactor = 1.0, vmax=None, plotFunc=_PaperMakePlotDpDe):
	filename = GetInputFilename(experiment, e0)
	baseEnergy = PaperGetBaseEnergy(experiment)
	photonEnergies = PaperGetDpDOmegaEnergies(experiment)
	return plotFunc(filename, baseEnergy, e0, photonEnergies, normFactor = normFactor, vmax=vmax)


def RepositionDpDe(fig, axMain, axLeft, axBottom, axColorbar):
	rectBound = (.15, .15, .80, .80)
	singleWidth = .025
	singleSpacing = 0.005
	energyLim = (0,14.5)

	#reposition
	box = GetOptimalAxesPosition(fig, axMain, [axLeft, axBottom, axColorbar])
	#box = GetOptimalAxesPosition(fig, axMain, padding=0)
	axMain.set_position(box)

	boxLeft = Bbox([[box.xmin - (singleWidth), box.ymin],
	                [box.xmin - (singleSpacing), box.ymax]])
	axLeft.set_position(boxLeft)

	boxBottom = Bbox([[box.xmin, box.ymin - (singleWidth)],
	                 [box.xmax, box.ymin - (singleSpacing)]])
	axBottom.set_position(boxBottom)

	boxColorbar = Bbox([[box.xmax + (3 * singleSpacing), box.ymin],
	                [box.xmax + 0.5, box.ymax]])
	axColorbar.set_position(boxColorbar)
	draw()

	return fig


def PaperMakePlotDpDeFixedEnergy(experiment, e0, photon):
	baseEnergy = PaperGetBaseEnergy(experiment)
	en1 = r_[0.05:photon*5 - baseEnergy:0.05]
	en2 = photon*5-baseEnergy - en1
	en, dp = GetDpDeDoubleFromFile(GetInputFilename(experiment, e0))
	dpInterp = scipy.interpolate.RectBivariateSpline(en, en, dp)
	z = array([dpInterp(x1, y1) for x1, y1 in zip(en1, en2)]).flatten()
	gca().plot(en1/(en1+en2), z)
	maxIdx = argmax(z)
	print "Peak at E1, E2 = %f, %f" % (en1[maxIdx], en2[maxIdx], )
	draw()

#-------------------------------------------------------------------------------------


def PaperMakePlotDpDomegaScan(experiment, doClose=True):
	interactive = rcParams["interactive"]
	rcParams["interactive"] = False

	try:
		e0list = [1,10,20]
		for photonName, E1, E2 in PaperGetDpDOmegaEnergies(experiment):
		
			folder = PaperGetFolder()
			PaperFigureSettings(FigWidth, FigWidth)
	
			for e0 in e0list:
				fig = PaperMakePlotDpDomega(experiment, e0, E1, E2)
				fig.savefig(os.path.join(folder, "%s_dpdomega_e0_%i_%s_equalenergy.eps" % (experiment, e0, photonName)), dpi=300)
				fig.savefig(os.path.join(folder, "%s_dpdomega_e0_%i_%s_equalenergy.pdf" % (experiment, e0, photonName)), dpi=300)
				if doClose: close(fig)
	
	finally:
		rcParams["interactive"] = interactive


def PaperMakePlotDpDomegaPolarScan(experiment, doClose=True, energyAverage=False):
	interactive = rcParams["interactive"]
	rcParams["interactive"] = False
	thCut = 0
	try:
		for photonName, E1, E2 in PaperGetDpDOmegaEnergies(experiment):
		
			folder = PaperGetFolder()
			PaperFigureSettings(FigWidth, FigWidth)

			fig = PaperMakePlotDpDomegaPolar(experiment, thCut, E1, E2, energyAverage=energyAverage, phiEvalType='avg')
			fig2 = PaperMakePlotDpDomegaPolar(experiment, thCut, E1, E2, energyAverage=energyAverage, phiEvalType='coplanar')

			outputPrefix = ""
			if energyAverage:
				outputPrefix = "energyavg_"
			fig.savefig(os.path.join(folder, "%s%s_dpdomega_polar_theta_%s_%s_E1_%2.2f_E2_%2.2f.eps" % (outputPrefix, experiment, thCut, photonName, E1, E2)), dpi=300)
			fig.savefig(os.path.join(folder, "%s%s_dpdomega_polar_theta_%s_%s_E1_%2.2f_E2_%2.2f.pdf" % (outputPrefix, experiment, thCut, photonName, E1, E2)), dpi=300)
			fig2.savefig(os.path.join(folder, "%s%s_dpdomega_polar_theta_%s_%s_E1_%2.2f_E2_%2.2f_coplanar.eps" % (outputPrefix, experiment, thCut, photonName, E1, E2)), dpi=300)
			fig2.savefig(os.path.join(folder, "%s%s_dpdomega_polar_theta_%s_%s_E1_%2.2f_E2_%2.2f_coplanar.pdf" % (outputPrefix, experiment, thCut, photonName, E1, E2)), dpi=300)
			if doClose: close(fig)
	
			if doClose: close(fig)
	
	finally:
		rcParams["interactive"] = interactive


def PaperMakePlotDpDeScan(experiment):
	#e0list = [1,5,10,15,20]
	e0list = [1,10,20]
	folder = PaperGetFolder()
	#single, double = PaperGetIonization(experiment, e0list[0])

	PaperFigureSettings(FigWidth*1.165, FigWidth)

	interactive = rcParams["interactive"]
	rcParams["interactive"] = False

	normFactors = {'1s1s': 0.00072561123476237374/5.0}

	try:
		for e0 in e0list:
			#fig, axMain, ax1, ax2 = PaperMakePlotDpDe(experiment, e0, normFactor = double / 5.0)
			fig, axMain = PaperMakePlotDpDe(experiment, e0, normFactor = normFactors[experiment], plotFunc= _PaperMakePlotDpDeAlt)
			fig.savefig(os.path.join(folder, "%s_dpde_e0_%i.eps" % (experiment, e0,)), dpi=300)
			fig.savefig(os.path.join(folder, "%s_dpde_e0_%i.pdf" % (experiment, e0,)), dpi=300)
			#close(fig)
	
	finally:
		rcParams["interactive"] = interactive

def PaperMakePlotIonization(experiment):
	fileList = GetAllExperimentFiles(experiment)
	symmetry = GetSymmetry(experiment)
	e0, single, double = GetIonizationProbabilityScan(fileList, symmetry)
	#filter
	idx = find(array(e0) <= PaperMaxE0)
	e0, single, double = array(e0)[idx], single[idx], double[idx]

	#add zeros to start and end to get proper closed polys
	e0 = array([0] + list(e0) + [e0[-1]+1])
	single = array([0] + list(single) + [0])
	double = array([0] + list(double) + [0])

	#sae
	e0_sae, single_sae, double_sae = PaperGetSAEIonization(experiment)
	single_sae[:] *= 1.0
	double_sae[:] *= 1.0
	#filter
	idx = find(e0_sae <= PaperMaxE0)
	e0_sae, single_sae, double_sae = e0_sae[idx], single_sae[idx], double_sae[idx]

	inter = isinteractive()
	ioff()
	try:
		PaperFigureSettings(FigWidth*1.165, FigWidth)
		fig = figure()
		ax = fig.gca()
		ax.fill(e0, single+double, facecolor=UiB_Green, linewidth=LineWidth)
		ax.fill(e0, single, facecolor=UiB_Red, linewidth=LineWidth)
		if PaperOutput == PlotOutputIcpeac:
			ax.plot(e0_sae, single_sae+double_sae, "k--")
			ax.plot(e0_sae, single_sae, "k:")
		ax.set_xlim(0,PaperMaxE0)
		ax.set_ylim(0,1)
		ax.set_xlabel("Field Strength (a.u.)")
		ax.set_ylabel("Ionization Probability")
		ax.text(10, .2, "Single ionization")
		ax.text(10, .6, "Double ionization")

		PaperUpdateFigure(fig)
		draw()
		ax.set_position(GetOptimalAxesPosition(fig, ax))
	finally:
		interactive(inter)
	draw()
	folder = PaperGetFolder()
	savefig(os.path.join(folder, "%s_ionization_probability.eps" % experiment), dpi=300)
	savefig(os.path.join(folder, "%s_ionization_probability.pdf" % experiment), dpi=300)


def PaperMakePlotPartialIonization(experiment):
	fileList = GetAllExperimentFiles(experiment)
	limits = [0] + [ 5 * (i+0.5) - PaperGetBaseEnergy(experiment) for i in range(1,5)]
	#print limits
	e0, partial = CalculatePartialIonizationProbabilityScan(fileList,limits)

	#filter
	idx = find(array(e0) <= PaperMaxE0)
	e0, partial = array(e0)[idx], partial[idx,:]

	yMax = (numpy.max(partial) < 0.15 and 0.1) or 1.0

	inter = isinteractive()
	ioff()
	try:
		PaperFigureSettings(FigWidth, FigWidth)
		fig = figure()
		ax = fig.gca()
		ax.plot(e0, sum(partial, axis=1), color=UiB_Blue, label="Total prob.", linestyle="-", linewidth=2*GraphLineWidth)
		ax.plot(e0, partial[:,0], color=UiB_Black, label="1 photon", linestyle="-", linewidth=GraphLineWidth)
		ax.plot(e0, partial[:,1], color=UiB_Green, label="2 photon", linestyle="--", linewidth=GraphLineWidth)
		ax.plot(e0, partial[:,2], color=UiB_Red, label="3 photon", linestyle="-.", linewidth=GraphLineWidth)
		ax.plot(e0, partial[:,3], color=UiB_Blue, label="4 photon", linestyle=":", linewidth=GraphLineWidth)
		legend()
		ax.get_legend().legendPatch.set_visible(False)
		ax.set_xlim(0,PaperMaxE0)
		ax.set_xlabel("Field Strength (a.u.)")
		ax.set_ylabel("Ionization Probability")
		ax.set_ylim(0, yMax)
		
		PaperUpdateFigure(fig)
		ax.set_position(GetOptimalAxesPosition(fig, ax))
	finally:
		interactive(inter)
	draw()

	folder = PaperGetFolder()
	savefig(os.path.join(folder, "%s_ionization_partial_probability.eps" % experiment), dpi=300)
	savefig(os.path.join(folder, "%s_ionization_partial_probability.pdf" % experiment), dpi=300)


def PaperMakePlotPartialIonizationWithModel(experiment):
	fileList = GetAllExperimentFiles(experiment)
	limits = [0] + [ 5 * (i+0.5) - PaperGetBaseEnergy(experiment) for i in range(1,5)]
	#print limits
	e0, partial = CalculatePartialIonizationProbabilityScan(fileList,limits)

	#get SAE data
	curDir = os.getcwd()
	os.chdir("../helium_stabilization_sae")
	e0Model = e0.copy()
	partialModel = helium_stabilization_sae.CalculatePartialIonizationProbabilityScan(experiment, e0Model, limits)
	os.chdir(curDir)

	#filter
	idx = find(array(e0) <= PaperMaxE0)
	e0, partial = array(e0)[idx], partial[idx,:]

	yMax = (numpy.max(partial) < 0.15 and 0.1) or 1.0

	inter = isinteractive()
	ioff()
	try:
		PaperFigureSettings(FigWidth*1.165, FigWidth)
		fig = figure()
		ax = fig.gca()
		ax.plot(e0, sum(partial, axis=1), color=UiB_Blue, label="Total", linestyle="-", linewidth=2*GraphLineWidth)
		ax.plot(e0, partial[:,0], color=UiB_Green, label="1 photon", linestyle="-", linewidth=GraphLineWidth)
		ax.plot(e0, partial[:,1], color=UiB_Black, label="2 photon", linestyle="-", linewidth=GraphLineWidth)
		ax.plot(e0, partial[:,2], color=UiB_Red, label="3 photon", linestyle="-", linewidth=GraphLineWidth)
		#ax.plot(e0, partial[:,3], color=UiB_Blue, label="4 photon", linestyle=":", linewidth=GraphLineWidth)

		#scaleFactor = max(partial[:,1]) / max(partialModel[:,0])
		scaleFactor = max(partial[:,2]) / max(partialModel[:,1])
		#partialModel[:] *= scaleFactor
		partialModel[:] *= 1.5
		saeMarker = ""
		saeLine = ":"
		ax.plot(e0Model, sum(partialModel, axis=1), color=UiB_Blue, label="Total (SAE)", linestyle=saeLine, marker=saeMarker, linewidth=2*GraphLineWidth, markersize=2)
		ax.plot(e0Model, partialModel[:,0], color=UiB_Black, label="2 photon (SAE)", marker=saeMarker, linewidth=GraphLineWidth, linestyle=saeLine, markersize=1)
		ax.plot(e0Model, partialModel[:,1], color=UiB_Red, label="3 photon (SAE)", marker=saeMarker, linestyle=saeLine, linewidth=GraphLineWidth, markersize=1)

		ax.legend(loc=(.5, .48), numpoints = 2)
		ax.get_legend().legendPatch.set_visible(False)
		#cbarFont = matplotlib.font_manager.FontProperties(size=5.0)
		ax.set_xlim(0,PaperMaxE0)
		ax.set_xlabel("Field Strength (a.u.)")
		ax.set_ylabel("Ionization Probability")
		ax.set_ylim(0, yMax)
		
		PaperUpdateFigure(fig)
		draw()
		#ax.set_position(GetOptimalAxesPosition(fig, ax))

		#these are from total ionization figure
		ax.set_position(Bbox(array([[ 0.15806867,  0.16926719],[ 0.95060392,  0.94708333]])))
	finally:
		interactive(inter)
	draw()

	folder = PaperGetFolder()
	savefig(os.path.join(folder, "%s_ionization_partial_probability.eps" % experiment), dpi=300)
	savefig(os.path.join(folder, "%s_ionization_partial_probability.pdf" % experiment), dpi=300)

	print "Scaling factor for SAE model results: %s" % scaleFactor


def PaperMakePlotsExperiment(experiment):
	#PaperMakePlotDpDomegaScan(experiment)
	PaperMakePlotDpDomegaPolarScan(experiment)
	PaperMakePlotDpDeScan(experiment)
	#PaperMakePlotIonization(experiment)
	PaperMakePlotPartialIonization(experiment)

def PaperMakePlots():
	"""
	Make all plots for paper
	"""
	PaperMakePlotsExperiment("1s1s")
	PaperMakePlotsExperiment("1s1s_cycle10")
	PaperMakePlotsExperiment("1s2s")
	PaperMakePlotsExperiment("1s2s_1")
	PaperMakePlotsExperiment("1s2s_1_cycle10")
	#PaperMakePlotsExperiment("1s2p")

