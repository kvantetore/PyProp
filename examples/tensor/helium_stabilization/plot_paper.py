execfile("plot_paper_common.py")
execfile("example.py")
execfile("plots.py")
execfile("plot.py")

InputFolder_1s1s = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s1s"
InputFolder_1s2s = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s2s"
InputFolder_1s2p = "output/stabilization/freq_5.0_grid_exponentiallinear_xmax80_xsize80_order5_xpartition20_gamma2.0_angular_lmax5_L0-6_M0_1s2p"

PlotToScreen = False
if PlotToScreen:
	FigWidth *= 2
	FigWidthLarge *= 2
	LineWidth = 1
	GraphLineWidth = 1.5

def GetInputFilename(experiment, e0):
	if isinstance(e0, str):
		return e0

	if experiment == "1s1s":
		return os.path.join(InputFolder_1s1s, "stabilization_I_%i_kb20_dt_1e-02.h5" % e0)
	elif experiment == "1s2s":
		return os.path.join(InputFolder_1s2s, "stabilization_I_%i_kb20_dt_1e-02_T_11.3.h5" % e0)
	elif experiment == "1s2p":
		return os.path.join(InputFolder_1s2p, "stabilization_I_%i_kb20_dt_1e-02_T_11.3.h5" % e0)
	else:
		raise Exception("Unknown experiment %s" % experiment)


def GetSymmetry(experiment):
	if experiment == "1s1s":
		return "sym"
	elif experiment == "1s2s":
		return "antisym"
	elif experiment == "1s2p":
		return "antisym"
	else:
		raise Exception("Unknown experiment %s" % experiment)


def GetAllExperimentFiles(experiment):
	e0list = r_[1:100]
	fileList = ["%s" % fname for fname in [GetInputFilename(experiment, e0) for e0 in e0list] if os.path.exists(fname)]
	return fileList


def PaperGetBaseEnergy(experiment):
	if experiment == "1s1s":
		return 2.903
	elif experiment == "1s2s":
		return 2.17
	elif experiment == "1s2p":
		return 2.13
	else:
		raise Exception("Unknown experiment %s" % experiment)


def PaperGetDpDOmegaEnergies(experiment):
	E = array([5*i - PaperGetBaseEnergy(experiment) for i in [1,2]])
	if experiment == "1s1s":
		return zip(["1photon", "2photon"], E/2, E/2)
	elif experiment == "1s2s":
		return [("2photon", 4.55, 3.28)]
	elif experiment == "1s2p":
		return [("2photon", 4.45, 3.42)]
	else:
		raise Exception("Unknown experiment %s" % experiment)
		

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


def PaperMakePlotDpDe(experiment, e0, vmax=None):
	"""

	"""

	#setup bounds
	rectBound = (.15, .15, .80, .80)
	singleWidth = .025
	singleSpacing = 0.005
	energyLim = (0,14.5)

	#get data
	filename = GetInputFilename(experiment, e0)
	energy_double, dpde_double = GetDpDeDoubleFromFile(filename)
	energy_single, dpde_single = GetDpDeSingleFromFile(filename)

	#interpolate to get smoother plot
	energy_double2 = linspace(energy_double[0], energy_double[-1], 512)
	interp = scipy.interpolate.RectBivariateSpline(energy_double, energy_double, dpde_double)
	dpde_double2 = interp(energy_double2, energy_double2)
	
	fig = figure()
	
	if vmax == None:
		vmax = numpy.max(dpde_double)
	print "vmax = %s" % (vmax,)

	#load color map
	gradientFile = GetGradientFile()
	cmap = LoadColormap(gradientFile, reverse=True)

	#plot double dpde
	rectMain = (rectBound[0]+singleWidth, rectBound[1]+singleWidth, rectBound[2]-singleWidth, rectBound[3]-singleWidth)
	axMain = fig.add_axes(rectMain)
	axMain.pcolorfast(energy_double2, energy_double2, dpde_double2, vmin=0, vmax=vmax, cmap=cmap)
	axMain.set_xticks([])
	axMain.set_yticks([])
	axMain.set_xlim(energyLim)
	axMain.set_ylim(energyLim)

	baseEnergy = PaperGetBaseEnergy(experiment)
	for curE in [(5*i - baseEnergy) for i in range(1,4)]:
		lin = Line2D([0,curE,], [curE, 0], color=UiB_Black)
		axMain.add_artist(lin)

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

	fig.figurePatch.set_alpha(1.0)
	PaperUpdateFigure(fig)
	RepositionDpDe(fig, axMain, axLeft, axBottom)
	draw()

	return fig, axMain, axLeft, axBottom


def RepositionDpDe(fig, axMain, axLeft, axBottom):
	rectBound = (.15, .15, .80, .80)
	singleWidth = .025
	singleSpacing = 0.005
	energyLim = (0,14.5)

	#reposition
	box = GetOptimalAxesPosition(fig, axMain, [axMain, axLeft, axBottom])
	#box = GetOptimalAxesPosition(fig, axMain, padding=0)
	axMain.set_position(box)

	boxLeft = Bbox([[box.xmin - (singleWidth), box.ymin],
	                [box.xmin - (singleSpacing), box.ymax]])
	axLeft.set_position(boxLeft)

	boxBottom = Bbox([[box.xmin, box.ymin - (singleWidth)],
	                 [box.xmax, box.ymin - (singleSpacing)]])
	axBottom.set_position(boxBottom)
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
			PaperFigureSettings(FigWidth+0.52, FigWidth)
	
			for e0 in e0list:
				fig = PaperMakePlotDpDomega(experiment, e0, E1, E2)
				fig.savefig(os.path.join(folder, "%s_dpdomega_e0_%i_%s_equalenergy.eps" % (experiment, e0, photonName)), dpi=300)
				fig.savefig(os.path.join(folder, "%s_dpdomega_e0_%i_%s_equalenergy.pdf" % (experiment, e0, photonName)), dpi=300)
				if doClose: close(fig)
	
	finally:
		rcParams["interactive"] = interactive


def PaperMakePlotDpDeScan(experiment):
	e0list = [1,5,10,15,20]
	folder = PaperGetFolder()

	PaperFigureSettings(FigWidth, FigWidth)

	interactive = rcParams["interactive"]
	rcParams["interactive"] = False
	try:
		for e0 in e0list:
			fig, axMain, ax1, ax2 = PaperMakePlotDpDe(experiment, e0)
			fig.savefig(os.path.join(folder, "%s_dpde_e0_%i.eps" % (experiment, e0,)), dpi=300)
			fig.savefig(os.path.join(folder, "%s_dpde_e0_%i.pdf" % (experiment, e0,)), dpi=300)
			#close(fig)
	
	finally:
		rcParams["interactive"] = interactive

def PaperMakePlotIonization(experiment):
	fileList = GetAllExperimentFiles(experiment)
	symmetry = GetSymmetry(experiment)
	e0, single, double = GetIonizationProbabilityScan(fileList, symmetry)

	#add zeros to start and end to get proper closed polys
	e0 = array([0] + list(e0) + [e0[-1]+1])
	single = array([0] + list(single) + [0])
	double = array([0] + list(double) + [0])

	inter = isinteractive()
	ioff()
	try:
		PaperFigureSettings(FigWidthLarge, FigWidthLarge/1.33)
		fig = figure()
		ax = fig.gca()
		ax.fill(e0, single+double, facecolor=UiB_Green, linewidth=LineWidth)
		ax.fill(e0, single, facecolor=UiB_Red, linewidth=LineWidth)
		ax.set_xlim(0,35)
		ax.set_ylim(0,1)
		ax.set_xlabel("Field Strength (a.u.)")
		ax.set_ylabel("Ionization Probability")
	
		PaperUpdateFigure(fig)
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
	print limits
	e0, partial = CalculatePartialIonizationProbabilityScan(fileList,limits)

	inter = isinteractive()
	ioff()
	try:
		PaperFigureSettings(FigWidthLarge, FigWidthLarge/1.33)
		fig = figure()
		ax = fig.gca()
		ax.plot(e0, partial[:,0], color=UiB_Black, label="1 photon", linestyle="-", linewidth=GraphLineWidth)
		ax.plot(e0, partial[:,1], color=UiB_Green, label="2 photon", linestyle="--", linewidth=GraphLineWidth)
		ax.plot(e0, partial[:,2], color=UiB_Red, label="3 photon", linestyle="-.", linewidth=GraphLineWidth)
		ax.plot(e0, partial[:,3], color=UiB_Blue, label="4 photon", linestyle=":", linewidth=GraphLineWidth)
		legend()
		ax.get_legend().legendPatch.set_visible(False)
		ax.set_xlim(0,35)
		ax.set_xlabel("Field Strength (a.u.)")
		ax.set_ylabel("Ionization Probability")
		ax.set_ylim(0,1)
		
		PaperUpdateFigure(fig)
		ax.set_position(GetOptimalAxesPosition(fig, ax))
	finally:
		interactive(inter)
	draw()

	folder = PaperGetFolder()
	savefig(os.path.join(folder, "%s_ionization_partial_probability.eps" % experiment), dpi=300)
	savefig(os.path.join(folder, "%s_ionization_partial_probability.pdf" % experiment), dpi=300)

def PaperMakePlotsExperiment(experiment):
	PaperMakePlotDpDomegaScan(experiment)
	PaperMakePlotDpDeScan(experiment)
	PaperMakePlotIonization(experiment)
	PaperMakePlotPartialIonization(experiment)

def PaperMakePlots():
	"""
	Make all plots for paper
	"""
	PaperMakePlotsExperiment("1s1s")
	PaperMakePlotsExperiment("1s2s")
	PaperMakePlotsExperiment("1s2p")

