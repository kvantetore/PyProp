execfile("plot_paper_common.py")
execfile("example.py")



#------------------ Run Stabilization -------------------------

def PaperGetExperimentArgs(experiment, e0):
	"""
	Sets up two argument dicts for propagating a product of 
	two 1e systems. I.e. if experiment == "1s2p",
	the first argument dict corresponds to a 1s he+ system
	and the second argument dict corresponds to a 2p h system

	if experiment == "1s1s" two idential arg dicts are returned
	containing a 1s he2 (Z=1.7) model
	"""
	freq = 5.
	arg = {
		"frequency": freq,
		"amplitude": e0/freq,
		"config": "config_paper.ini",
		"experiment": experiment,
		"e0": e0,
		"initialN": 1,
		"initialL": 0,
		"outputCount": 10,
		}

	if experiment == "1s1s":
		arg["model"] = "he2"
		arg1 = arg2 = arg

	elif experiment == "1s2s":
		arg1 = dict(arg)
		arg2 = dict(arg)
		arg1["model"] = "he+"
		arg2["model"] = "h"
		arg2["initialN"] = 2
		arg2["initialL"] = 0

	elif experiment == "1s2s_2":
		arg1 = dict(arg)
		arg2 = dict(arg)
		arg1["model"] = "he+"
		arg2["model"] = "he"
		arg2["initialN"] = 2
		arg2["initialL"] = 0

	elif experiment == "1s2p":
		arg1 = dict(arg)
		arg2 = dict(arg)
		arg1["model"] = "he+"
		arg2["model"] = "h"
		arg2["initialN"] = 2
		arg2["initialL"] = 1

	elif experiment == "1s2p_2":
		arg1 = dict(arg)
		arg2 = dict(arg)
		arg1["model"] = "he+"
		arg2["model"] = "he"
		arg2["initialN"] = 2
		arg2["initialL"] = 1
	

	else:
		raise Exception("Experiment %s not found" % experiment)

	return arg1, arg2


if not "ComputedProblems" in globals():
	ComputedProblems = {}
def PaperGetProblems(experiment, e0):
	global ComputedProblems
	key = "%s_%f" % (experiment, e0)
	if not key in ComputedProblems:
		#get arguments for this experiment
		arg1, arg2 = PaperGetExperimentArgs(experiment, e0)
		print arg1

		#do stabilization runs for this experiment
		prop1 = RunStabilization(**arg1)
		if arg2 == arg1:
			prop2 = prop1
		else:
			prop2 = RunStabilization(**arg2)
		ComputedProblems[key] = (prop1, prop2)

	return ComputedProblems[key]

if not "AnalysisData" in globals():
	AnalysisData = {}
def GetAnalysisData(experiment):
	global AnalysisData
	key = "%s" % (experiment, )
	if not key in AnalysisData:
		#get arguments for this experiment
		arg1, arg2 = PaperGetExperimentArgs(experiment, 0)

		def getAnalysisData(**args):
			prop = SetupProblem(**args)
			S = SetupOverlapMatrix(prop)
			E, V = SetupRadialEigenstates(prop)
			return S, E, V

		AnalysisData[key] = (getAnalysisData(**arg1), getAnalysisData(**arg2))

	return AnalysisData[key]


def PaperGetDpDeDouble(experiment, e0):
	#x = r_[:10:0.1]
	#return x, outer(sin(x), sin(x))

	prop1, prop2 = PaperGetProblems(experiment, e0)
	(S1, E1, V1), (S2, E2, V2) = GetAnalysisData(experiment)

	en1, dp1 = CalculateEnergyDistribution(prop1.psi, E1, V1, S1, 0.01)
	en2, dp2 = CalculateEnergyDistribution(prop1.psi, E2, V2, S2, 0.01)

	return en1, en2, outer(dp1,dp2)


def PaperGetDpDeSingle(experiment, e0):
	#x = r_[:10:0.1]
	#return x, sin(x)

	prop1, prop2 = PaperGetProblems(experiment, e0)
	(S1, E1, V1), (S2, E2, V2) = GetAnalysisData(experiment)

	en, dp = CalculateEnergyDistribution(prop1.psi, E1, V1, S1, 0.005)
	return en, dp
	 

def GetIonizationProbabilityScan(experiment):
	e0List = r_[1:35:1]
	
	def getIonizationProbability(e0):
		prop1, prop2 = PaperGetProblems(experiment, e0)
		Ion1 = 1-prop1.BoundTotalList[-1]
		Ion2 = 1-prop2.BoundTotalList[-1]
		print Ion1, Ion2
		single = Ion1*(1-Ion2) + Ion2*(1-Ion1)
		double = Ion1*Ion2
		return single, double

	singleList, doubleList = zip(*map(getIonizationProbability, e0List))
	return e0List, singleList, doubleList


def SaveIonizationResults(experiment):
	e0, single, double = GetIonizationProbabilityScan(experiment)
	folder = PaperGetFolder()
	f = tables.openFile(os.path.join(folder, "%s.h5" % experiment), "w")
	try:
		f.createArray(f.root, "Amplitude", array(e0))
		f.createArray(f.root, "SingleIonization", array(single))
		f.createArray(f.root, "DoubleIonization", array(double))

	finally:
		f.close()

	

#------------------ Run Stabilization -------------------------

#def PaperMakePlotDpDomega(e0, energy1, energy2=None, vmin=None, vmax=None):
#	if energy2 == None:
#		energy2 = energy1
#
#	#get data
#	if type(e0) == str:
#		filename = e0
#	else:
#		filename = "output/stabilization_freq5/stabilization_I_%i_kb20_dt_1e-02.h5" % e0
#	en, th, dp = GetDpDomegaDoubleFromFile(filename)
#	
#	#slice at energies
#	idx1 = argmin(abs(en-energy1))
#	idx2 = argmin(abs(en-energy2))
#	dpSlice = dp[:,:,idx1, idx2]
#
#	#interpolate
#	th2 = linspace(0,pi,512)
#	dp2 = scipy.interpolate.RectBivariateSpline(th, th, dpSlice)(th2, th2)
#
#	#load color map
#	cmap = LoadColormap(GradientFile, reverse=True)
#
#	#plot
#	fig = figure()
#
#	ax = gca()
#	ax.pcolorfast(th2/pi, th2/pi, dp2, cmap=cmap, vmin=vmin, vmax=vmax)
#	ax.set_xlabel("$\\theta_1$")
#	ax.set_ylabel("$\\theta_2$")
#	xticks([0,0.25,0.5,0.75,1], ("$0$", "$\pi/4$", "$\pi/2$", "$\pi\  3/4$", "$\pi$"))
#	yticks([0,0.25,0.5,0.75,1], ("$0$", "$\pi/4$", "$\pi/2$", "$\pi\  3/4$", "$\pi$"))
#	ax.set_position(Bbox(array([[0.25, 0.25], [0.95, 0.95]])))
#	draw()
#
#
#	return fig
#

def PaperMakePlotDpDe(experiment, e0, vmax=None):
	"""

	"""

	#setup bounds
	rectBound = (.15, .15, .80, .80)
	singleWidth = .025
	singleSpacing = 0.005
	energyLim = (0,14.5)

	#get data
	energy_double1, energy_double2, dpde_double = PaperGetDpDeDouble(experiment, e0)
	energy_single, dpde_single = PaperGetDpDeSingle(experiment, e0)

	#interpolate to get smoother plot
	energy_interp = linspace(energy_double1[0], energy_double1[-1], 512)
	interp = scipy.interpolate.RectBivariateSpline(energy_double1, energy_double2, dpde_double)
	dpde_double_interp = interp(energy_interp, energy_interp)
	
	fig = figure()
	
	if vmax == None:
		vmax = numpy.max(dpde_double)
	print "vmax = %s" % (vmax,)

	#load color map
	cmap = LoadColormap(GradientFile, reverse=True)

	#plot double dpde
	rectMain = (rectBound[0]+singleWidth, rectBound[1]+singleWidth, rectBound[2]-singleWidth, rectBound[3]-singleWidth)
	axMain = fig.add_axes(rectMain)
	axMain.pcolorfast(energy_interp, energy_interp, dpde_double_interp, vmin=0, vmax=vmax, cmap=cmap)
	axMain.set_xticks([])
	axMain.set_yticks([])
	axMain.set_xlim(energyLim)
	axMain.set_ylim(energyLim)

	for curE in [(5*i - 2.9) for i in range(1,4)]:
	#for curE in [(5*i - 2.17) for i in range(1,4)]:
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
	draw()

	fig.figurePatch.set_alpha(1.0)
	PaperUpdateFigure(fig)
	RepositionDpDe(fig, axMain, axLeft, axBottom)

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


	

#-------------------------------------------------------------------------------------

#
#def PaperMakePlotDpDomegaScan(doClose=True):
#	e0list = [1,15,30]
#	E1 = (5 - 2.9)/2
#	E2 = (2*5 - 2.9)/2
#
#	folder = PaperGetFolder()
#	
#
#	interactive = rcParams["interactive"]
#	rcParams["interactive"] = False
#	try:
#		for e0 in e0list:
#			PaperFigureSettings(FigWidth, FigWidth)
#			fig = PaperMakePlotDpDomega(e0, E1, E1)
#			fig.figurePatch.set_alpha(1.0)
#			PaperUpdateFigure(fig)
#			fig.savefig(os.path.join(folder, "dpdomega_e0_%i_1photon_equalenergy.eps" % (e0,)), dpi=300)
#			fig.savefig(os.path.join(folder, "dpdomega_e0_%i_1photon_equalenergy.pdf" % (e0,)), dpi=300)
#			if doClose: close(fig)
#		
#			PaperFigureSettings(FigWidth, FigWidth)
#			fig = PaperMakePlotDpDomega(e0, E2, E2)
#			fig.figurePatch.set_alpha(1.0)
#			PaperUpdateFigure(fig)
#			fig.savefig(os.path.join(folder, "dpdomega_e0_%i_2photon_equalenergy.eps" % (e0,)), dpi=300)
#			fig.savefig(os.path.join(folder, "dpdomega_e0_%i_2photon_equalenergy.pdf" % (e0,)), dpi=300)
#			if doClose: close(fig)
#	
#	finally:
#		rcParams["interactive"] = interactive
#

def PaperMakePlotDpDeScan(experiment):
	e0list = [1,5,10,20]
	E1 = (5 - 2.9)/2
	E2 = (2*5 - 2.9)/2

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
	e0, single, double = GetIonizationProbabilityScan(experiment)

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
		ax.axis([0,35,0,1])
		ax.set_xlabel("Field Strength (a.u.)")
		ax.set_ylabel("Ionization Probability")
	
		ax.set_position(Bbox(array([[0.21, 0.21], [0.92, 0.92]])))
	
		PaperUpdateFigure(fig)
	finally:
		interactive(inter)
	
	draw()
	folder = PaperGetFolder()
	savefig(os.path.join(folder, "%s_ionization_probability.eps" % (experiment, )), dpi=300)
	savefig(os.path.join(folder, "%s_ionization_probability.pdf" % (experiment, )), dpi=300)


#def PaperMakePlotPartialIonization():
#	limits, e0, partial = CalculatePartialIonizationProbabilityScan()
#	
#	inter = isinteractive()
#	ioff()
#	try:
#		PaperFigureSettings()
#		fig = figure()
#		ax = fig.gca()
#		ax.plot(e0, partial[:,0], color=UiB_Black, label="1 photon", linestyle="-")
#		ax.plot(e0, partial[:,1], color=UiB_Green, label="2 photon", linestyle="--")
#		ax.plot(e0, partial[:,2], color=UiB_Red, label="3 photon", linestyle="-.")
#		ax.plot(e0, partial[:,3], color=UiB_Blue, label="4 photon", linestyle=":")
#		legend()
#		ax.get_legend().legendPatch.set_visible(False)
#		ax.set_xlim(0,35)
#		ax.set_xlabel("Field Strength (a.u.)")
#		ax.set_ylabel("Ionization Probability")
#		ax.set_ylim(0,1)
#		
#		ax.set_position(Bbox(array([[0.21, 0.21], [0.92, 0.92]])))
#		PaperUpdateFigure(fig)
#	finally:
#		interactive(inter)
#	draw()
#
#	folder = PaperGetFolder()
#	savefig(os.path.join(folder, "%s_ionization_partial_probability.eps" % experiment), dpi=300)
#	savefig(os.path.join(folder, "%s_ionization_partial_probability.pdf" % experiment), dpi=300)
#
def PaperMakePlots():
	"""
	Make all plots for paper
	"""
	PaperMakePlotDpDomegaScan()
	PaperMakePlotDpDeScan()
	PaperMakePlotIonization()
	PaperMakePlotPartialIonization()


