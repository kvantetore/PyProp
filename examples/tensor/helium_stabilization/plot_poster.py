from matplotlib.transforms import Bbox

execfile("example.py")
execfile("plots.py")
execfile("plot.py")

UiB_Orange = "#d95900"
UiB_Green  = "#77af00"
UiB_Red    = "#aa0000"
UiB_Blue   = "#005473"
UiB_Black  = "#000000"


SandbjergFigureFont = matplotlib.font_manager.FontProperties(family="Times New Roman", size=16)
def SandbjergGetFont():
	return SandbjergFigureFont

def SandbjergSetAllFonts(parent):
	if hasattr(parent, "get_children"):
		for chld in parent.get_children():
			if hasattr(chld, "set_fontproperties"):
				chld.set_fontproperties(SandbjergGetFont())
			SandbjergSetAllFonts(chld)
	

def SandbjergFigureSettings(fig_width=16, fig_height=12):
	cm_to_inch = 0.393700787
	fig_width *= cm_to_inch
	fig_height *= cm_to_inch
	fig_size =  [fig_width,fig_height]
	params = {
			  'axes.labelsize': 16,
			  'text.fontsize': 16,
			  'legend.fontsize': 14,
			  'xtick.labelsize': 14,
			  'ytick.labelsize': 14,
			  'text.usetex': False,
			  'figure.figsize': fig_size,
			  'font.family' : 'serif',
			  'font.serif' : 'Times New Roman'}
	matplotlib.rcParams.update(params)

def SandbjergUpdateFigure(fig):
	SandbjergSetAllFonts(fig)
	fig.figurePatch.set_visible(False)
	for ax in fig.axes:
		ax.axesPatch.set_visible(False)
	draw()


def SandbjergGetFolder():
	folder = "figs/sandbjerg"
	if not os.path.exists(folder):
		os.makedirs(folder)
	return folder



def SandbjergMakePlotDpDomega(e0, energy1, energy2=None, vmin=None, vmax=None):
	if energy2 == None:
		energy2 = energy1

	#get data
	if type(e0) == str:
		filename = e0
	else:
		filename = "output/stabilization_freq5/stabilization_I_%i_kb20_dt_1e-02.h5" % e0
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

	pcolormesh(th2/pi, th2/pi, dp2, cmap=cmap, vmin=vmin, vmax=vmax)
	xticks([0,0.25,0.5,0.75,1], ("$0$", "$\pi/4$", "$\pi/2$", "$\pi\  3/4$", "$\pi$"))
	yticks([0,0.25,0.5,0.75,1], ("$0$", "$\pi/4$", "$\pi/2$", "$\pi\  3/4$", "$\pi$"))


	return fig


def SandbjergMakePlotDpDe(e0, vmax=None):
	"""

	"""

	#setup bounds
	rectBound = (.15, .15, .80, .80)
	singleWidth = .025
	singleSpacing = 0.005
	energyLim = (0,14.5)

	#get data
	if type(e0) == str:
		filename = e0
	else:
		filename = "output/stabilization_freq5/stabilization_I_%i_kb20_dt_1e-02.h5" % e0
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
	axMain.pcolormesh(energy_double2, energy_double2, dpde_double2, vmin=0, vmax=vmax, cmap=cmap)
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
	axLeft.pcolormesh(array([0,1]), energy_single, array([dpde_single, dpde_single]).transpose(), vmin=0, vmax=vmax, cmap=cmap)
	axLeft.set_xticks([])
	axLeft.set_ylim(energyLim)
	ylabel("Energy (a.u.)")


	#plot bottom single dpde
	rectBottom = (rectBound[0]+singleWidth, rectBound[1], rectBound[2]-singleWidth, singleWidth-singleSpacing)
	axBottom = fig.add_axes(rectBottom)
	axBottom.pcolormesh(energy_single, array([0,1]), array([dpde_single, dpde_single]), vmin=0, vmax=vmax, cmap=cmap)
	axBottom.set_yticks([])
	axBottom.set_xlim(energyLim)
	xlabel("Energy (a.u.)")

	draw()

	return fig


	

#-------------------------------------------------------------------------------------


def SandbjergMakePlotDpDomegaScan():
	e0list = [1,15,30]
	E1 = (5 - 2.9)/2
	E2 = (2*5 - 2.9)/2

	folder = SandbjergGetFolder()
	

	interactive = rcParams["interactive"]
	rcParams["interactive"] = False
	try:
		for e0 in e0list:
			SandbjergFigureSettings(11.5, 11.5)
			fig = SandbjergMakePlotDpDomega(e0, E1, E1)
			fig.figurePatch.set_alpha(1.0)
			SandbjergUpdateFigure(fig)
			fig.savefig(os.path.join(folder, "dpdomega_e0_%i_1photon_equalenergy.png" % (e0,)), dpi=300)
			close(fig)
		
			SandbjergFigureSettings(11.5, 11.5)
			fig = SandbjergMakePlotDpDomega(e0, E2, E2)
			fig.figurePatch.set_alpha(1.0)
			SandbjergUpdateFigure(fig)
			fig.savefig(os.path.join(folder, "dpdomega_e0_%i_2photon_equalenergy.png" % (e0,)), dpi=300)
			close(fig)
	
	finally:
		rcParams["interactive"] = interactive


def SandbjergMakePlotDpDeScan():
	e0list = [1,15,30]
	E1 = (5 - 2.9)/2
	E2 = (2*5 - 2.9)/2

	folder = SandbjergGetFolder()

	SandbjergFigureSettings(11.5, 11.5)

	interactive = rcParams["interactive"]
	rcParams["interactive"] = False
	try:
		for e0 in e0list:
			fig = SandbjergMakePlotDpDe(e0)
			fig.figurePatch.set_alpha(1.0)
			SandbjergUpdateFigure(fig)
			fig.savefig(os.path.join(folder, "dpde_e0_%i.png" % (e0,)), dpi=300)
			#close(fig)
	
	finally:
		rcParams["interactive"] = interactive


def SandbjergMakePlotIonization():
	e0, single, double = GetIonizationProbabilityScan()

	#add zeros to start and end to get proper closed polys
	e0 = array([0] + list(e0) + [e0[-1]+1])
	single = array([0] + list(single) + [0])
	double = array([0] + list(double) + [0])

	SandbjergFigureSettings(16, 12)
	fig = figure()
	ax = fig.gca()
	ax.fill(e0, single+double, facecolor=UiB_Green)
	ax.fill(e0, single, facecolor=UiB_Red)
	ax.axis([0,35,0,1])
	ax.set_xlabel("Field Strength (a.u.)")
	ax.set_ylabel("Ionization Probability")

	ax.set_position(Bbox(array([[0.125, 0.125], [0.9, 0.9]])))

	SandbjergUpdateFigure(fig)
	draw()
	folder = SandbjergGetFolder()
	savefig(os.path.join(folder, "ionization_probability.png"), dpi=300)
	savefig(os.path.join(folder, "ionization_probability.svg"))


def SandbjergMakePlotPartialIonization():
	limits, e0, partial = CalculatePartialIonizationProbabilityScan()
	
	folder = SandbjergGetFolder()

	SandbjergFigureSettings(16, 12)
	fig = figure()
	ax = fig.gca()
	ax.plot(e0, partial[:,0], color=UiB_Black, label="1 photon")
	ax.plot(e0, partial[:,1], color=UiB_Green, label="2 photon")
	ax.plot(e0, partial[:,2], color=UiB_Red, label="3 photon")
	ax.plot(e0, partial[:,3], color=UiB_Blue, label="4 photon")
	legend()
	ax.get_legend().legendPatch.set_visible(False)
	ax.set_xlim(0,35)
	ax.set_xlabel("Field Strength (a.u.)")
	ax.set_ylabel("Ionization Probability")
	ax.set_ylim(0,1)
	
	ax.set_position(Bbox(array([[0.125, 0.125], [0.9, 0.9]])))
	
	SandbjergUpdateFigure(fig)
	folder = SandbjergGetFolder()
	savefig(os.path.join(folder, "ionization_partial_probability.png"), dpi=300)
	savefig(os.path.join(folder, "ionization_partial_probability.svg"))

def SandbjergMakePlots():
	"""
	Make all plots for sandbjerg poster
	"""
	SandbjergMakePlotDpDomegaScan()
	SandbjergMakePlotDpDeScan()
	SandbjergMakePlotIonization()
	SandbjergMakePlotPartialIonization()

