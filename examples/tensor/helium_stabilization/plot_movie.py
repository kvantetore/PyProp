execfile("plot.py")
execfile("plots.py")

def InterpolateRadialDensity(oldGrid, newGrid, radialDensity):
	"""
	Interpolate a 2D radial density from an old mesh (oldGrid, oldGrid) to a new
	mesh (newGrid, newGrid).
	"""
	densityInterp = scipy.interpolate.RectBivariateSpline(oldGrid, oldGrid, radialDensity, kx=3, ky=3)
	densityNew = densityInterp(newGrid, newGrid)
	return densityNew


def CreateIonizationMovieFrame(evalGrid, radialDensity, singleIon, doubleIon, vmax = None):
	"""
	"""

	#Settings
	rmin = evalGrid[0]
	rmax = evalGrid[-1]
	interpGrid = linspace(rmin, rmax, 1000)

	#Figsize settings
	rectBound = (.1, .1, .80, .80)
	spacing = .05
	densityWidth = .8
	densityHeight = .6
	panelWidth = .2

	#load color map
	gradientFile = GetGradientFile()
	cmap = LoadColormap(gradientFile, reverse=True)

	#Set up figure
	fig = figure()
	if vmax == None:
		vmax = numpy.max(radialDensity)
	print "vmax = %s" % (vmax,)

	#Plot radial density
	radialDensityInterp = InterpolateRadialDensity(evalGrid, interpGrid, radialDensity)
	rectDensity = (rectBound[0], rectBound[1], densityWidth, densityHeight)
	axDensity = fig.add_axes(rectDensity)
	axDensity.pcolorfast(interpGrid, interpGrid, radialDensityInterp, vmin=0, vmax=vmax, cmap=cmap)
	#axLeft.set_xticks([])
	#axLeft.set_ylim(energyLim)
	#setp(axLeft.get_yticklines(), "markersize", 2)

	axDensity.plot([0,5], [5,5], 'k-', linewidth=1)
	axDensity.plot([5,5], [0,5], 'k-', linewidth=1)
	axDensity.plot([5,5], [5,25], 'k--', linewidth=.6)
	axDensity.plot([5,25], [5,5], 'k--', linewidth=.6)
	axDensity.text(2.5, 22, "SI", ha="center", va="center")
	axDensity.text(22, 2.5, "SI", ha="center", va="center")
	axDensity.text(22, 22, "DI", ha="center", va="center")

	#Set up info panel
	rectPanel = (rectBound[0], rectBound[1]+densityHeight+spacing, panelWidth, rectBound[3]-densityHeight-spacing)
	axPanel = fig.add_axes(rectPanel)
	axPanel.set_xticks([])
	axPanel.set_yticks([])
	axPanel.bar([1,2,3], [1, singleIon, doubleIon], width=.7, facecolor=UiB_Green, linewidth=.3, align="center")
	axPanel.set_ylim([0,1.1])
	axPanel.set_xlim([.3,3.7])
	axPanel.text(1, 1.3, "B", ha="center", va="center")
	axPanel.text(2, 1.3, "SI", ha="center", va="center")
	axPanel.text(3, 1.3, "DI", ha="center", va="center")

	#Set up laser panel
	rectLaser = (rectBound[0]+panelWidth+spacing, rectBound[1]+densityHeight+spacing, rectBound[2]-panelWidth-spacing, rectBound[3]-densityHeight-spacing)
	axLaser = fig.add_axes(rectLaser)
	axLaser.set_xticks([])
	axLaser.set_yticks([])

	sampleTimes = linspace(0, 2*pi*6, 200)
	laserPulse = sin(pi/sampleTimes[-1]*sampleTimes)**2 * cos(sampleTimes)
	axLaser.plot(sampleTimes, laserPulse, UiB_Green, linewidth=.8)
	axLaser.set_xlim([sampleTimes[0], sampleTimes[-1]])
	axLaser.set_ylim([-1.1, 1.1])
	axLaser.plot([sampleTimes[10]]*2, [laserPulse[10]]*2, "o", color=UiB_Red, markeredgecolor=UiB_Red)

	return fig


def CreateIonizationMovieFrameFromFile(filename, masterFile, radialGrid=None):

	if radialGrid == None:
		radialGrid = linspace(0, 25, 200)
	radialDensity = CreateRadialDensityFromFile(filename, radialGrid)

	CreateIonizationMovieFrame(evalGrid, radialDensity, singleIon, doubleIon, vmax = None)
