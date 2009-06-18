execfile("plot_paper_common.py")
execfile("example.py")
execfile("plot.py")
execfile("plots.py")
from pyprop import DefaultCoupledIndexIterator

def InterpolateRadialDensity(oldGrid, newGrid, radialDensity):
	"""
	Interpolate a 2D radial density from an old mesh (oldGrid, oldGrid) to a new
	mesh (newGrid, newGrid).
	"""
	densityInterp = scipy.interpolate.RectBivariateSpline(oldGrid, oldGrid, radialDensity, kx=3, ky=3)
	densityNew = densityInterp(newGrid, newGrid)
	return densityNew


def CreateIonizationMovieFrame(evalGrid, radialDensity, singleIon, doubleIon, timeIdx, sampleTimes, conf, vmin = None, vmax = None):
	"""
	"""

	#Settings
	rmin = evalGrid[0]
	rmax = evalGrid[-1]
	interpGrid = linspace(rmin, rmax, 1000)
	singleIonBoxSize = 10

	#Figsize settings
	rectBound = (.1, .1, .80, 0.8)
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
	if vmin == None:
		vmin = 0.0
		print "vmin = %s" % (vmin,)

	#Plot radial density
	radialDensityInterp = InterpolateRadialDensity(evalGrid, interpGrid, radialDensity)
	rectDensity = (rectBound[0], rectBound[1], densityWidth, densityHeight)
	axDensity = fig.add_axes(rectDensity)
	axDensity.pcolorfast(interpGrid, interpGrid, log(radialDensityInterp), vmin=vmin, vmax=vmax, cmap=cmap)
	#axLeft.set_xticks([])
	#axLeft.set_ylim(energyLim)
	#setp(axLeft.get_yticklines(), "markersize", 2)

	axDensity.plot([0,singleIonBoxSize], [singleIonBoxSize,singleIonBoxSize], 'k-', linewidth=1)
	axDensity.plot([singleIonBoxSize,singleIonBoxSize], [0,singleIonBoxSize], 'k-', linewidth=1)
	axDensity.plot([singleIonBoxSize,singleIonBoxSize], [singleIonBoxSize,rmax], 'k--', linewidth=.6)
	axDensity.plot([singleIonBoxSize,rmax], [singleIonBoxSize,singleIonBoxSize], 'k--', linewidth=.6)
	axDensity.text(singleIonBoxSize/2, rmax-2, "SI", ha="center", va="center")
	axDensity.text(rmax-2, singleIonBoxSize/2, "SI", ha="center", va="center")
	axDensity.text(rmax-2, rmax-2, "DI", ha="center", va="center")

	#Set up info panel
	rectPanel = (rectBound[0], rectBound[1]+densityHeight+spacing, panelWidth, rectBound[3]-densityHeight-spacing)
	axPanel = fig.add_axes(rectPanel)
	axPanel.set_xticks([])
	axPanel.set_yticks([1])
	axPanel.bar([1], [1-(singleIon+doubleIon)], width=.7, facecolor=UiB_Green, linewidth=.3, align="center")
	axPanel.bar([2], [singleIon], width=.7, facecolor=UiB_Green, linewidth=.3, align="center")
	axPanel.bar([3], [doubleIon], width=.7, facecolor=UiB_Green, linewidth=.3, align="center")
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

	confSec = conf.LaserPotentialVelocityBase
	laserPulse = [confSec.time_function(confSec, t) for t in sampleTimes]
	axLaser.plot(sampleTimes, laserPulse, UiB_Green, linewidth=.8)
	axLaser.set_xlim([sampleTimes[0], sampleTimes[-1]])
	axLaser.set_ylim([-confSec.amplitude, confSec.amplitude])
	axLaser.plot([sampleTimes[timeIdx]]*2, [laserPulse[timeIdx]]*2, "o", color=UiB_Red, markeredgecolor=UiB_Red)
	draw()
	return fig


def CreateIonizationMovieFrameFromFile(filename, masterFile, radialGrid=None):

	if radialGrid == None:
		radialGrid = linspace(0, 30, 250)
	radialDensity = CreateRadialDensityFromFile(filename, radialGrid)

	CreateIonizationMovieFrame(radialGrid, radialDensity, singleIon, doubleIon, 280, sampleTimes, conf)


def CreateIonizationMovie(fileList, masterFile, offset=0):
	interactive = rcParams["interactive"]
	outputDir = "figs/stabilization_movie"
	radialGrid = linspace(0, 30, 250)

	PaperFigureSettings(8,10)

	#Get misc data from main (master) file
	sampleTimes, singleIon, doubleIon, conf = GetMasterData(masterFile)

	rcParams["interactive"] = False
	for idx, file in enumerate(fileList):
		curIdx = idx + offset
		print "Creating frame %03i/%03i..." % (curIdx, len(fileList)+offset)
		sys.stdout.flush()
		curFigname = "%s/helium_stabilization_freq_5_cycle_6_I_20_%03i.png" % (outputDir, curIdx)
		radialDensity = CreateRadialDensityFromFile(file, radialGrid)
		fig = CreateIonizationMovieFrame(radialGrid, radialDensity, singleIon[curIdx].real, doubleIon[curIdx].real, curIdx, sampleTimes, conf, vmax = 1.0)
		fig.savefig(curFigname, dpi=500)	
		close(fig)

	rcParams["interactive"] = interactive


def GetMasterData(masterFile):
	sampleTimes = GetArrayFromFile(masterFile, "SampleTimes")
	singleIon = GetArrayFromFile(masterFile, "SingleIonization")
	doubleIon = GetArrayFromFile(masterFile, "DoubleIonization")
	conf = pyprop.Config(pyprop.serialization.GetConfigFromHDF5(masterFile))

	return sampleTimes, singleIon, doubleIon, conf
