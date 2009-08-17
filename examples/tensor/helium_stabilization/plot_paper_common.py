import matplotlib
#matplotlib.use("qtagg")
matplotlib.use("agg")
from matplotlib.transforms import Bbox
import os

execfile("load_cmap.py")

UiB_Orange = "#d95900"
UiB_Green  = "#77af00"
UiB_Red    = "#aa0000"
UiB_Blue   = "#005473"
UiB_Black  = "#000000"

FigWidth = None
FigWidthLarge = None
LineWidth = None
GraphLineWidth = None
FontSize = None
PaperMaxE0 = None
PaperFigureFont = None 
FontFamily = "Times New Roman"
GradientFile = "gradient_uib.txt"

PlotOutputPaper = "paper"
PlotOutputIcpeac = "icpeac"

PaperOutput = None
def PaperSetOutput(plotOutput, polar=False):
	global FigWidth, FigWidthLarge, LineWidth, GraphLineWidth
	global FontSize, PaperOutput, PaperMaxE0, PaperFigureFont

	if plotOutput == PlotOutputPaper:
		FigWidth = 5
		FigWidthLarge = 8
		LineWidth = 0.4
		GraphLineWidth = 0.5
		FontSize = 6
		if polar:
			FontSize = 6
		PaperMaxE0 = 35
		PaperFigureFont = matplotlib.font_manager.FontProperties(family=FontFamily, size=FontSize)


	elif plotOutput == PlotOutputIcpeac:
		FigWidth = 12
		FigWidthLarge = 16
		LineWidth = 1
		GraphLineWidth = 1
		FontSize = 16
		PaperMaxE0 = 25
		PaperFigureFont = matplotlib.font_manager.FontProperties(family=FontFamily, size=FontSize)

	else:
		raise Exception("Unknown output")

	PaperOutput = plotOutput


def WalkChildren(parent):
	if hasattr(parent, "get_children"):
		for child in parent.get_children():
			yield child
			for subChild in WalkChildren(child):
				yield subChild


def HideUnusedTickLabels(parent):
	for artist in WalkChildren(parent):
		if isinstance(artist, matplotlib.axis.Tick):
			if not artist.label1On:
				artist.label1.set_visible(False)
			if not artist.label2On:
				artist.label2.set_visible(False)

def GetOptimalAxesPosition(fig, ax, additionalArtists=[], padding=5):
	rendr = fig.canvas.get_renderer()
	#import cairo
	##w_inch, h_inch = fig.get_size_inches()
	##width, height = 72*w_inch, 72*h_inch
	#print fig.dpi
	#width, height = fig.canvas.get_width_height()
	#rendr = matplotlib.backends.backend_cairo.RendererCairo(fig.dpi)
	#rendr.set_width_height(width, height)
	#surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
	#rendr.set_ctx_from_surface(surface)

	axBounds = ax.bbox
	getArea = lambda bbox: bbox.width * bbox.height
	def getExtent(artist):
		try:
			return artist.get_window_extent(rendr)
		except:
			return artist.get_window_extent()

	hasExtent = lambda artist: hasattr(artist, "get_window_extent") and getArea(getExtent(artist)) > 0
	isVisible = lambda artist: not hasattr(artist, "get_visible") or artist.get_visible()

	HideUnusedTickLabels(fig)
	axBoxes = map(getExtent, filter(isVisible, filter(hasExtent, WalkChildren(ax))))
	addBoxes = [map(getExtent, filter(isVisible,filter(hasExtent, WalkChildren(artist)))) for artist in additionalArtists] 
	totalBounds = Bbox.union(axBoxes + map(Bbox.union, addBoxes))

	width = rendr.width
	height = rendr.height
	xmin = ((axBounds.xmin - totalBounds.xmin) + padding) 
	ymin = ((axBounds.ymin - totalBounds.ymin) + padding) 
	xmax = (width  + (axBounds.xmax - totalBounds.xmax) - padding)  
	ymax = (height + (axBounds.ymax - totalBounds.ymax) - padding) 
	t = matplotlib.transforms.BboxTransformFrom(fig.bbox)
	newPos = Bbox(t.transform(Bbox([[xmin, ymin], [xmax, ymax]])))

	return newPos			



#------------------ Setup Matplotlib -------------------------

def PaperGetFont():
	return PaperFigureFont

def PaperSetAllFonts(parent):
	if hasattr(parent, "get_children"):
		for chld in parent.get_children():
			if hasattr(chld, "set_fontproperties"):
				chld.set_fontproperties(PaperGetFont())
			PaperSetAllFonts(chld)


def PaperFigureSettings(fig_width, fig_height):
	cm_to_inch = 0.393700787
	fig_width *= cm_to_inch
	fig_height *= cm_to_inch
	fig_size =  [fig_width,fig_height]
	params = {
			  'toolbar': "none",
			  'axes.linewidth': LineWidth,
			  'lines.linewidth': LineWidth,
			  'text.usetex': False,
			  'figure.figsize': fig_size,
			  'font.family' : 'serif',
			  'font.serif' : 'Arial'}
	matplotlib.rcParams.update(params)

def PaperUpdateFigure(fig):
	PaperSetAllFonts(fig)
	fig.figurePatch.set_visible(False)
	for ax in fig.axes:
		ax.axesPatch.set_visible(False)
		for ln in ax.get_xticklines() + ax.get_yticklines():
			ln.set_linewidth(LineWidth)
	draw()

def PaperUpdatePolarFigure(fig):
	#PolarAxColor = "#c5c5c5"
	PolarAxColor = "#bababa"
	#PolarAxColor = "#767676"
	#PolarAxColor = "gray"
	PolarAxLineWidth = 0.4
	
	PaperSetAllFonts(fig)

	fig.figurePatch.set_visible(False)
	for ax in fig.axes:
		ax.axesPatch.set_visible(False)
		for ln in ax.get_xgridlines() + ax.get_ygridlines():
			ln.set_color(PolarAxColor)
			ln.set_linewidth(PolarAxLineWidth)
			ln.set_linestyle("-")
	draw()

def PaperGetFolder():
	if not PaperOutput:
		raise Exception("Please call PaperSetOutput before plotting")
	folder = "figs/%s" % PaperOutput
	if not os.path.exists(folder):
		os.makedirs(folder)
	return folder



