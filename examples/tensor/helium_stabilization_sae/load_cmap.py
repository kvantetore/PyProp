import matplotlib
import pylab
from numpy import linspace

def LoadColormap(filename, reverse=False):
	data = pylab.load(filename)

	samples = len(data)/4
	t = linspace(0,1,samples)
	r = list(data[0::4])
	g = list(data[1::4])
	b = list(data[2::4])

	if reverse:
		r.reverse()
		g.reverse()
		b.reverse()

	red = []
	green = []
	blue = []

	for i in range(samples):

		red.append((t[i], r[i], r[i]))
		green.append((t[i], g[i], g[i]))
		blue.append((t[i], b[i], b[i]))

	cdict = { "red": red, "green": green, "blue": blue }
	cmap = matplotlib.colors.LinearSegmentedColormap("my_colors", cdict, 1024)

	return cmap

def LoadColormapMirrored(filename):
	data = pylab.load(filename)

	samples = len(data)/2
	t = linspace(0,1,samples)
	r = list(data[0::4])
	g = list(data[1::4])
	b = list(data[2::4])

	r.reverse()
	g.reverse()
	b.reverse()

	r = list(reversed(b)) + r 
	g = list(reversed(g)) + g  
	b = list(reversed(r)) + b 

	red = []
	green = []
	blue = []

	for i in range(samples):

		red.append((t[i], r[i], r[i]))
		green.append((t[i], g[i], g[i]))
		blue.append((t[i], b[i], b[i]))

	cdict = { "red": red, "green": green, "blue": blue }
	cmap = matplotlib.colors.LinearSegmentedColormap("my_colors", cdict, 1024)

	return cmap


