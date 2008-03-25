import matplotlib
import pylab

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


