
#-------------------------------------------------------------------------------------------
#                                         Figure 2
#-------------------------------------------------------------------------------------------

def MakeFigure2Top():
	#tControl, cControl = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_all_phase_0.00pi_pump_5fs_5e13_scaling_2.h5", partitionCount=0)

	rcParams["text.usetex"] = False

	fig = figure(figsize=(8,3))

	cmap = LoadColormap("gradient.txt", True)
	color = "#0045a2" #"#ff8e16" 
	maxState = 10

	stateShift = 0.5
	c = cControl[::2,:maxState+1]
	t = tControl[::2]
	states = arange(maxState+1) - stateShift

	pcolormesh(t, states, c.transpose(), cmap=cmap, shading="flat", vmin=0.02, vmax=0.35)
	axis((0,650,-stateShift,maxState-stateShift))
	yticks([0,1,2,3,4,5,6,7,8,9])
	title("Vibrational Distribution")
	ylabel("Vibrational State")
	xlabel("Delay time (fs)")
	a = axis()
	cb = colorbar(ticks=[0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35])
	cb.ax.set_position([0.91, 0.16, 0.014, 0.73]);
	fig.subplots_adjust(left=0.08, bottom=0.16, right=0.89, top=0.89)
	text(30,1.5,"(a)", color="w")
	draw()

	fig.savefig("figure2/figure2-top.png", dpi=400)
	
	fig = figure(figsize=(8,3))
	axis(a)

	title("Vibrational Distribution")
	ylabel("Vibrational State")
	xlabel("Delay time (fs)")
	axis((0,650,-stateShift,maxState-stateShift))
	yticks([0,1,2,3,4,5,6,7,8,9])
	cb = colorbar(ticks=[0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35])
	cb.ax.set_position([0.91, 0.16, 0.014, 0.73]);
	fig.subplots_adjust(left=0.1, bottom=0.18, right=0.89, top=0.89)
	text(30,1.5,"(a)", color="w")
	
	draw()
	fig.savefig("figure2/figure2-top.svg")


def MakeFigure2Bottom():
	#tControl, cControl = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_all_phase_0.00pi_pump_5fs_5e13_scaling_2.h5", partitionCount=0)
	
	rcParams["text.usetex"] = False

	fig = figure(figsize=(8,2.5))
	cmap = LoadColormap("gradient.txt", True)
	color = "#0045a2" #"#ff8e16" 
	maxState = 12

	fcPop = (abs(GetFranckCondonPopulation(molecule="d2+", radialScaling=2))**2)[:maxState+1]

	t1 = where(tControl == 293)[0][0]
	c1 = cControl[t1, :maxState+1]
	t2 = where(tControl == 306)[0][0]
	c2 = cControl[t2, :maxState+1]

	ax = subplot(1,3,1)
	col = [cmap(v/0.35)[:3] for v in fcPop]
	CorrelationBarPlot(fcPop, ax, 0.3, colors=col)
	ax.text(0, 0.23, "(b)")
	ax.text(8.5, 0.23, "Initial")

	ax.set_ylabel("Probability")

	ax = subplot(1,3,2)
	col = [cmap(v/0.35)[:3] for v in c1]
	CorrelationBarPlot(c1, ax, 0.3, colors=col)
	ax.set_yticklabels([])
	ax.text(0, 0.23, "(c)")
	ax.text(8.5, 0.23, "293fs")


	ax.set_xlabel("Vibrational Eigenstate")

	ax = subplot(1,3,3)
	col = [cmap(v/0.35)[:3] for v in c2]
	CorrelationBarPlot(c2, ax, 0.3, colors=col)
	ax.set_yticklabels([])
	ax.text(0, 0.23, "(d)")
	ax.text(8.5, 0.23, "306fs")

	fig.subplots_adjust(left=0.1, bottom=0.20, right=0.89, top=0.89, wspace=0)
	draw()

	fig.savefig("figure2/figure2-bottom.svg")


def MakeFourierPlotCheckerboard():
	#tControl, cControl = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_all_phase_0.00pi_pump_5fs_5e13_scaling_2.h5", partitionCount=0)

	rcParams["text.usetex"] = True
	
	dT = diff(tControl)
	minDT = min(dT)
	maxDT = max(dT)
	if maxDT - minDT > 1e-10:
		raise Exception("Not equidistant time grid")
	dt = dT[0]

	#start = 500
	#end = 1800
	start = 0
	end = len(tControl)
	s = s_[start:end]
	N = len(tControl[s])
	k = fft.helper.fftfreq(N, dt*femtosec_to_au)[:N/2] * 2*pi
	W = hanning(N)

	fig = figure()
	ax = gca()

	for i in range(0,8):
		f = fft.fft(W * cControl[s,i])[:N/2]
		ax.plot(k, abs(f), label="$|%i\\rangle$" % i)

	lineColor = "#666666"
	E, V = LoadBoundEigenstates(molecule="d2+", radialScaling=2)
	for dE in diff(E): ax.axvline(dE, color=lineColor, linestyle=":", alpha="0.5")

	axis([0.004,0.008,0,80])
	legend()

	draw()


#-------------------------------------------------------------------------------------------
#                                         Figure 3
#-------------------------------------------------------------------------------------------

def MakeFourierPlotEvenOdd():
	rcParams["text.usetex"] = True

	#tControl, cControl = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_all_phase_0.00pi_pump_5fs_5e13_scaling_2.h5", partitionCount=0)
	#t1, c1 = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_control_all_control_293.00fs_5e13_probe_5fs_40e13_scaling_2.h5", partitionCount=0)
	#t2, c2 = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_control_all_control_306.00fs_5e13_probe_5fs_40e13_scaling_2.h5", partitionCount=0)

	s1 = s_[:]
	d1 = 1 - sum(c1[s1], axis=1)
	N1 = len(d1)
	k1 = fft.helper.fftfreq(N1, 1.0*femtosec_to_au)[:N1/2] * 2*pi

	s2 = s_[:]
	d2 = 1 - sum(c2[s2], axis=1)
	N2 = len(d2)
	k2 = fft.helper.fftfreq(N2, 1.0*femtosec_to_au)[:N2/2] * 2*pi

	E, V = LoadBoundEigenstates(molecule="d2+", radialScaling=2)

	#Convert atomic units of frequency (energy) to THz
	k1 = k1 / 2*pi * femtosec_to_au * 1000
	k2 = k2 / 2*pi * femtosec_to_au * 1000
	freq = E / 2*pi * femtosec_to_au * 1000
	
	fig = figure(figsize=(6,4))
	c = [0, 4000, 0, 1]
	b = [-1, 15, 0, 0.4]
	a = [500, 1000, 0, 100]
	shift = 2
	lineColor = "#666666"
	fig.subplots_adjust(wspace=0.05, hspace=0.05, left=0.06, bottom=0.10, right=0.95, top=0.92)

	ticks = diff(freq)
	labels = ["$v_{%i} - v_{%i}$" % (i+1, i) for i,e in enumerate(ticks)]

	#293fs
	ax = subplot(2, 1, 1)
	f1 = fft.fft(d1)[:N1/2]
	for dE in freq[2:] - freq[:-2]: ax.axvline(dE, color=lineColor, linestyle="--", alpha="0.5")
	for i, dE in enumerate((freq[2:] - freq[:-2])[:10]): ax.text(dE - shift, 90, "$\omega_{%i} - \omega_{%i}$" % (i+2, i), rotation="90", horizontalalignment="right", verticalalignment="top", fontsize=14)
	for dE in diff(freq): ax.axvline(dE, color=lineColor, linestyle=":", alpha="0.5")
	ax.plot(k1, abs(f1), color="k", label="293fs (states 2,3)")
	ax.set_yticks([])
	#xticks(ticks, labels, rotation=45, fontsize=12)
	ax.axis(a)
	ax.set_xticklabels([])
	ax.set_ylabel("Control @ 293fs")
	ax.set_title("Fourier Transformed Dissociation Yield")
	ax.text(0.0145, 80, "(a)")

	#306fs
	ax = subplot(2, 1, 2)
	f2 = fft.fft(d2)[:N1/2]
	for dE in freq[2:] - freq[:-2]: ax.axvline(dE, color=lineColor, linestyle="--", alpha="0.5")
	for dE in diff(freq): ax.axvline(dE, color=lineColor, linestyle=":", alpha="0.5")
	ax.plot(k2, abs(f2), color="k", label="306fs (states 3,4)")
	yticks([])
	ax.axis(a)
	#xticks(ticks, labels, rotation=45, fontsize=12)
	ax.set_ylabel("Control @ 306fs")
	ax.set_xlabel("Frequency (THz)")
	ax.text(0.0145, 80, "(b)")

	show()

	fig.savefig("figure3/figure3.eps")
	fig.savefig("figure3/figure3.pdf")


#-------------------------------------------------------------------------------------------
#                                         Figure 4
#-------------------------------------------------------------------------------------------

def DrawBigClock(theta, drawLines=True, drawTicks=True):
	ax = gca()

	winWidth = ax.get_window_extent().width()
	winHeight = ax.get_window_extent().height()

	if winHeight > winWidth:
		ax.axis([-1,1,-1*winHeight/winWidth, 1*winHeight/winWidth])
	else:
		ax.axis([-1*winWidth/winHeight, 1*winWidth/winHeight, -1,1])

	radius = 0.8
	space = 0.05
	
	ell = matplotlib.patches.Ellipse(xy=[0,0], width=2*radius, height=2*radius, facecolor="w")
	ax.add_artist(ell)

	if drawLines:
		lin1 = matplotlib.patches.lines.Line2D(xdata=[-radius,radius], ydata=[0,0], color="#aaaaaa")
		ax.add_artist(lin1)
		lin2 = matplotlib.patches.lines.Line2D(xdata=[0, 0], ydata=[-radius, radius], color="#aaaaaa")
		ax.add_artist(lin2)

	if drawTicks:
		ax.text(0, radius+space, "$0$", verticalalignment="bottom", horizontalalignment="center")
		ax.text(0, -(radius+space), "$\\pi$", verticalalignment="top", horizontalalignment="center")
		
		ax.text(radius+space, 0, "$\\pi/2$", verticalalignment="center", horizontalalignment="left")
		ax.text(-(radius+space), 0, "$3\\pi/2$", verticalalignment="center", horizontalalignment="right")

	dx = radius * cos(+pi/2 - theta)
	dy = radius * sin(+pi/2 - theta)
	lin3 = matplotlib.patches.lines.Line2D(xdata=[0,dx], ydata=[0,dy], color="b")
	ax.add_artist(lin3)
		
	ax.set_axis_off()

	draw_if_interactive()

	

def CorrelationBarPlotPhaseClockFig(corr, corr2, ax=None, axisHeight=None):
	if ax == None:
		ax = gca()
	
	basisCount = len(corr)
	left = r_[:basisCount] - 0.4

	#barHeight = abs(corr)**2)
	barHeight = abs(corr)**2
	ax.bar(left, barHeight)
	ax.bar(left[4], barHeight[4], color="red")

	if axisHeight == None:
		axisHeight = ax.axis()[3]
	axisWidth = len(corr) + 1

	winWidth = ax.get_window_extent().width()
	winHeight = ax.get_window_extent().height()

	scale = winWidth * axisHeight / (winHeight * axisWidth)

	clockHeight = scale / (1 - scale)
	yStart = - 2* clockHeight
	yEnd = axisHeight
	print yEnd
	
	xStart = -1
	xEnd = len(corr)

	ax.axis([xStart, xEnd, yStart, yEnd])

	#y-position of the clocks
	yCenter1 = -clockHeight / 2
	yCenter2 = -clockHeight - clockHeight / 2

	for i, curCorr in enumerate(corr[3:6]):
			i = i + 3
		#if abs(corr[i])**2 > 1e-2:
			theta1 = arctan2(corr[i].imag, corr[i].real)
			theta2 = arctan2(corr2[i].imag, corr2[i].real)
		
			xCenter1 = i
			xy = [xCenter1, yCenter1]
		
			width = 0.7
			height = 2 * abs(yCenter1) * 0.7
			ell = matplotlib.patches.Ellipse(xy=xy, width=width, height=height, facecolor="w")
			ax.add_artist(ell)
		
			dx = width/2. * cos(+pi/2 - theta1)
			dy = height/2. * sin(+pi/2 - theta1)
			lin = matplotlib.patches.lines.Line2D(xdata=[xCenter1, xCenter1+dx], ydata=[yCenter1, yCenter1+dy])
			ax.add_artist(lin)

			xCenter2 = i
			xy = [xCenter2, yCenter2]
			ell = matplotlib.patches.Ellipse(xy=xy, width=width, height=height, facecolor="w")
			ax.add_artist(ell)
		
			dx = width/2. * cos(+pi/2 - theta2)
			dy = height/2. * sin(+pi/2 - theta2)
			lin = matplotlib.patches.lines.Line2D(xdata=[xCenter2, xCenter2+dx], ydata=[yCenter2, yCenter2+dy])

			ax.add_artist(lin)

	ax.set_xticks(ax.get_xticks()[1:-1])
	yt = ax.get_yticks()

	#yt = yt[where(yt>0)]
	#yticks([yCenter1, yCenter2] + list(yt), ["293fs", "306fs"] + list(yt))
	#ax.axis([xStart, xEnd, yStart, yEnd])
	
	yt = yt[where(yt>=0)]
	yticks(yt)
	ax.axis([xStart, xEnd, yStart, yEnd])
	text(-0.35, yCenter1, "293fs", verticalalignment="center")
	text(-0.35, yCenter2, "306fs", verticalalignment="center")

	draw_if_interactive()
	
	return ax.axis()
		




def MakeClockFig():
	rcParams["text.usetex"] = True
	
	initStates = [3,4,5]
	pulseDelay1 = 293*femtosec_to_au
	pulseDelay2 = 306*femtosec_to_au
	args = {}
	args["radialScaling"] = 2
	args["pulseIntensity"] = 5e13
	args["pulseDuration"] = 5*sqrt(2)*femtosec_to_au
	args["molecule"] = "d2+"
	args["duration"] = 350*femtosec_to_au
	args["fastForward"] = 250*femtosec_to_au
	args["outputCount"] = 1
	args["relativePhase"] = "follow"
	args["silent"] = True
	args["configSilent"] = True
	args["initStatesWeight"] = "one"
	
	#args["dt"] = 30
	
	fig = figure(figsize=(6,5))
	fig.subplots_adjust(left=0.10, right=0.96, wspace=0.20, hspace=0.30)
	subplots = [2,3,4]
	labels = ["(a) $|3\\rangle$", "(b) $|4\\rangle$", "(c) $|5\\rangle$"]
	
	#fig.text(0.5,0.95,"Single State Propagations",verticalalignment="center",horizontalalignment="center", fontsize="medium")
		
	ax = subplot(2,2,1)
	DrawBigClock(theta=pi/4)
	draw()

	corr1 = []
	corr2 = []
	
	for i, state in enumerate(initStates):
		ax = subplot(2,2,subplots[i])
		args["initStates"] = [state]
		args["pulseDelay"] = pulseDelay1
		t, c, X, Y, Z, W = Propagate(**args)
		args["pulseDelay"] = pulseDelay2
		t, c2, X, Y, Z, W = Propagate(**args)
		CorrelationBarPlotPhaseClockFig(c[-1,:9], c2[-1,:9], ax=ax, axisHeight=1.)
		corr1.append(c[-1,:])
		corr2.append(c2[-1,:])
		#ax.set_xlabel("Vibrational Level")
		#ax.set_ylabel("Distribution Probability and Phase-Shift")
		#ax.set_title("Final distribution for initial v=%i. Pulse @ %ifs" % (state, pulseDelay / femtosec_to_au))
	
		text(8.3,0.8,labels[i], verticalalignment="center", horizontalalignment="right")

	subplot(2,2,2)
	xlabel("Vibrational State")
	ylabel("Probability")
	subplot(2,2,3)
	ylabel("Probability")
	xlabel("Vibrational State")
	subplot(2,2,4)
	xlabel("Vibrational State")

	savefig("clockfig/figure3.eps")
	savefig("clockfig/figure3.pdf")
	draw()


	return corr1, corr2



def CalculateClockFig2():
	initStates = [3,4,5]
	pulseDelay1 = 293*femtosec_to_au
	pulseDelay2 = 306*femtosec_to_au
	args = {}
	args["radialScaling"] = 2
	args["pulseIntensity"] = 5e13
	args["pulseDuration"] = 5*sqrt(2)*femtosec_to_au
	args["molecule"] = "d2+"
	args["duration"] = 350*femtosec_to_au
	args["fastForward"] = 250*femtosec_to_au
	args["outputCount"] = 1
	args["relativePhase"] = "follow"
	args["silent"] = True
	args["configSilent"] = True
	args["initStatesWeight"] = "one"
	
	#args["dt"] = 30

	delayRange = r_[280:310]
	corr = zeros((len(delayRange), 10), dtype=complex)
	
	for i, pulseDelay in enumerate(delayRange):
		print "pulseDelay = %ifs" % pulseDelay
		args["initStates"] = [5]
		args["pulseDelay"] = pulseDelay * femtosec_to_au
		t, c, X, Y, Z, W = Propagate(**args)
		corr[i,:] = c[-1,:10]
		
	return corr

def MakeClockFig2(corr=None):
	if corr == None:
		corr = CalculateClockFig2()

	fig = figure(figsize=(16,3))	
	fig.subplots_adjust(wspace=0, hspace=0, bottom=0, left=0, top=1, right=1)

	delayCount = corr.shape[0]
	stateCount = corr.shape[1]

	interactive = rcParams["interactive"]
	try:
		rcParams["interactive"] = False

		theta = arctan2(corr.imag, corr.real)

		for i, delay in enumerate(r_[280:310]):
			for j in range(stateCount):
				ax = subplot(stateCount, delayCount, (delayCount*j) + i + 1)
				DrawBigClock(theta=theta[i,j], drawLines=False, drawTicks=False)
				ax.set_axis_off()
			ax = subplot(stateCount, delayCount, i + 1)
			ax.text(0,1,"%ifs" % delay, verticalalignment="center", horizontalalignment="center")


	finally:
		rcParams["interactive"] = interactive

	draw()


def MakeClockFig3():

	E, V = LoadRotatedBoundEigenstates(molecule="d2+", radialScaling=2)
	E = diff(E)
	states = r_[0:10]
	#delayTimes = r_[280:295:1]
	delayTimes = [280, 289, 293, 306]

	fig = figure(figsize=(16,3))	
	subplot(len(delayTimes), len(states), 1)

	interactive = rcParams["interactive"]
	try:
		rcParams["interactive"] = False

		for j, curT in enumerate(delayTimes):
			for i, curE in enumerate(E[states]):
				theta = - curE * curT * femtosec_to_au
				print curE
				print (len(states)*j) + i + 1	
				ax = subplot(len(delayTimes), len(states), (len(states)*j) + i + 1)
				DrawBigClock(theta=theta, drawLines=False, drawTicks=False)
				ax.set_axis_off()


	finally:
		rcParams["interactive"] = interactive

	fig.subplots_adjust(wspace=0, hspace=0, bottom=0, left=0, top=1, right=1)
	draw()




#-------------------------------------------------------------------------------------------
#                          Coupling Figure
#-------------------------------------------------------------------------------------------


def MakeCouplingFig(**args):
	rcParams["text.usetex"] = True

	args["molecule"] = "d2+"
	args["radialScaling"] = 2

	matrix = real(GetTransitionMatrixElements(**args))
	figure(figsize=(6,4.5))
	
	cmap = LoadColormapMirrored("gradient.txt")
	imshow(matrix[:8,:8], interpolation="nearest", cmap=cmap, vmin=-2, vmax=2)
	colorbar()

	#title(r"Matrix elements $\langle i | c^2 | j \rangle$")
	xlabel("Vibrational State $| j \\rangle$")
	ylabel("Vibrational State $| i \\rangle$")

	savefig("figure4/figure4.pdf")
	savefig("figure4/figure4.eps")


