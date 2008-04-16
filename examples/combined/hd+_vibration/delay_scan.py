import commands
import pyprop.utilities.submitpbs as submit
import numpy

cm_to_inch = 1./2.5
figureSizePaper = 16*cm_to_inch
figureSizeScreen = 32*cm_to_inch

try: 
	import scipy.interpolate
	ScipyAvailable=True
except: 
	print "Scipy not available"
	ScipyAvailable=False


def PrintAxisPosition(event):
	if event.inaxes:
		print "Clicked at (%.4f, %.4f)"  % (event.xdata, event.ydata)
		figure()

#-----------------------------------------------------------------
#            Parse input files
#-----------------------------------------------------------------

def MakeFourierPlotEvenOdd():
	#tControl, cControl = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_all_phase_0.00pi_pump_5fs_5e13_scaling_2.h5", partitionCount=0)
	#t1, c1 = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_control_all_control_293.00fs_5e13_probe_5fs_30e13_scaling_2.h5", partitionCount=0)
	#t2, c2 = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_control_all_control_306.00fs_5e13_probe_5fs_30e13_scaling_2.h5", partitionCount=0)

	s1 = s_[:]
	d1 = 1 - sum(c1[s1], axis=1)
	N1 = len(d1)
	k1 = fft.helper.fftfreq(N1, 1.0*femtosec_to_au) * 2*pi

	s2 = s_[:]
	d2 = 1 - sum(c2[s2], axis=1)
	N2 = len(d2)
	k2 = fft.helper.fftfreq(N2, 1.0*femtosec_to_au) * 2*pi

	E, V = LoadBoundEigenstates(molecule="d2+", radialScaling=2)
	
	figure(figsize=(11,11))
	b = [-1, 15, 0, 0.4]
	a = [0.005, 0.016, 0, 230]
	subplots_adjust(wspace=0.05)

	ticks = diff(E)
	labels = ["$v_{%i} - v_{%i}$" % (i+1, i) for i,e in enumerate(ticks)]

	#293fs
	ax = subplot(2, 2, 0 + 2*0 + 1)
	tIndex = where(tControl == 293)[0][0]
	CorrelationBarPlot(cControl[tIndex,:], ax)
	ax.axis(b)
	xticks(fontsize=12)
	ylabel("Control Pulse @ 293fs", fontsize=12)
	title("Vibrational Distribution", fontsize=12)

	ax = subplot(2, 2, 1 + 2*0 + 1)
	ax.plot(k1[:N1/2], abs(fft.fft(d1)[:N1/2]), label="293fs (states 2,3)")
	for dE in E[2:] - E[:-2]: ax.axvline(dE, color="r")
	for dE in diff(E): ax.axvline(dE, color="k")
	yticks([])
	#xticks(ticks, labels, rotation=45, fontsize=12)
	ax.axis(a)
	title("Fourier Transformed Dissociation Yield", fontsize=12)

	#306fs
	ax = subplot(2, 2, 0 + 2*1 + 1)
	tIndex = where(tControl == 306)[0][0]
	CorrelationBarPlot(cControl[tIndex,:], ax)
	ax.axis(b)
	xticks(fontsize=12)
	ylabel("Control Pulse @ 306fs", fontsize=12)
	
	ax = subplot(2, 2, 1 + 2*1 + 1)
	ax.plot(k2[:N2/2], abs(fft.fft(d2)[:N2/2]), label="306fs (states 3,4)")
	for dE in E[2:] - E[:-2]: ax.axvline(dE, color="r")
	for dE in diff(E): ax.axvline(dE, color="k")
	yticks([])
	#xticks(ticks, labels, rotation=45, fontsize=12)
	ax.axis(a)

	
	#Second figure
	figure()
	plot(t1, 1-sum(c1, axis=1), "r", label="22fs (states 2,3)")
	#plot(t1[s1], d1, "r,")

	plot(t2, 1-sum(c2, axis=1), "g", label="24.5fs (states 3,4)")
	#plot(t2[s2], d2, "g,")

	legend()
	show()


def MakeFourierPlot():
	#tControl, cControl = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_all_phase_0.00pi_pump_5fs_20e13_scaling_2.h5", partitionCount=0)
	#t1, c1 = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_control_all_control_22.00fs_20e13_probe_5fs_30e13_scaling_2.h5", partitionCount=0)
	#t2, c2 = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_control_all_control_24.50fs_20e13_probe_5fs_30e13_scaling_2.h5", partitionCount=0)
	#t3, c3 = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_control_all_control_26.20fs_20e13_probe_5fs_30e13_scaling_2.h5", partitionCount=0)

	s1 = s_[:]
	d1 = 1 - sum(c1[s1], axis=1)
	N1 = len(d1)
	k1 = fft.helper.fftfreq(N1, 1.0*femtosec_to_au) * 2*pi

	s2 = s_[:]
	d2 = 1 - sum(c2[s2], axis=1)
	N2 = len(d2)
	k2 = fft.helper.fftfreq(N2, 1.0*femtosec_to_au) * 2*pi

	s3 = s_[:]
	d3 = 1 - sum(c3[s3], axis=1)
	N3 = len(d3)
	k3 = fft.helper.fftfreq(N3, 1.0*femtosec_to_au) * 2*pi


	E, V = LoadBoundEigenstates(molecule="d2+", radialScaling=2)
	
	figure(figsize=(11,11))
	b = [-1, 15, 0, 0.4]
	a = [0.00480, 0.0080, 0, 230]
	subplots_adjust(wspace=0.05)

	ticks = diff(E)
	labels = ["$v_{%i} - v_{%i}$" % (i+1, i) for i,e in enumerate(ticks)]

	#22fs
	ax = subplot(3, 2, 0 + 2*0 + 1)
	tIndex = where(tControl == 22)[0][0]
	CorrelationBarPlot(cControl[tIndex,:], ax)
	ax.axis(b)
	xticks(fontsize=12)
	ylabel("Control Pulse @ 22fs", fontsize=12)
	title("Vibrational Distribution", fontsize=12)

	ax = subplot(3, 2, 1 + 2*0 + 1)
	ax.plot(k1[:N1/2], abs(fft.fft(d1)[:N1/2]), label="22fs (states 2,3)")
	for dE in E[2:] - E[:-2]: ax.axvline(dE, color="r")
	#for dE in diff(E): ax.axvline(dE, color="k")
	yticks([])
	xticks(ticks, labels, rotation=45, fontsize=12)
	ax.axis(a)
	title("Fourier Transformed Dissociation Yield", fontsize=12)

	#24.5fs
	ax = subplot(3, 2, 0 + 2*1 + 1)
	tIndex = where(tControl == 24.5)[0][0]
	CorrelationBarPlot(cControl[tIndex,:], ax)
	ax.axis(b)
	xticks(fontsize=12)
	ylabel("Control Pulse @ 24.5fs", fontsize=12)
	
	ax = subplot(3, 2, 1 + 2*1 + 1)
	ax.plot(k2[:N2/2], abs(fft.fft(d2)[:N2/2]), label="24.5fs (states 3,4)")
	for dE in E[2:] - E[:-2]: ax.axvline(dE, color="r")
	#for dE in diff(E): ax.axvline(dE, color="k")
	yticks([])
	xticks(ticks, labels, rotation=45, fontsize=12)
	ax.axis(a)

	#26.2fs
	ax = subplot(3, 2, 0 + 2*2 + 1)
	tIndex = where(tControl == 26.5)[0][0]
	CorrelationBarPlot(cControl[tIndex,:], ax)
	ax.axis(b)
	xticks(fontsize=12)
	ylabel("Control Pulse @ 26.2fs", fontsize=12)
		
	ax = subplot(3, 2, 1 + 2*2 + 1)
	ax.plot(k3[:N3/2], abs(fft.fft(d3)[:N3/2]), label="26.2fs (states 4,5)")
	for dE in E[2:] - E[:-2]: ax.axvline(dE, color="r")
	#for dE in diff(E): ax.axvline(dE, color="k")
	yticks(E)
	xticks(ticks, labels, rotation=45, fontsize=12)
	ax.axis(a)

	figure()
	plot(t1, 1-sum(c1, axis=1), "r", label="22fs (states 2,3)")
	#plot(t1[s1], d1, "r,")

	plot(t2, 1-sum(c2, axis=1), "g", label="24.5fs (states 3,4)")
	#plot(t2[s2], d2, "g,")

	plot(t3, 1-sum(c3, axis=1), "b", label="26.2fs (states 4,15)")
	#plot(t3[s3], d3, "b,")
	legend()
	show()

class DelayScanPlot():
	def __init__(self, molecule, filename, partitionCount=0, figureSize=30*cm_to_inch, batchMode=False, radialScaling=1, useExistingFigure=False, figureTitle=None, cmap=None, color=None):
		#Disable interactive when plotting
		interactive = rcParams["interactive"]
		rcParams["interactive"] = False

		try:
		
			path, name = os.path.split(filename)
			root, ext = os.path.splitext(name)
			
			#Default color and colormap
			if cmap == None:
				cmap = LoadColormap("gradient.txt", True)
			if color == None:
				color = "#0045a2" #"#ff8e16" 
			self.ColorMap = cmap
			self.Color = color
		
			#Create New/Use Existing figure
			if useExistingFigure:
				curFig = gcf()
			else:
				curFig = figure(figsize=(25*cm_to_inch, 15*cm_to_inch))
			subplot(3,1,1)
			if figureTitle == None:
				figureTitle = "%s %s" % (molecule, root)
			title(figureTitle)
			self.Figure = curFig
		
			if not batchMode:
				pylab.connect("button_press_event", self.MouseClick)
			
			#load data
			t, eigenstateDistribution = GetScanDelayCorrelation(outputfile=filename, partitionCount=partitionCount)
			t, E, energyDistribution = GetScanDelayEnergyDistribution(molecule=molecule, outputfile=filename, partitionCount=partitionCount, radialScaling=radialScaling, recalculate=True)
			self.Time = t
			self.Energy = E
			self.EnergyDistribution = energyDistribution
			self.EigenstateDistribution = eigenstateDistribution
		
			#Plot
			self.PlotEigenstateDistribution()
			self.PlotDissociation()
			self.PlotEnergyDistribution()
			self.AdjustFigure()

		finally:
			#restore interactive state
			rcParams["interactive"] = interactive
			if not batchMode:
				show()

	
	def PlotEigenstateDistribution(self):
		t = self.Time

		#Plot Final Correlation
		subplot(3,1,3)
		maxState = 10
		pcolormesh(t, arange(maxState+1), self.EigenstateDistribution[:,:maxState+1].transpose(), cmap=self.ColorMap, shading="flat")
		axis((0,800,0,maxState))
		ylabel("Vibrational State")
		xlabel("Delay time (fs)")
		self.AxesEigenstateDistrib = gca()
	
	def PlotDissociation(self):
		t = self.Time

		#Plot Dissociation
		boundStateThreshold = 25
		dissociation = 1 - sum(self.EigenstateDistribution[:,:boundStateThreshold], axis=1)
		subplot(3,1,1)
		cla()
		plot(t, dissociation, color=self.Color)
		axis([0,800,0,1])
		ylabel("Dissoc. Probability")
		self.AxesDissoc = gca()
	
	def PlotEnergyDistribution(self):
		E = self.Energy
		t = self.Time

		#Plot Energy Distribution
		subplot(3,1,2)
		Emax = 2
		Emin = 0.0
		eIndex = where(E >Emin)[0]
		eIndex = eIndex[where(E[eIndex] <= Emax)[0]][::1]
		pcolormesh(t, E[eIndex], sum(self.EnergyDistribution[:,:,eIndex], axis=1).transpose(), shading="flat", cmap=self.ColorMap, vmin=0, vmax=20)
		axis((0,800,0,Emax))
		ylabel("Projectile Energy")
		self.AxesEnergyDistrib = gca()

	def AdjustFigure(self):
		#Fixup
		subplots_adjust(hspace=0)
		subplot(3,1,1)
		ticPos, ticLabel = xticks()
		ticLabel = ["" for i in range(len(ticLabel))]
		xticks(ticPos, ticLabel)
		
		subplot(3,1,2)
		ticPos, ticLabel = xticks()
		ticLabel = ["" for i in range(len(ticLabel))]
		xticks(ticPos, ticLabel)
		ticPos, ticLabel = yticks()
		yticks(ticPos[:-1])
		
		subplot(3,1,3)
		ticPos, ticLabel = yticks()
		yticks(ticPos[:-1]+0.5, ["   %i" % i for i in ticPos[:-1]])

	def MouseClick(self, event):
		interactive = rcParams["interactive"]
		rcParams["interactive"] = False

		if event.inaxes == self.AxesDissoc:
			self.PlotDissociation()

		if event.inaxes == self.AxesEnergyDistrib:
			clickedTime = event.xdata
			timeDiff = abs(self.Time - clickedTime)
			idx = where(timeDiff == numpy.min(timeDiff))[0][0]
			
			subplot(3,1,1)
			cla()
			hold(True)
			plot(self.Energy, self.EnergyDistribution[idx, 0, :], "r-", label="Energy distrib. t=%.2ffs binding" % self.Time[idx])
			plot(self.Energy, self.EnergyDistribution[idx, 1, :], "g--", label="Energy distrib. t=%.2ffs unbinding" % self.Time[idx])
			legend()
			axis((0,2,0,30))

		if event.inaxes == self.AxesEigenstateDistrib:
			clickedTime = event.xdata
			timeDiff = abs(self.Time - clickedTime)
			idx = where(timeDiff == numpy.min(timeDiff))[0][0]
	
			subplot(3,1,1)
			cla()
			hold(True)
			c = self.EigenstateDistribution
			bar(r_[0:c.shape[1]].copy(), c[idx, :].copy())
			title("t = %f" % self.Time[idx])

		rcParams["interactive"] = True	
		draw_if_interactive()

	

def MakeDelayScanPlot(molecule, filename, partitionCount=0, figureSize=30*cm_to_inch, batchMode=False, radialScaling=1, useExistingFigure=False, figureTitle=None, cmap=None, color=None):
	#Disable interactive when plotting
	interactive = rcParams["interactive"]
	rcParams["interactive"] = False

	path, name = os.path.split(filename)
	root, ext = os.path.splitext(name)

	#Default color and colormap
	if cmap == None:
		cmap = LoadColormap("gradient.txt", True)
	if color == None:
		color = "#0045a2" #"#ff8e16" 

	if not useExistingFigure:
		figure(figsize=(25*cm_to_inch, 15*cm_to_inch))
	subplot(3,1,1)
	if figureTitle == None:
		figureTitle = "%s %s" % (molecule, root)
	title(figureTitle)

	#Get norm
	#t, n = GetScanDelayNorm(outputfile=filename, partitionCount=partitionCount)

	#Plot Final Correlation
	t, c = GetScanDelayCorrelation(outputfile=filename, partitionCount=partitionCount)
	##pcolor plot
	subplot(3,1,3)
	maxState = 10
	pcolormesh(t, arange(maxState+1), c[:,:maxState+1].transpose(), cmap=cmap, shading="flat")
	axis((0,800,0,maxState))
	ylabel("Vibrational State")
	xlabel("Delay time (fs)")

	#Plot Ionization
	boundStateThreshold = 25
	ionization = 1 - sum(c[:,:boundStateThreshold], axis=1)
	#figure(figsize=(figureSize,figureSize*3./4.))
	subplot(3,1,1)
	plot(t, ionization, color=color)
	axis([0,800,0,1])
	ylabel("Dissoc. Probability")
	#xlabel("Pulse Delay")
	#title("Disociation yield %s %s" % (molecule, root))

	#Plot Energy Distribution
	t, E, c = GetScanDelayEnergyDistribution(molecule=molecule, outputfile=filename, partitionCount=partitionCount, radialScaling=radialScaling, recalculate=True)
	#figure(figsize=(figureSize, figureSize*3./4))
	subplot(3,1,2)
	Emax = 2
	Emin = 0.0
	eIndex = where(E >Emin)[0]
	eIndex = eIndex[where(E[eIndex] <= Emax)[0]][::1]
	pcolormesh(t, E[eIndex], sum(c[:,:,eIndex], axis=1).transpose(), shading="flat", cmap=cmap, vmin=0, vmax=20)
	#xlabel("Pulse Delay")
	axis((0,800,0,Emax))
	ylabel("Projectile Energy")
	#title("Energy distribution %s %s" % (molecule, root))

	#Fixup
	subplots_adjust(hspace=0)
	subplot(3,1,1)
	ticPos, ticLabel = xticks()
	ticLabel = ["" for i in range(len(ticLabel))]
	xticks(ticPos, ticLabel)

	subplot(3,1,2)
	ticPos, ticLabel = xticks()
	ticLabel = ["" for i in range(len(ticLabel))]
	xticks(ticPos, ticLabel)
	ticPos, ticLabel = yticks()
	yticks(ticPos[:-1])

	subplot(3,1,3)
	ticPos, ticLabel = yticks()
	yticks(ticPos[:-1]+0.5, ["   %i" % i for i in ticPos[:-1]])

	#restore interactive state
	rcParams["interactive"] = interactive
	if not batchMode:
		show()

	

def MakeDelayScanPlots(save=False):
	interactive = rcParams["interactive"]
	rcParams["interactive"] = False
	size = 16 * cm_to_inch

	def Plot(filename):
		t, c = GetScanDelayCorrelation(outputfile=filename, partitionCount=8)
		figure(figsize=(size,size*3./4.))
		ionization = 1 - sum(c[:,:20], axis=1)
		plot(t, ionization, label="Ionization")
		for i in range(10):
			if i < 5:
				label = "$n = %i$" % i
			else:
				label = "_nolegend_"
			plot(t, c[:,i], label=label)
		axis([20, 100, 0,  0.8])
		legend()
		xlabel("Pulse Delay")
		ylabel("Proabability")

	def PlotEnergy(molecule, filename):
		t, E, c = GetScanDelayEnergyDistribution(molecule=molecule, outputfile=filename, partitionCount=8)
		figure(figsize=(size, size*3./4))
		pcolormesh(t, E, sum(c, axis=1).transpose(), shading="flat")
		xlabel("Pulse Delay")
		ylabel("Energy")

		
	PlotEnergy("hd+", "outputfiles/hd+/delay_%i.h5")
	title("$HD^+$")
	Plot("outputfiles/hd+/delay_%i.h5")
	title("$HD^+$")
	if save:
		pylab.savefig("doc/fig/delay_scan_hdp.eps")
	
	PlotEnergy("h2+", "outputfiles/h2+/delay_%i.h5")
	title("$H_2^+$")
	Plot("outputfiles/h2+/delay_%i.h5")
	title("$H_2^+$")
	if save:
		pylab.savefig("doc/fig/delay_scan_h2p.eps")
	
	PlotEnergy("d2+", "outputfiles/d2+/delay_%i.h5")
	title("$D_2^+$")
	Plot("outputfiles/d2+/delay_%i.h5")
	title("$D_2^+$")
	if save:
		pylab.savefig("doc/fig/delay_scan_d2p.eps")
	
	#Plot("outputfiles/hd+/nodipole_delay_%i.h5")
	#title("$HD^+$ without static dipole moment")
	#if save:
	#	pylab.savefig("doc/fig/delay_scan_hdp_nodipole.eps")
	#
	##Plot the difference between HD+ with and without the static dipole term
	#figure(figsize=(size,size*3./4.))
	#t, c = GetScanDelayCorrelation(outputfile="outputfiles/hd+/delay_%i.h5", partitionCount=1)
	#t2, c2 = GetScanDelayCorrelation(outputfile="outputfiles/hd+/nodipole_delay_%i.h5", partitionCount=1)
	#plot(t, c - c2)
	#ax = axis()
	#ax[0] = 20
	#ax[1] = 100
	#axis(ax)
	#title("Difference between HD+ with and without\nstatic dipole moment")
	#xlabel("Pulse Delay")
	#ylabel("Proabability difference")
	#if save:
	#	pylab.savefig("doc/fig/delay_scan_nodipole_diff.eps")


	#restore previous interactive setting
	rcParams["interactive"] = interactive
	show()


#------------------------------------------------------------------------------	
#                Parse output files
#------------------------------------------------------------------------------	

def IterateScanNodes(outputfile, partitionCount=0):
	"""
	Iterate through all HDF files
		[outputfile % i for i in range(partitionCount)]

	For each file list through all root-level nodes that
	starts with "delay_", and yield that node along with the
	number behind the "delay_"

	This function is used to iterate through all outputfiles
	generated by SubmitScanDelay
	"""

	if partitionCount == 0:
		filenames = [outputfile]
	else:
		filenames = [outputfile % i for i in range(partitionCount)]

	for filename in filenames:
		f = tables.openFile(filename, "r")
		try:
			for node in f.listNodes(f.root):
				name = node._v_name
				if name.startswith("delay_"):
					try:
						t = float(name[len("delay_"):])
						yield t, node
					except:
						print "Could not process %s" % name
		finally:
			f.close()



def GetScanDelayEnergyDistribution(**args):
	"""
	Uses IterateScanNodes to iterate through all outputfiles, and calculates
	the energy distribution of the last wavefunction at each delay time

	It returns three variables
	1) delay times
	2) energy values
	3) probability[delay, energy]

	"""
	molecule = args["molecule"]
	outputfile = args["outputfile"]
	partitionCount = args["partitionCount"]

	E1, V1, E2, V2 = LoadContinuumEigenstates(**args)

	recalculate = False
	if "recalculate" in args:
		recalculate = args["recalculate"]

	#To add up the energies we must interpolate the values from 
	#the E1 and E2 grids to a common E-grid
	dE = average(diff(E2))
	E = r_[-0.5:0:dE]

	projList = []
	timeList = []

	for t, node in IterateScanNodes(outputfile, partitionCount):
		if recalculate:
			psi = node.wavefunction[-1,:,:]
			proj1, proj2 = CalculateEnergyDistribution(psi, E, E1, V1, E2, V2) 
		else:
			proj1 = node.energyDistribution[-1,0,:]
			proj2 = node.energyDistribution[-1,1,:]

		projList.append([proj1, proj2])
		timeList.append(t)

	sortedIndex = argsort(timeList)
	time = array(timeList)[sortedIndex]
	proj= array(projList)[sortedIndex,:] 

	#If the probability density is less than 0, it is an interpolation error,
	#To make the plots more consistent, we clamp it to 0
	proj[where(proj<0.0)] = 0.0

	E = (E + 0.5) / 2 / eV_to_au

	return time, E, proj




def CalculateEnergyDistribution(psi, outputEnergies, E1, V1, E2, V2):
	"""
	Calculates the energy distribution of the array psi, by projecting 
	on to the basis function sets V1 and V2.
	V1 is the eigenstates for the binding potential
	V2 is the eigenstates for the unbinding potential

	the energy distribution is then interpolated (with cubic splines)
	over outputEnergies
	"""
	proj1 = abs(dot(V1[:-1,:], psi[:,0]))**2 / diff(E1)
	proj2 = abs(dot(V2[:-1,:], psi[:,1]))**2 / diff(E2)

	if ScipyAvailable:
		tick1 = scipy.interpolate.interp1d(E1[:-1], proj1, bounds_error=False, fill_value=0)
		tick2 = scipy.interpolate.interp1d(E2[:-1], proj2, bounds_error=False, fill_value=0)
		#outDistrib1 = scipy.interpolate.splev(outputEnergies, tick1)
		#outDistrib2 = scipy.interpolate.splev(outputEnergies, tick2)
		outDistrib1 = tick1(outputEnergies)
		outDistrib2 = tick2(outputEnergies)
	
	else:
		interp1 = spline.Interpolator(E1[:-1], proj1)
		interp2 = spline.Interpolator(E2[:-1], proj2)
		outDistrib1 = array([interp1.Evaluate(e) for e in outputEnergies])
		outDistrib2 = array([interp2.Evaluate(e) for e in outputEnergies])

	return outDistrib1, outDistrib2



def RepackDelayScan(**args):
	outputfile = args["outputfile"]
	partitionCount = args["partitionCount"]
	repackFile = args["repackFile"]
	
	print "Repacking correlation"
	time, proj = GetScanDelayCorrelation(**args)
	print "Repacking energy distribution"
	t, E, corr = GetScanDelayEnergyDistribution(**args)
	print "Repacking norm"
	t, norm = GetScanDelayNorm(**args)

	
	output = tables.openFile(repackFile, "w")
	try:
		SaveArray(output, "/", "pulse_delay", time)
		SaveArray(output, "/", "energy", E)

		SaveArray(output, "/", "final_correlation", proj)
		SaveArray(output, "/", "energy_distribution", corr)
		SaveArray(output, "/", "norm", norm)

	finally:
		output.close()

	
def GetScanDelayCorrelation(**args):
	outputfile = args["outputfile"]
	partitionCount = args["partitionCount"]

	projList = []
	timeList = []

	for t, node in IterateScanNodes(outputfile, partitionCount):
		proj = node.eigenstateProjection[-1,:]
		projList.append(proj)
		timeList.append(t)

	sortedIndex = argsort(timeList)
	time = array(timeList)[sortedIndex]
	proj= array(projList)[sortedIndex]

	return time, proj

def GetScanDelayNorm(**args):
	outputfile = args["outputfile"]
	partitionCount = args["partitionCount"]

	projList = []
	timeList = []

	for t, node in IterateScanNodes(outputfile, partitionCount):
		proj = node.norm[-1]
		projList.append(proj)
		timeList.append(t)

	sortedIndex = argsort(timeList)
	time = array(timeList)[sortedIndex]
	proj= array(projList)[sortedIndex]

	return time, proj


#------------------------------------------------------------------------------	
#                Submit scan delay jobs
#------------------------------------------------------------------------------	

#ipython1:
DefaultControllerHost = "localhost"
DefaultControllerPort = 61001

try:
	import ipython1.kernel.api as kernel
except:
	print "Could not load IPython1, fancy submitting wil be unavailable"

def GetStalloEngineCount():
	controllerHost = DefaultControllerHost
	controllerPort = DefaultControllerPort

	#Create connection to stallo
	rc = kernel.RemoteController((controllerHost, controllerPort))
	return rc.getIDs()

def SubmitAll():
	args = {}
	args["delayList"] = r_[0:800:2]
	args["radialScaling"] = 2

	molecules = ["d2+", "hd+", "h2+"]
	pulseDurations = [8, 12, 18]
	controlIntensities = [0, 0.5e14, 1e14, 2e14]
	probeIntensities = [0.5e14, 1e14, 3e14, 4e14]
	controlDelays = [30]

	count = 0
	for m in molecules:
		for p in pulseDurations:
			for ci in controlIntensities:
				for pi in probeIntensities:
					if pi > ci:
						for cd in controlDelays:
							if cd == controlDelays[0] or ci > 0:
								args["molecule"] = m
								args["pulseDuration"] = p*femtosec_to_au
								args["pulseIntensity"] = pi
								args["controlDuration"]  = p*femtosec_to_au
								args["controlIntensity"] = ci
								args["controlDelay"] = cd*femtosec_to_au
								SubmitControlPumpExperiment(**args)
								count += 1

	print "Count= ", count
					

def SubmitAll2():
	args = {}
	args["delayList"] = r_[0:800:2]
	args["radialScaling"] = 2

	molecules = ["d2+", "hd+", "h2+"]
	pulseDurations = [5, 12, 18]
	probeIntensities = [0.5e14, 1e14, 4e14]
	phases = [0, pi/4., pi/2.]

	count = 0
	for m in molecules:
		for p in pulseDurations:
			for intensity in probeIntensities:
				for phase in phases:
					args["molecule"] = m
					args["pulseDuration"] = p*femtosec_to_au
					args["pulseIntensity"] = intensity
					args["pulsePhase"] = phase
					SubmitPhaseExperiment(**args)
					count += 1

	print "Count= ", count

def SubmitAllFinal():
	duration5fs = 5 * sqrt(2) * femtosec_to_au
	duration12fs = 12 * sqrt(2) * femtosec_to_au

	#5fs (intensity) 
	SubmitFinalExperiment(radialScaling=1, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration5fs, pulseIntensity=5e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration5fs, pulseIntensity=5e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration5fs, pulseIntensity=10e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration5fs, pulseIntensity=20e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration5fs, pulseIntensity=40e13, pulsePhase=0.0)
	
	#12fs (intensity) 
	SubmitFinalExperiment(radialScaling=1, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration12fs, pulseIntensity=5e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration12fs, pulseIntensity=5e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration12fs, pulseIntensity=10e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration12fs, pulseIntensity=20e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration12fs, pulseIntensity=40e13, pulsePhase=0.0)
	
	#CEP
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration5fs, pulseIntensity=20e13, pulsePhase=pi/2)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", pulseDuration=duration5fs, pulseIntensity=40e13, pulsePhase=pi/2)
	
	#h2+
	#5fs (intensity) 
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration5fs, pulseIntensity=5e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration5fs, pulseIntensity=10e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration5fs, pulseIntensity=20e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration5fs, pulseIntensity=40e13, pulsePhase=0.0)
	
	#12fs (intensity) 
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration12fs, pulseIntensity=5e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration12fs, pulseIntensity=10e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration12fs, pulseIntensity=20e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration12fs, pulseIntensity=40e13, pulsePhase=0.0)

	#CEP
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration5fs, pulseIntensity=20e13, pulsePhase=pi/2)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="h2+", pulseDuration=duration5fs, pulseIntensity=40e13, pulsePhase=pi/2)

	#hd+
	#5fs (intensity) 
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration5fs, pulseIntensity=5e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration5fs, pulseIntensity=10e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration5fs, pulseIntensity=20e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration5fs, pulseIntensity=40e13, pulsePhase=0.0)

	#12fs (intensity) 
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration12fs, pulseIntensity=5e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration12fs, pulseIntensity=10e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration12fs, pulseIntensity=20e13, pulsePhase=0.0)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration12fs, pulseIntensity=40e13, pulsePhase=0.0)

	#CEP
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration5fs, pulseIntensity=20e13, pulsePhase=pi/2)
	SubmitFinalExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="hd+", pulseDuration=duration5fs, pulseIntensity=40e13, pulsePhase=pi/2)

def SubmitAllFinalControl():
	duration5fs = 5 * sqrt(2) * femtosec_to_au
	duration12fs = 12 * sqrt(2) * femtosec_to_au

	"""
	#5fs (intensity) 
	#26.5fs delay 
	SubmitFinalControlExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", controlDuration=duration5fs, \
		controlIntensity=20e13, controlDelay=26.5*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=30e13)
	SubmitFinalControlExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", controlDuration=duration5fs, \
		controlIntensity=20e13, controlDelay=26.5*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=40e13)

	#493fs delay
	SubmitFinalControlExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", controlDuration=duration5fs, \
		controlIntensity=20e13, controlDelay=493*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=30e13)
	SubmitFinalControlExperiment(radialScaling=2, delayList=r_[0:800:0.25], molecule="d2+", controlDuration=duration5fs, \
		controlIntensity=20e13, controlDelay=493*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=40e13)
	"""
	#293fs delay (even states)
	#SubmitFinalControlExperiment(radialScaling=2, delayList=r_[310:3000:1], molecule="d2+", controlDuration=duration5fs, controlIntensity=5e13, controlDelay=293*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=30e13)
	SubmitFinalControlExperiment(radialScaling=2, delayList=r_[310:4000], molecule="d2+", controlDuration=duration5fs, controlIntensity=5e13, controlDelay=293*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=40e13)

	#306fs delay (odd states)
	#SubmitFinalControlExperiment(radialScaling=2, delayList=r_[310:3000:1], molecule="d2+", controlDuration=duration5fs, controlIntensity=5e13, controlDelay=306*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=30e13)
	SubmitFinalControlExperiment(radialScaling=2, delayList=r_[310:4000], molecule="d2+", controlDuration=duration5fs, controlIntensity=5e13, controlDelay=306*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=40e13)

	SubmitFinalControlExperiment(radialScaling=2, delayList=r_[0:4000:1], molecule="d2+", controlDuration=duration5fs, \
		controlIntensity=20e13, controlDelay=22*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=40e13)
	SubmitFinalControlExperiment(radialScaling=2, delayList=r_[0:4000:1], molecule="d2+", controlDuration=duration5fs, \
		controlIntensity=20e13, controlDelay=24.5*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=40e13)
	SubmitFinalControlExperiment(radialScaling=2, delayList=r_[0:4000:1], molecule="d2+", controlDuration=duration5fs, \
		controlIntensity=20e13, controlDelay=26.2*femtosec_to_au, pulseDuration=duration5fs, pulseIntensity=40e13)

def SubmitFinalControlExperiment(**args):
	args["outputfile"] = "outputfiles/%s/final_control_%s_control_%.2ffs_%ie13_probe_%ifs_%ie13_scaling_%i.h5" % \
		( \
		args["molecule"], \
		"%i", \
		args["controlDelay"]/femtosec_to_au, \
		args["controlIntensity"]/1e13, \
		args["pulseDuration"]/(femtosec_to_au*sqrt(2)), 
		args["pulseIntensity"]/1e13, \
		args["radialScaling"] \
		)

	print args["outputfile"]

	#Submit:
	args["partitionCount"] = args.get("partitionCount", 64)
	SubmitDelayScan(**args)

	##Make Plots
	#datafile = args["outputfile"].replace("%i", "all")
	#MakeDelayScanPlot(args["molecule"], datafile, radialScaling=args["radialScaling"], batchMode=True)
	#path, name = os.path.split(datafile)
	#root, ext = os.path.splitext(name)
	#pylab.savefig("figures/%s_%s.png" % (args["molecule"], root))
	#pylab.close("all")



def SubmitFinalExperiment(**args):
	args["outputfile"] = "outputfiles/%s/final_%s_phase_%.2fpi_pump_%ifs_%ie13_scaling_%i.h5" % \
		( \
		args["molecule"], \
		"%i", \
		args["pulsePhase"]/pi, \
		args["pulseDuration"]/(femtosec_to_au*sqrt(2)), 
		args["pulseIntensity"]/1e13, \
		args["radialScaling"] \
		)

	print args["outputfile"]

	#Submit:
	args["partitionCount"] = args.get("partitionCount", 64)
	SubmitDelayScan(**args)

	##Make Plots
	#datafile = args["outputfile"].replace("%i", "all")
	#MakeDelayScanPlot(args["molecule"], datafile, radialScaling=args["radialScaling"], batchMode=True)
	#path, name = os.path.split(datafile)
	#root, ext = os.path.splitext(name)
	#pylab.savefig("figures/%s_%s.png" % (args["molecule"], root))
	#pylab.close("all")


					
def SubmitPhaseExperiment(**args):
	args["outputfile"] = "outputfiles/%s/phase_%s_%.2fpi_pump_%ifs_%ie13.h5" % \
		( \
		args["molecule"], \
		"%i", \
		args["pulsePhase"]/pi, \
		args["pulseDuration"]/femtosec_to_au, \
		args["pulseIntensity"]/1e13 \
		)

	print args["outputfile"]
	args["partitionCount"] = args.get("partitionCount", 64)
	SubmitDelayScan(**args)



def SubmitControlPumpExperiment(**args):
	args["outputfile"] = "outputfiles/%s/control_%s_%ifs_%ifs_%ie13_pump_%ifs_%ie13.h5" % \
		( \
		args["molecule"], \
		"%i", \
		args["controlDuration"]/femtosec_to_au, \
		args["controlDelay"]/femtosec_to_au, \
		args["controlIntensity"]/1e13, \
		args["pulseDuration"]/femtosec_to_au, \
		args["pulseIntensity"]/1e13 \
		)

	print args["outputfile"]
	args["partitionCount"] = args.get("partitionCount", 64)
	SubmitDelayScan(**args)


def SubmitDelayScanIPython1(**args):
	delayList = args["delayList"]
	outputfile = args["outputfile"]
	molecule = args["molecule"]

	controllerHost = DefaultControllerHost
	controllerPort = DefaultControllerPort
	if "controllerHost" in args:
		controllerHost = args["controllerHost"]
	if "controllerPort" in args:
		controllerPort = args["controllerPort"]

	#Create connection to stallo
	print "Connecting to ipython1 controller..."
	rc = kernel.RemoteController((controllerHost, controllerPort))
	partitionCount = len(rc.getIDs())

	if partitionCount == 0:
		raise Exception("No engines connected to controller @ stallo.")

	#Make sure pyprop is loaded
	print "Loading pyprop..."
	rc.executeAll('import os')
	rc.executeAll('os.environ["PYPROP_SINGLEPROC"] = "1"')
	rc.executeAll('execfile("example.py")')

	#scatter delay list
	print "Distributing jobs..."
	rc.scatterAll("delayList", delayList)
	rc.scatterAll("partitionId", r_[:partitionCount])
	rc.pushAll(args=args)

	#run
	print "Running jobs..."
	rc.executeAll('args["delayList"] = delayList')
	rc.executeAll('args["outputfile"] = args["outputfile"] % partitionId[0]')
	rc.executeAll('RunDelayScan(**args)')

	#gather all files into one
	print "Gathering all outputfiles into one..."
	filenames = [outputfile % i for i in range(partitionCount)]
	combinedFile = outputfile.replace("%i", "all")
	rc.push(1, filenames=filenames, combinedFile=combinedFile)
	rc.execute(1, 'ConcatenateHDF5(filenames, combinedFile)')
	#remove all single proc files
	rc.executeAll('os.unlink(args["outputfile"])')

	print "Done."


import pyprop.utilities.submitpbs_stallo as submitpbs

def SubmitDelayScan(**args):
	delayList = args["delayList"]
	outputfile = args["outputfile"]
	molecule = args["molecule"]
	partitionCount = args["partitionCount"]

	partitionSize = int(ceil(len(delayList)/float(partitionCount)))
	partitionStep = delayList[1] - delayList[0] #Assume equidistant delayList

	plist = []
	
	for i in range(partitionCount):
		start = min(i*partitionSize, len(delayList))
		end = min((i+1)*partitionSize, len(delayList))
		delaySlice = slice(delayList[start-1]+partitionStep, delayList[end-1]+partitionStep, partitionStep)
		args["delayList"] = delaySlice
		args["outputfile"] = outputfile % i

		script = submitpbs.SubmitScript()
		script.interconnect = None
		script.executable = './run_delay_scan.py'
		script.parameters = commands.mkarg(repr(args))  
		script.ppn = 1
		script.nodes = 1
		script.walltime = submitpbs.timedelta(hours=2)
		print "\n".join(script.CreateScript())
		#script.Submit()



#------------------------------------------------------------------------------	
#                Propagate scan delay
#------------------------------------------------------------------------------	

def RunDelayScan(**args):
	#Required parameters
	delayList = args["delayList"]
	outputfile = args["outputfile"]
	molecule = args["molecule"]
	args["silent"] = True
	args["configSilent"] = True

	print "USING OUTPUTFILE ", outputfile

	for delay in delayList:
		args["pulseDelay"] = delay*femtosec_to_au
		args["outputpath"] = "/delay_%.4f" % (delay)
		print "Propagating %s with pulse delay %ffs" % (molecule, delay)
		PropagateDelayScan(**args)

def PropagateDelayScan(**args):
	"""
	Propagates one run for the scan pulse delay experiment, and
	saves the result to disk.

	The following arguments are required
	molecule: The molecule to simulate "hd+" or "d2+" or "h2+"
	outputfile: The hdf5 file to save the results to (i.e. pulsedelay_30.h5)
	outputpath: The group inside outputfile where (i.e. "/delay_30fs")

	"""	
	args['config'] = "config.ini"
	args['imtime'] = False
	inputfile = GetInputFile(**args)
	outputfile = args["outputfile"]
	outputpath = args["outputpath"]

	conf = SetupConfig(**args)
	prop = SetupProblem(**args)

	pumpFrequency = conf.InitialCondition.pump_frequency
	pumpTransitionProbability = conf.InitialCondition.pump_probability
	pumpCount = conf.InitialCondition.pump_count
	pumpStateEnergy = GetInitialStateEnergy(inputfile)
	pumpTimes, pumpProb, pumpPhase = GetPump(pumpFrequency, pumpCount, pumpTransitionProbability, pumpStateEnergy)

	LoadInitialState(prop, **args)
	initPsi = prop.psi.Copy()

	boundE, boundV = LoadBoundEigenstates(**args)
	contE1, contV1, contE2, contV2 = LoadContinuumEigenstates(**args)

	r = prop.psi.GetRepresentation().GetLocalGrid(0)
	timeList = []
	initCorrList = []
	corrList = []
	normList = []
	psiTimeList = []
	timeList = []
	energyDistribution = []

	dE = average(diff(contE2))
	E = r_[-0.5:0:dE]

	pulseStart = conf.ProbePulsePotential.delay
	pulseDuration = conf.ProbePulsePotential.duration
	
	output = tables.openFile(outputfile, "a")
	try:
		RemoveNodeIfExists(output, outputpath, "wavefunction")
		atom = tables.ComplexAtom(16)
		shape = (0,) + prop.psi.GetData().shape
		psiList = output.createEArray(outputpath, "wavefunction", atom, shape, createparents=True)
		
		def checkpoint():
			timeList.append(prop.PropagatedTime)
			initCorrList.append(prop.psi.InnerProduct(initPsi))
			normList.append(prop.psi.GetNorm())
			psiTimeList.append(prop.PropagatedTime)
			psiList.append(prop.psi.GetData().reshape((1,) + prop.psi.GetData().shape))
			energyDistribution.append(array(CalculateEnergyDistribution(prop.psi.GetData(), E, contE1, contV1, contE2, contV2)))
			corrList.append(abs(dot(boundV, prop.psi.GetData()[:,0]))**2)
		
		#-1) Pump in the initial wavepacket
		prop.psi.GetData()[:] = 0
		for i in range(len(pumpTimes)):
			#Add a franck-condon with the correct phase wavepacket to our wavefunction.
			prop.psi.GetData()[:] += pumpProb[i] * pumpPhase[i] * initPsi.GetData()
			if i<len(pumpTimes)-1:
				#propagate until the next pump time
				prop.Duration = pumpTimes[i+1]
				for t in prop.Advance(False): pass
		#Normalize after the pumping is complete
		prop.psi.Normalize()

		#0) Store the initial packet
		#checkpoint()

		#1) Propagate until the pulse starts
		prop.Duration = pulseStart - 2*pulseDuration
		for t in prop.Advance(minimum(20, prop.Duration)): 
			#checkpoint()
			pass
		
		#1a) Save the wavepacket before the pulse to see how much has flowed out
		#checkpoint()	
		
		#2) Propagate until the end of the pulse
		prop.Duration = pulseStart + 2*pulseDuration
		index = 0
		for t in prop.Advance(True):
			#checkpoint()
			pass
		
		#2a) Save the wavepacket at the end 
		checkpoint()
		
		SaveArray(output, outputpath, "time", array(timeList))
		SaveArray(output, outputpath, "eigenstateProjection", array(corrList))
		SaveArray(output, outputpath, "initstateProjection", array(initCorrList))
		SaveArray(output, outputpath, "norm", array(normList))
		SaveArray(output, outputpath, "timeWavefunction", array(psiTimeList))
		SaveArray(output, outputpath, "energyDistribution", array(energyDistribution))
		SaveArray(output, "/", "energy", E)
		psiList.close()

	finally:
		output.close()

	return array(timeList), array(corrList)

