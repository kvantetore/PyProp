
#-------------------------------------------------------------------------------------------
#                                         Chessboard Figure
#-------------------------------------------------------------------------------------------

def MakePosterAbstractFigure():
	#tControl, cControl = GetScanDelayCorrelation(outputfile="outputfiles/d2+/final_all_phase_0.00pi_pump_5fs_5e13_scaling_2.h5", partitionCount=0)

	rcParams["text.usetex"] = True
	rcParams["font.serif"] = ["Times New Roman"]
	rcParams["font.family"] = "serif"
	rcParams["figure.dpi"] = "80"
	rcParams["font.size"] = 10
	rcParams["axes.titlesize"] = "medium"
	rcParams["axes.labelsize"] = "medium"
	rcParams["xtick.labelsize"] = "small"
	rcParams["ytick.labelsize"] = "small"

	fig = figure(figsize=(17*cm_to_inch, 17./1.6*cm_to_inch))

	cmap = LoadColormap("gradient.txt", True)
	color = "#0045a2" #"#ff8e16" 
	maxState = 10

	c = cControl[::2,:maxState+1]
	t = tControl[::2]

	pcolormesh(t, arange(maxState+1), c.transpose(), cmap=cmap, shading="flat", vmin=0, vmax=0.35)
	axis((0,650,0,maxState))
	title("Vibrational Distribution")
	ylabel("Vibrational State")
	xlabel("Delay time (fs)")
	a = axis()
	cb = colorbar(ticks=[0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35])
	cb.ax.set_position([0.91, 0.16, 0.014, 0.73]);
	fig.subplots_adjust(left=0.08, bottom=0.16, right=0.89, top=0.89)
	draw()

	fig.savefig("poster/abstractfig.png", dpi=300)



