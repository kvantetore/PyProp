
def Run(prop, numPlot):
	index = 0
	plotStep = int((prop.Config.Propagation.duration / prop.TimeStep) / numPlot)
	while prop.PropagatedTime < prop.Config.Propagation.duration:
		prop.AdvanceStep()
		index += 1
		if index % plotStep == 0:
			if IsMaster():
				print "t = ", prop.PropagatedTime
			Plot2DFull(prop)

