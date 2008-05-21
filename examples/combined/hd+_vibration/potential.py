

def GetVibrationalPotential(psi, config, potential):
	"""
	Updates a StaticPotential with the interpolated values for the
	electronic energy states for diatomic ionic molecules as a function
	of internuclear separation.

	if config.species == "neutral", the electronic ground state energies for
	the neutral molecule are returned,
	if config.species == "ion", the electronic energies for the two lowest energies
	of the ionic molecule are returned
	"""


	#Decide which potential	to use
	species = config.species.lower()
	if species == "neutral":	
		grid = potential_data.GridNeutral
		data1 = potential_data.PotentialNeutral
		data2 = potential_data.PotentialNeutral
	elif species == "ion":
		grid = potential_data.GridIon
		data1 = potential_data.Potential_1s_sigma_G
		data2 = potential_data.Potential_2p_sigma_U
	else:
		raise Exception("Unknown species '%s'" % config.species)

	#mirror round r=0
	grid = array(list(-grid[::-1]) + list(grid))
	data1 = array(list(data1[::-1]) + list(data1))
	data2 = array(list(data2[::-1]) + list(data2))

	#Interpolate the potential
	interp1 = spline.Interpolator(grid, data1)
	interp2 = spline.Interpolator(grid, data2)

	r = psi.GetRepresentation().GetLocalGrid(0)
	dr = r[1] - r[0]
	if psi.GetRank() == 1:
		potentialSlope = 1
		if hasattr(config, "potential_slope"):
			potentialSlope = config.potential_slope
		if potentialSlope == 1:
			interp = interp1
		else:
			interp = interp2

		innerIndex = where(abs(r)<=max(grid))
		outerIndex = where(abs(r)>max(grid))
		potential[innerIndex] = [interp.Evaluate(x) +  minimum(2./dr, 1./ abs(x)) for x in r[innerIndex]]
		potential[outerIndex] = [-1 for x in r[outerIndex]]

	if psi.GetRank() == 2:
		potential[:,0] = [interp1.Evaluate(x) +  minimum(2./dr, 1./ abs(x)) for x in r]
		potential[:,1] = [interp2.Evaluate(x) +  minimum(2./dr, 1./ abs(x)) for x in r]

	#plot(r, potential)

def GetElectronicCouplingBates(psi, config):
	"""
	Returns the coupling between the 1s sigma g state and the 
	2p sigma u electronic states, as modeled by Bates (1951).

	Implementation is translated to python from Dohmnalls code
	"""
	r = psi.GetRepresentation().GetLocalGrid(0)
	rho = (1 + abs(r) + r**2/3.) * exp(-abs(r)) 

	coupling = -1/(2 + 1.4*abs(r)) + r/(2*sqrt(1 - rho**2))
	coupling[where(r==0)] = 0
	return coupling


def GetElectronicCouplingBunkinTugov(psi, config):
	"""
	Returns the coupling between the 1s sigma g state and the 
	2p sigma u electronic states, as modeled by Bunkin and Tugov.

	Implementation is translated to python from Dohmnalls code
	"""
	r = psi.GetRepresentation().GetLocalGrid(0)
	r_zero = 2.
	alpha = 0.72
	d = 1.07
	d_prime = 0.396
	xp = -0.055
	
	return d + d_prime*(abs(r)-r_zero) - alpha*xp*d_prime*(abs(r)-r_zero)**2

def GetElectronicCoupling(psi, config):
	"""
	Returns the timeindependent part of the electronic coupling potential: the matrix
	elements connecting the binding and unbinding electronic states.

	Two models for the coupling have been implemented, one by Bates(1951), and another
	by Bunkin and Tugov (...). Currently, only the Bates model is used, the otherone is implemented
	because it was available in Dohmnalls code
	"""
	if config.coupling_model == "Bates":
		coupling = GetElectronicCouplingBates(psi, config)
	elif config.coupling_model == "BunkinTugov":
		coupling = GetElectronicCouplingBunkinTugov(psi, config)
	else:
		raise "Invalid coupling model %s" % config.coupling_model

	shape = psi.GetData().shape
	if psi.GetRank() == 1:
		raise Exception("Must be a 2D wavefunction in order to use the coupling")
	if shape[1] != 2:
		raise Exception("Electronic coupling rank (2) must be of size 2")

	r = psi.GetRepresentation().GetLocalGrid(0)
	matrix = zeros((shape[0], 2, 2), dtype=complex)
	matrix[:,0,0] = config.static_dipole_moment * abs(r) / 2.0
	matrix[:,0,1] = coupling
	matrix[:,1,0] = conj(coupling)
	matrix[:,1,1] = config.static_dipole_moment * abs(r) / 2.0

	return matrix
	
def GetLaserField(config, t):
	"""
	Returns the timedependent part of the electronic coupling potential: the scalar value
	of the laser field at a time t
	"""
	E0 = 2744 * volt_per_metre_to_au * sqrt(config.intensity) 
	t0 = config.delay
	scaling = 4 * log(2) 
		
	envelope = exp( - scaling * (t - t0)**2 / config.duration**2 )
	field = E0 * cos(config.frequency*(t - t0) + config.phase) * envelope
	return field

def GetSquareField(config, t):
	"""
	Timedependent part for a square pulse
	"""
	E0 = 2744 * volt_per_metre_to_au * sqrt(abs(config.intensity)) * sign(config.intensity)
	t0 = config.delay
	w = config.duration
	if t0 - w/2 <= t < t0 + w/2:
		field = E0
	else:
		field = 0
	return field
	
def GetLaserPeaks(w, phase, t0, halfwidth):
	"""
	Returns the time for the significant peaks of a pulse given 
	frequency w, phase, delay t0 and halfwidth

	NOT IN USE
	"""
	phaseTime = - phase / w

	#The separation (in time) between two peaks. We need both positive
	#and negative peaks, so we use pi instead of 2*pi here
	peakSeparation = pi / w
	halfwidthCyclesPositive = floor((halfwidth/2 - phaseTime)  / peakSeparation)
	halfwidthCyclesNegative = floor((halfwidth/2 + phaseTime)  / peakSeparation)
	
	start = -halfwidthCyclesNegative * peakSeparation + phaseTime + t0
	count = halfwidthCyclesNegative + halfwidthCyclesPositive

	peakTimes = start + r_[0:count]*peakSeparation
	return peakTimes

def GetPump(laserFrequency, peakCount, transitionProbability, pumpstateEnergy):
	"""
	peakTimes, peakProbabilities = GetPump(...) 

	Returns the times and transition probabilities for pumping with a laser
	at a frequency laserFrequency, with peakCount number of peaks (counting
	both positive and negative peaks). Each transition has an probability
	given by transitionProbability


	"""
	peakSeparation = pi / laserFrequency

	peakTimes = r_[0:peakCount]*peakSeparation
	#Modify the transition probability with the probability for transition
	#to have occured at one of the earlier peaks
	peakProbabilities = zeros(peakCount, dtype=double)
	for i in range(peakCount):
		peakProbabilities[i] = (1 - sum(peakProbabilities)) * transitionProbability
	peakPhases = exp(-1.0j * pumpstateEnergy * peakTimes)
	
	return peakTimes, peakProbabilities, peakPhases

