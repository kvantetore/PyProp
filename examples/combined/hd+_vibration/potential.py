

def GetVibrationalPotential(psi, config, potential):
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
		potential[:] = [interp1.Evaluate(x) +  minimum(2./dr, 1./ abs(x)) for x in r]
	
	if psi.GetRank() == 2:
		potential[:,0] = [interp1.Evaluate(x) +  minimum(2./dr, 1./ abs(x)) for x in r]
		potential[:,1] = [interp2.Evaluate(x) +  minimum(2./dr, 1./ abs(x)) for x in r]

	#plot(r, potential)

def GetElectronicCouplingBates(psi, config):
	r = psi.GetRepresentation().GetLocalGrid(0)
	rho = (1 + abs(r) + r**2/3.) * exp(-abs(r)) 

	coupling = -1/(2 + 1.4*abs(r)) + r/(2*sqrt(1 - rho**2))
	coupling[where(r==0)] = 0
	return coupling


def GetElectronicCouplingBunkinTugov(psi, config):
	r = psi.GetRepresentation().GetLocalGrid(0)
	r_zero = 2.
	alpha = 0.72
	d = 1.07
	d_prime = 0.396
	xp = -0.055
	
	return d + d_prime*(abs(r)-r_zero) - alpha*xp*d_prime*(abs(r)-r_zero)**2

def GetElectronicCoupling(psi, config):
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
	E0 = 2744 * volt_per_metre_to_au * sqrt(config.intensity) 
	t0 = config.delay
	scaling = 4 * log(2) 

	envelope = exp( - scaling * (t - t0)**2 / config.duration**2 )
	field = E0 * cos(config.frequency*(t - t0)) * envelope

	return field
	
