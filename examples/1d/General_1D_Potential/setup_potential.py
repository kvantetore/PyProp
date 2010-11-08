import numpy

class General1DPotential(object):
    """
    Class for using a general potential, given on a grid, with pyprop.
    This is slower than an C++ implementation, but decent for 1D problems.
    """

    def Setup(self, psi):
	pass


    def ApplyConfigSection(self, conf):
	"""
	Getting the information from the relevant config section.
	"""

	self.filename = conf.f_name
	self.do_interp = conf.do_interp


    def UpdatePotentialData(self, data, psi, t, dt):
	"""
	Returns the potential data on the integration grid.
	"""

	#If the potential is stored on a different grid from the integration grid.
	if(self.do_interp):	
	    
	    #Loading potential data, stored in <f_name>, given in the config section.
	    A = loadtxt(self.filename)

	    #Potential in 2nd column, grid in 1st column, please.
	    my_pot = A[:,1]
	    my_grid = A[:,0]
	    
	    #Interpolation
	    #-----------------
	    
	    #Finds integration grid.
	    spline_object = psi.GetRepresentation().GetRepresentation(0).GetBSplineObject()
	    in_grid = spline_object.GetQuadratureGridGlobal()
	    
	    #Interpolating potental onto this grid.
	    interp_pot = numpy.interp(in_grid, my_grid, my_pot.real)
	    #-----------------

	    data[:] = interp_pot
	
	
	#If the potential IS stored on the intergation grid.
	else: 
	    #Loading potential data, stored in <f_name>, given in the config section.
	    A = loadtxt(self.filename)
	    #Potential in 2nd column.
	    my_pot = A[:,1]
	    
	    #OBS: Assumes the saved data is on the integration grid!
	    if data.shape[0] != my_pot.shape[0]:
		raise Exception("Potential %s not compatible with grid %s"%(my_pot.shape, data.shape))
	    

	    data[:] = my_pot.real	


