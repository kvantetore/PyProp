execfile("setup_methods.py")
execfile("setup_potential.py")

def diagonalize(number_to_keep, config = "config.ini"):
    """
    E, wf, prop = diagonalize(number_to_keep, config = "config.ini")

    Diagonalizes the 1D hamiltonian supplied in the config file.

    Parametres
    ----------
    number_to_keep : int, number of states that are to be returned, from lowest energy.
    config : string, containing name of config file, defaultly set to "config.ini".

    Returns
    -------
    E : 1D array, containing the 'number_to_keep' first energies.
    wf : 2D array, containing 'number_to_keep' columns of wavefunctions.
    prop : the corresponding PyProp object.
    """

    #Creating config object.
    conf = pyprop.Load(config)

    #Generating pyprop object.
    prop = SetupProblem(conf)

    #Finding the Hamiltonian Matrix.
    H,S = construct_hamiltonian_matrix(prop, return_overlap = True)

    #Solving the eigenvalue problem (calling eig).
    #   Last parametre decides how many of the states that are kept.
    E, wf = Solve(H.real, S, prop, number_to_keep)

    return E, wf, prop

def diagonalize_to_grid(my_grid, number_to_keep, project_name, config = "config.ini"):
    """
    diagonalize_to_grid(my_grid, number_to_keep, project_name, config = "config.ini")
    
    Diagonalizes the 1D hamiltonian supplied in the config file. 
    Transforms it onto my_grid. Saves the result to disk in 'TXT/<project_name>'.

    Parametres
    ----------
    my_grid : 1D float array, containing the grid onto which the wavefunctions should be stored. 
    number_to_keep : int, number of states that are to be returned, from lowest energy.
    project_name : string, which gives name to the subfolder of TXT that should contain the result.
	If the folder does not exist, it is created.
    config : string, containing name of config file, defaultly set to "config.ini".

    """
    
    #Diagonalize.
    E, wf, prop = diagonalize(number_to_keep, config = config)

    #Get spline info.
    spline_object = prop.psi.GetRepresentation().GetRepresentation(0).GetBSplineObject() 

    #Initialise output array.
    WF = zeros([number_to_keep, len(my_grid)])

    for i in range(number_to_keep):
	#Assuring correct format for the transformation.
	temp = zeros([wf.shape[0], 2])
	temp[:, 0] = wf[:, i]
	temp.dtype = 'complex128'

	#Doing the transformation.
	WF[i, :] = spline_object.ConstructFunctionFromBSplineExpansion(temp[:, 0], my_grid)
    
    #Saving to txt-files.
    if not os.path.exists("TXT/%s"%project_name):
	#Creating folder.
	os.system("mkdir TXT/%s"%project_name)
    
    #Saving to mentioned folder.
    savetxt("TXT/%s/Energies.txt"%project_name, E)
    savetxt("TXT/%s/Wavefunctions.txt"%project_name, WF)
    savetxt("TXT/%s/Grid.txt"%project_name, my_grid)

