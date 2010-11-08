#==============================================================================
#General methods
#==============================================================================

import sys
import os
import scipy.linalg

sys.path.append("./pyprop")
import pyprop
#All of this necessary?
pyprop = reload(pyprop)
pyprop.ProjectNamespace = globals()
#

from pyprop import PrintOut
from pylab import *
from numpy import *
from libpotential import *


def SetupProblem(conf):
    """
    prop = SetupProblem(conf)

    Makes a pyprop propagation object.

    Parametres
    ----------
    conf : config object, created using pyprop.Load("some_config_file.ini")

    Returns
    -------
    prop : pyprop propagation object.

    Example
    -------
    prop = SetupProblem(conf)
    for t in prop.Advance(50)
	print prop.PropegationTime
    """
    prop = pyprop.Problem(conf)
    prop.SetupStep()

    return prop

def construct_hamiltonian_matrix(prop,return_overlap = False):
    """
    H, (S) = construct_hamiltonian_matrix(prop)

    Returns the full Hamiltonian matrix (H), which is implicit in the pyprop object (prop).
    Can also return the overlap matrix, if prompted.

    Parametres
    ----------
    prop : pyprop propagator object
    return_overlap : boolean, if True, the overlap matrix, S is returned as well.

    Returns
    -------
    H : Hamiltonian matrix
    S : Overlap matrix

    Notes
    -----
    The program assumes that the entire potential is contained in PotentialList[0],
    and that the problem is one dimensional.

    Examples
    --------
    Setting up a problem, constructing the Hamiltonian, and finding the 
    eigenvalues and eigenvectors:
    
    >>> execfile('example.py')
    >>> prop = SetupProblem(config = 'config.ini')
    >>> H = construct_hamiltonian_matrix(prop)
    >>> w,v = numpy.linalg.eigh(H)

    """

    #Retrieving the potential.
    pot = prop.Propagator.BasePropagator.PotentialList[0]
    basis_pairs = pot.BasisPairs[0]
    zip_data = pot.PotentialData
 
    #Initializing the Hamiltonian matrix.
    nr_bsplines = pot.GeometryList[0].RankCount
    H = zeros([nr_bsplines,nr_bsplines],dtype=complex128)
   
    #Filling the Hamiltonian matrix.
    for i, (row,col) in enumerate(basis_pairs):
	H[row,col] = zip_data[i]
    
    if(return_overlap):
	#Setting up an "overlap potential"
	overlap = prop.Propagator.BasePropagator.GeneratePotential(prop.Config.OverlapMatrixPotential)
	overlap.SetupStep(0.)
	
	#Initializing the overlap matrix
	S = zeros(shape(H))
	
	basis_pairs = overlap.BasisPairs[0]
	zip_data = overlap.PotentialData

	for i, (row,col) in enumerate(basis_pairs):
	    S[row,col] = zip_data[i]

	return H, S
    #Returning the Hamiltonian.
    return H

def Solve(H, S, prop, nr, degenerated = False):
    """
    E,V = Solve(H, S, prop, nr, degenerated = False)
    
    Solves the generalized eigenvalue problem Hc = eSc, where H is the hamiltonian,
    S is the overlap matrix (between B-splines) and e is the eigenvalue. 
    Returns the eigenvalues, E, and eigenvectors, V.    

    Parametres
    ----------
    H : Hamiltonian matrix.
    S : Overlap matrix between B-splines.
    prop : pyprop object in which the B-splines are defined.
    nr : Integer. The number of states to be kept.
    degenerated : Boolean, defaultly False, determines if the wavefunction is 
	nonsymmetric, and therefor has to be transformed into symmetric and antisymmetric.

    Returns
    -------
    E : Vector containing the sorted eigenvalues.
    V : Array containing the normalized eigenvectors.

    Notes
    -----
    The degenerate function only works if the B-spline knotpoints are symmetrically 
    distributed about zero.
    """

    #Finding the eigenvalues and eigenvectors.
    E, V = scipy.linalg.eig(H, b=S)
    
    #Sorting them according to rising eigenvalues.
    I = argsort(E)
    E = E[I][:nr]
    V = V[:,I][:,:nr]
    
    if degenerated:
	#Temporary variable.
	U = zeros(nr)	
	for i in r_[:nr:2]:
	    #Symmetric:
	    #Fs =  (F1(x) + F2(x)) + (F1(-x) + F2(-x))
	    #Antisymmetric:
	    #Fa =  (F1(x) + F2(x)) - (F1(-x) + F2(-x))
	    U = V[:,i] + V[:,i+1]
	    V[:,i] = U + U[-1::-1] 
	    V[:,i+1] = U - U[-1::-1] 

    #Normalizing.
    V = V/normalisation_factor(V, S, prop)
    
    return E,V

def normalisation_factor(V, S, prop):
    """
    N = normalizaton_factor(V, S, prop)

    Calculates the inner product <fun,fun>, for all eigenfunctions.
    Simply divide the eigenfunctions by the factors N.

    Parametres
    ----------
    V : 2D array containing the unnormalized eigenvectors.
    S : Overlap matrix between B-splines.
    prop : pyprop object in which the B-splines are defined.

    Returns
    -------
    N : Vector containing the normalization factors.

    Examples
    --------
    >>> E, V = scipy.linalg.eig(self.H, b=self.S)
    >>> V = V/normalisation_factor(V, S, prop)
    """

    #B-spline basis
    nr_bsplines = prop.Propagator.BasePropagator.PotentialList[0].GeometryList[0].RankCount
    nr_used = V.shape[1]
    spline_order = prop.Config.BsplineRepresentation.order

    #Normalization constants
    N = zeros([1,nr_used])
    for i in range(nr_bsplines): 
	    for j in range(i,nr_bsplines):
		    #Only the nonzero integrands
		    if abs(i - j) < spline_order:    
			    #Exploiting the symmetry
			    if i == j:
				    N += V[i] * V[j] * S[i,j]
			    else:
				    N += 2 * V[i] * V[j] * S[i,j]
    return  sqrt(N) 



