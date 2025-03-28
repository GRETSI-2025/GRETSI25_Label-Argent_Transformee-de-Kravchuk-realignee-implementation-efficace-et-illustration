# LOAD NECESSARY PYTHON LIBRARIES

import numpy         as np
import scipy.special as sps

def the_ptransform(x, p = 0.5, thetas = np.linspace(1e-10,np.pi,500), phis = np.linspace(0,2*np.pi,500)):

    """
    Compute the p-Kravchuk transform of a discrete signal.

    Args:
        - x (numpy.ndarray): discrete signal, noisy or not, possibly complex valued.
        - p (float, optional): parameter of the transform between 0 and 1 (default 0.5 corresponding to the original transform).
        - thetas (numpy.ndarray, optional): polar angles at which to compute the transform.
        - phis (numpy.ndarray, optional): azimuthal angles at which to compute the transform.

    Returns:
        - Ks (numpy.ndarray): p-Kravchuk transform, complex-valued, evaluated at a discrete set of points on the sphere.
        - thetas (numpy.ndarray): polar angles at which the transform is computed.
        - phis (numpy.ndarray): azimuthal angles at which the transform is computed.
        
    """

    N  = x.shape[0] - 1
    q  = np.sqrt(p/(1-p))
    
    Nt = len(thetas)
    Np = len(phis)
    
    # duplicate the signal
    X = np.kron(np.conj(x),np.ones((Np,Nt,1))).T
    
    # maps of angles 
    Thetas, Phis = np.meshgrid(thetas,phis,indexing='ij')
    
    # map of phase space variable
    z = np.cos(Thetas/2)/np.sin(Thetas/2)*np.exp(1j*Phis)
    Z = np.kron(z,np.ones((N+1,1,1)))
    
    # indices
    ell = np.arange(N+1)
    L = np.kron(ell,np.ones((Np,Nt,1))).T
    
    # compute the transform using the generative function trick
    normZ  = np.sqrt((1+q**2)*(1+np.abs(Z)**2))
    sumand = X*np.sqrt(sps.comb(N,L))*((q-Z)/(normZ))**L*((1+q*Z)/(normZ))**(N-L)
    Kz = np.sum(sumand,0)
    
    return Kz, thetas, phis


    



    
    