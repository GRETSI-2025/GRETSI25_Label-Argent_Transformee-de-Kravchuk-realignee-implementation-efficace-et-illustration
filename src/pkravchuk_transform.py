import numpy         as np
import scipy.special as sps

def the_ptransform(x, p = 0.5, thetas = np.linspace(1e-10,np.pi,500), phis = np.linspace(0,2*np.pi,500)):

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


    



    
    