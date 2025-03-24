import numpy         as np
import scipy.special as sps
import pyssht        as ssht

def the_ssht_transform(x,method = 'MW',p = 1/2):
    '''
    Default method in SSH toolbox https://pypi.org/project/pyssht/ is 'MW' for McEwen & Wiaux sampling.
    '''

    N  = x.shape[0] - 1
    
    if  np.mod(N,2):
        print('Error: length of signal should be odd.')
    else:
        
        # largest index explored
        ell = int(N/2)

        # band-limit 
        L = ell + 1
        
        # spherical coordinates of the pixel centers used by the SSHT toolbox
        (thetas, phis) = ssht.sample_positions(L, Method=method)              # vectors
        (Thetas, Phis) = ssht.sample_positions(L, Method=method, Grid = True) # meshgrid
        thetas         = np.pi - thetas
        
        # store the spin-spherical harmonic coefficient
        flm  = np.zeros(L**2,dtype = 'complex128')

        # compute the angle beta associated to p
        beta = 2*np.arcsin(np.sqrt(p))

        for m in np.arange(-ell,ell+1):
            n          = ell+m 
            index      = ssht.elm2ind(ell,m)
            flm[index] = np.conj(x[n])
        
        slm        = ssht.rotate_flms(flm, 0,beta, 0,L)
        f          = ssht.inverse(slm,L,Spin = ell, Method=method)
        f          *= np.sqrt(4*np.pi/(N+1))
        f          *= np.exp(1j*N*Phis/2)
        f          *= (-1)**(N/2)

    return f,thetas, phis


def rotate_signal(x,theta = 0, phi = 0, method = 'MW'):
    '''
    Default method in SSH toolbox https://pypi.org/project/pyssht/ is 'MW' for McEwen & Wiaux sampling.
    '''

    N  = x.shape[0] - 1
    
    if  np.mod(N,2):
        print('Error: length of signal should be odd.')
    else:
        
        # largest index explored
        ell = int(N/2)

        # band-limit 
        L = ell + 1
        
        # store the spin-spherical harmonic coefficient
        flm  = np.zeros(L**2,dtype = 'complex128')

        for m in np.arange(-ell,ell+1):
            n          = ell+m 
            index      = ssht.elm2ind(ell,m)
            flm[index] = np.conj(x[n])
        
        slm        = ssht.rotate_flms(flm, phi,theta, 0,L)

        # store the reconstructed signal coefficients
        y = np.zeros((N+1,), dtype = 'complex128')
        for m in np.arange(-ell,ell+1):
            n          = ell+m 
            index      = ssht.elm2ind(ell,m)
            y[n]       = np.conj(slm[index])
        
    return y
    

def the_new_transform(x,method = 'MW'):
    '''
    Default method in SSH toolbox https://pypi.org/project/pyssht/ is 'MW' for McEwen & Wiaux sampling.
    '''

    N  = x.shape[0] - 1
    
    if  np.mod(N,2):
        print('Error: N+1 should be even')
    else:
        
        # largest index explored
        ell = int(N/2)

        # band-limit 
        L = ell + 1
        
        # spherical coordinates of the pixel centers used by the SSHT toolbox
        (thetas, phis) = ssht.sample_positions(L, Method=method)              # vectors
        (Thetas, Phis) = ssht.sample_positions(L, Method=method, Grid = True) # meshgrid
        
        # store the spin-spherical harmonic coefficients
        flm = np.zeros(L**2,dtype = 'complex128')
        for m in np.arange(-ell,ell+1):
            n          = ell+m 
            index      = ssht.elm2ind(ell,m)
            flm[index] = np.conj(x[n])

        # compute the inverse SSHT transform which coincides with the new Kravchuk transform
        f          = ssht.inverse(flm,L,Spin = ell, Method=method)

        # center the frequencies
        phis      += - np.pi
        
    return f,thetas, phis


def the_inverse_transform(f,N,method = 'MW'):
    '''
    Default method in SSH toolbox https://pypi.org/project/pyssht/ is 'MW' for McEwen & Wiaux sampling.
    '''

    if  np.mod(N,2):
        print('Error: N+1 should be even')
    else:
        
        # largest index explored
        ell = int(N/2)

        # band-limit 
        L = ell + 1
        
        # spherical coordinates of the pixel centers used by the SSHT toolbox
        (thetas, phis) = ssht.sample_positions(L, Method=method)              # vectors

        # compute the SSHT transform which coincides with the new Kravchuk inverse transform
        flm            = ssht.forward(f,L,Spin = ell, Method=method)

        # store the reconstructed signal coefficients
        x = np.zeros((N+1,), dtype = 'complex128')
        for m in np.arange(-ell,ell+1):
            n          = ell+m 
            index      = ssht.elm2ind(ell,m)
            x[n]       = np.conj(flm[index])
        
    return x, thetas, phis


def the_spherical_angles(N, method = 'MW'):
    '''
    Default method in SSH toolbox https://pypi.org/project/pyssht/ is 'MW' for McEwen & Wiaux sampling.
    '''
    
    # largest index explored
    ell = int(N/2)
    
    # band-limit 
    L   = ell + 1
    
    (thetas, phis) = ssht.sample_positions(L, Method=method)

    # adjust the azimuthal angle to fit the aligned transform definition
    phis           = phis - np.pi
    
    return thetas, phis
