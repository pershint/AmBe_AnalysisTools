#Program takes in a spectrum and information on the binning.  The code
#Then plays darts, picking random numbers on the spectrum's defined
#region.  If the number shot is below the spectrum, Store it in an array as
#an event.  When the number of events requested is collected, build the
#Experiment's histogram.

import numpy as np

def RandShoot_g0(mu, sigma,n):
    '''
    Returns an array of n numbers from a gaussian distribution of
    average mu and variance sigma. If less than zero, reshoot number.
    '''
    result = mu + (sigma * np.random.randn(n))
    for j,num in enumerate(result):
        if num < 0.0:
            isNegative = True
            while isNegative:
                result[j] = mu + (sigma * np.random.randn(1))
                if result[j] >= 0.0:
                    isNegative = False
    return result

def RandShoot(mu, sigma,n):
    '''
    Returns an array of n numbers from a gaussian distribution of
    average mu and variance sigma.
    '''
    result = mu + (sigma * np.random.randn(n))
    return result

def RandShoot_p(lamb, n):
    '''
    Returns an array of n numbers from a poisson distribution of
    average value lamb.
    '''
    result = np.random.poisson(lamb, n)
    return result
