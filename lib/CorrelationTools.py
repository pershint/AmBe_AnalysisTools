#This library contains functions relevant to producing the correlation matrix
#Based on the phi coefficient.  The phi coefficient gives you a measure
#of how correlated two binary values are.
import playDarts as pd

import ROOT
from ROOT import TChain, gDirectory, gROOT

import copy
import numpy as np
import scipy
import scipy.linalg
import numpy.linalg

class CovarianceMatrix(object):
    def __init__(self, cov_matrix=None, variables=None):
        '''
        Class takes in a covariance matrix and the variables used to generate
        it (listed from row 0 to row N, and col 0 to col N).
        Tools included for Cholesky decomposition or SVD decomposition, then random
        shooting pass/fails with correlations included.
        '''

        self.cov_matrix = cov_matrix
        self.correlation_shooter_matrix = None #Cholesky decomposed or SVD
        self.variables = variables

    def CholeskyDecompose(self):
        dimension = len(np.diagonal(self.cov_matrix))
        ain = scipy.array(self.cov_matrix)
        #icorr = (0.4 * np.identity(dimension))
        #print(icorr)
        #ain = ain + icorr
        eigenvals = scipy.linalg.eigvalsh(ain)
        print("EIGENVALUES: " + str(eigenvals))
        for val in eigenvals:
            if val < 0:
                print("matrix must be positive definite to cholesky" +\
                        "decompose.  If they're close, consider shifting a bit")
                #return None
        c = scipy.linalg.cholesky(ain, lower=True)
        #u = scipy.linalg.cholesky(ain, lower=false)
        self.correlation_shooter_matrix = c

    def SVDDecompose(self):
        dimension = len(np.diagonal(self.cov_matrix))
        ain = scipy.array(self.cov_matrix)
        U, V, S = numpy.linalg.svd(ain)
        self.correlation_shooter_matrix = U   #TODO: Do we only need U to random shoot probabilities?

    def ShootCorrVars(self,numsets=1):
        fired_variables = None
        if self.variables is None:
            print("Need to give the variable list to associate shots with\n" + \
                    "variables in the output dictionary.")
            return
        #First, shoot random numbers from a normal distribution
        fired_norms = pd.RandShoot(0,1,len(self.variables)*numsets)
        fired_norms = np.array(np.split(fired_norms,numsets))
        #now, multiply your cholesky decomposition to add correlations
        corr_vectors = []
        #FIXME: SLOW... NEED TO DO DOT PRODUCT ACROSS MULTIPLE VECTORS AT ONCE
        for i in xrange(len(fired_norms)):
            corr_vectors.append(np.dot(np.array(self.correlation_shooter_matrix),\
                                fired_norms[i]))
        return corr_vectors

