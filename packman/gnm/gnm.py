# -*- coding: utf-8 -*-
"""The 'GNM' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'GNM' object documentation for details.

Example::
    from packman.GNM import GNM
    help( GNM )

Notes:
    * Tutorial:

Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
    * Part Credits: Ambuj Kumar (ambuj@iastate.edu)
"""

import numpy
import logging

'''
##################################################################################################
#                                              GNM                                              #
##################################################################################################
'''

class GNM:
    """This class contains the functions essential to carry out the Gaussian Network Model.

        Notes:
        * Tutorial:

        Todo:
            * Add set_ functions
        
        Args:
            coords ([float])       : Two dimentional array of three dimentional points in the space.
            gamma (float, optional): Spring Constant Value.                                      Defaults to 1.0.
            dr (float, optional)   : Distance Cutoff.                                            Defaults to 7.3 (Yang et. al., Protein elastic network models and the ranges of cooperativity. (2009))
            power (int, optional)  : Power of distance (mainly useful in non-parametric mode).   Defaults to 0.
    """
    
    def __init__(self, atoms, gamma=1.0, dr=7.3, power=0):
        self.gamma   = gamma
        self.dr      = dr
        self.power   = power
        self.atoms   = [i for i in atoms]
        self.coords  = numpy.array([i.get_location() for i in self.atoms])
        if self.gamma <= 0:
            raise Exception("gamma value cannot be zero or negative")
        if self.dr <= 0:
            raise Exception("distance cutoff value cannot be zero or negative")

        self.fluctuations       = None

    '''Get Functions'''
    def get_kirchoff(self):
        """Get the Hessian Matrix of the ANM model.

        Notes:
            * Make sure that the GNM().calculate_kirchoff() is called before calling this function. (will return None otherwise)
        
        Returns:
            numpy.ndarray: Hessian matrix if successful; None otherwise
        """
        return self.kirchoff
    
    def get_eigenvalues(self):
        """Get the Eigenvalues obtained by decomposing the Hessian Matrix of the ANM model.
        
        Notes:
            * Make sure that the ANM().calculate_hessian() and ANM().calculate_decomposition() is called before calling this function. (will return None otherwise)

        Returns:
            numpy.ndarray: Eigenvalues if successful; None otherwise
        """
        return self.eigen_values
    
    def get_eigenvectors(self):
        """Get the Eigenvectors obtained by decomposing the Hessian Matrix of the ANM model.
        
        Notes:
            * Make sure that the ANM().calculate_hessian() and ANM().calculate_decomposition() is called before calling this function. (will return None otherwise)

        Returns:
            numpy.ndarray: Eigenvectors if successful; None otherwise
        """
        return self.eigen_vectors
    
    def get_fluctuations(self):
        """Get the Fluctuations obtained from Eigenvectors and Eigenvalues
        
        Notes:
            * Make sure that the ANM().calculate_hessian(), ANM().calculate_decomposition() and ANM().calculate_fluctuations() is called before calling this function. (will return None otherwise)

        Returns:
            numpy.ndarray: Eigenvectors if successful; None otherwise
        """
        return self.fluctuations

    '''Calculate Functions'''
    def calculate_kirchoff(self, gamma = 1.0):
        """Calculate the Gaussian Network Model (GNM) Kirchoff Matrix.

        The matrix is stored in the self.GNM_MAT variable.
        
        Returns:
            True if successful; None otherwise.
        """
        n_atoms=len(self.coords)
        self.kirchoff = numpy.zeros((n_atoms, n_atoms), float)
        distance_mat = numpy.ones((n_atoms, n_atoms), float)
        for i in range(len(self.coords)):
            diff = self.coords[i+1:, :] - self.coords[i]
            squared_diff = diff**2
            for j, s_ij in enumerate(squared_diff.sum(1)):
                if s_ij <= self.dr**2:
                    diff_coords = diff[j]
                    j = j + i + 1
                    self.kirchoff[i, j] = - gamma
                    self.kirchoff[j, i] = - gamma
                    self.kirchoff[i, i] = self.kirchoff[i, i] + gamma
                    self.kirchoff[j, j] = self.kirchoff[j, j] + gamma
        
        return True
    
    def calculate_decomposition(self):
        """Decompose the Hessian Matrix of the ANM model.
        
        Note:
            Eigen values and Eigen Vectors are calculated. use ANM().get_eigenvalues() and ANM().get_eigenvectors() to obtain them.
        """
        self.eigen_values,self.eigen_vectors=numpy.linalg.eigh(self.kirchoff)
        return True

    def calculate_fluctuations(self, endmode=None):
        """Calculate the Fluctuations of the ANM model.

        The fluctualtions/ theoretical b-factors are calculated using this method.
        
        Note:
            - Fluctuations are calculated. use ANM().get_fluctuations() to obtain the fluctuations.
            - Endmode needs to be put in the code if and when required.
        """
        '''
            inverse=inverse+(float(1)/eigen_values[i])*eigen_vectors[:,i]*eigen_vectors[:,i].transpose()
        return inverse.diagonal()
        '''
        #Initiate
        if(endmode==None):
            stop_at = len(self.eigen_values)
        else:
            try:
                stop_at = int(endmode)
            except:
                logging.warning('Please provide valid input for the "endmode" parameter.')


        pseduinverse = numpy.zeros(shape=(len(self.eigen_values),len(self.eigen_values)))
        for i in range(1, stop_at):
            pseduinverse = pseduinverse + ( (float(1) / self.eigen_values[i])*self.eigen_vectors[:,i]*self.eigen_vectors[:,i].transpose() )

        self.fluctuations = pseduinverse.diagonal()
        return True
    
    def calculate_crosscorrelation(self):
        """Calculate the cross-correlation. (Read the paper for more details)
        Returns:
            True if successful; None otherwise.
        """
        n = len(self.atoms)
        EVec=self.eigen_vectors.T
        self.hessian_inv= numpy.matmul( numpy.matmul( EVec[1:].transpose() , numpy.diag(1/self.eigen_values[1:]) ) , EVec[1:] )
        self.C = numpy.zeros((n , n))
        for i in range(n):
            for j in range(n):
                self.C[i, j] = float(self.hessian_inv[i, j])/numpy.sqrt(self.hessian_inv[i, i]*self.hessian_inv[j, j])
                
        return True