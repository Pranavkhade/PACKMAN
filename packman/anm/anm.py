# -*- coding: utf-8 -*-
"""The 'ANM' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'ANM' object documentation for details.

Example::
    from packman.anm import ANM
    help( ANM )

Notes:
    * Tutorial: https://jerniganlab.github.io/Software/PACKMAN/Tutorials/compliance
    * For more details about the parameters for compliance, or to site this, read the following paper: https://doi.org/10.1002/prot.25968

Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.
    * Make sure that parameter to the ANM is changed from [float] to packman.molecule.atom

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
    * Part Credits: Ambuj Kumar (ambuj@iastate.edu)
"""

import numpy

from packman.molecule import Protein, Model, Atom

'''
##################################################################################################
#                                              ANM                                              #
##################################################################################################
'''

class ANM:
    """This class contains the functions essential to carry out the Anisotropic Network Model and Compliance analysis.

        Notes:
        * Tutorial: https://jerniganlab.github.io/Software/PACKMAN/Tutorials/compliance
        * For more details about the parameters for compliance, or to site this, read the following paper:

        Todo:
            * fix init; coords to atoms
            * Change the pf in such a way that it is not confusing
            * Add set_ functions
        
        Args:
            coords ([float])       : Two dimentional array of three dimentional points in the space.
            gamma (float, optional): Spring Constant Value.                                      Defaults to 1.0.
            dr (float, optional)   : Distance Cutoff.                                            Defaults to 15.0.
            power (int, optional)  : Power of distance (mainly useful in non-parametric mode).   Defaults to 0.
            pf (None, optional)    : Parameter free model. (Check the dr and power params)       Defaults to None.
        """
    
    def __init__(self, atoms, gamma=1.0, dr=15.0, power=0, pf=None):
        self.gamma   = gamma
        self.dr      = dr
        self.power   = power
        self.pf      = pf
        self.atoms   = [i for i in atoms]
        self.coords  = numpy.array([i.get_location() for i in self.atoms])
        if self.pf != None and self.pf <= 0:
            raise Exception("pf value cannot be zero or negative")
        if self.gamma <= 0:
            raise Exception("gamma value cannot be zero or negative")
        if self.dr <= 0:
            raise Exception("distance cutoff value cannot be zero or negative")

        self.fluctuations       = None
        self.stiffness_map      = None
        self.compliance_map     = None
        self.stiffness_profile  = None
        self.compliance_profile = None 


    '''Get Functions'''
    def get_hessian(self):
        """Get the Hessian Matrix of the ANM model.

        Notes:
            * Make sure that the ANM().calculate_hessian() is called before calling this function. (will return None otherwise)
        
        Returns:
            numpy.ndarray: Hessian matrix if successful; None otherwise
        """
        return self.hessian
    
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
    
    def get_stiffness_map(self):
        """Get the Stiffness Map obtained from Stiffness and Compliance Analysis
        
        Notes:
            * Make sure that the ANM().calculate_hessian(), ANM().calculate_decomposition() and ANM().calculate_stiffness_compliance() is called before calling this function. (will return None otherwise)
            * Stiffness=1/Compliance
        Returns:
            numpy.ndarray: Stiffness Map if successful; None otherwise
        """
        return self.stiffness_map
    
    def get_compliance_map(self):
        """Get the Compliance Map obtained from Stiffness and Compliance Analysis
        
        Notes:
            * Make sure that the ANM().calculate_hessian(), ANM().calculate_decomposition() and ANM().calculate_stiffness_compliance() is called before calling this function. (will return None otherwise)
            * Stiffness=1/Compliance
        Returns:
            numpy.ndarray: Compliance Map if successful; None otherwise
        """
        return self.compliance_map
    
    def get_stiffness_profile(self):
        """Get the Stiffness profile obtained from Stiffness and Compliance Analysis
        
        Notes:
            * Make sure that the ANM().calculate_hessian(), ANM().calculate_decomposition() and ANM().calculate_stiffness_compliance() is called before calling this function. (will return None otherwise)
            * Stiffness=1/Compliance
        Returns:
            numpy.ndarray: Stiffness profile if successful; None otherwise
        """
        return self.stiffness_profile
    
    def get_compliance_profile(self):
        """Get the Compliance profile obtained from Stiffness and Compliance Analysis
        
        Notes:
            * Make sure that the ANM().calculate_hessian(), ANM().calculate_decomposition() and ANM().calculate_stiffness_compliance() is called before calling this function. (will return None otherwise)
            * Stiffness=1/Compliance
        Returns:
            numpy.ndarray: Compliance profile if successful; None otherwise
        """
        return self.compliance_profile
    

    '''Calculate Functions'''
    def calculate_hessian(self):
        """Build the Hessian Matrix of the ANM model.

        This is the most essential step for ANM/ Compliance analysis.
        
        Notes:
            * Hessian matrix is built; use ANM().get_hessian() to obtain the hessian matrix.
        """
        n_atoms=len(self.coords)
        hessian=numpy.zeros((n_atoms*3, n_atoms*3), float)
        distance_mat=numpy.ones((n_atoms*3, n_atoms*3), float)
        for i in range(len(self.coords)):
            diff = self.coords[i+1:, :] - self.coords[i] 
            squared_diff = diff**2
            for j, s_ij in enumerate(squared_diff.sum(1)):
                if s_ij <= self.dr**2:
                    diff_coords = diff[j]
                    j = j + i + 1
                    derivative = numpy.outer(diff_coords, diff_coords)*(float(-self.gamma)/numpy.sqrt(s_ij)**(2+self.power))
                    hessian[i*3:i*3+3, j*3:j*3+3] = derivative
                    hessian[j*3:j*3+3, i*3:i*3+3] = derivative
                    hessian[i*3:i*3+3, i*3:i*3+3] = hessian[i*3:i*3+3, i*3:i*3+3] - derivative
                    hessian[j*3:j*3+3, j*3:j*3+3] = hessian[j*3:j*3+3, j*3:j*3+3] - derivative
                    #abs added to avoid negative numbers
                    d = numpy.sqrt(s_ij)
                    lobj = [[d,d,d],[d,d,d], [d,d,d]]
                    dmat = numpy.array(lobj)
                    distance_mat[i*3:i*3+3, j*3:j*3+3] = dmat

        if self.pf != None:
            hessian = numpy.divide(hessian, distance_mat)
        
        self.hessian=hessian
        return True
    
    def calculate_decomposition(self):
        """Decompose the Hessian Matrix of the ANM model.
        
        Note:
            Eigen values and Eigen Vectors are calculated. use ANM().get_eigenvalues() and ANM().get_eigenvectors() to obtain them.
        """
        self.eigen_values,self.eigen_vectors=numpy.linalg.eigh(self.hessian)
        return True

    def calculate_fluctuations(self,endmode=None):
        """Calculate the Fluctuations of the ANM model.

        The fluctualtions/ theoretical b-factors are calculated using this method.
        
        Note:
            - Fluctuations are calculated. use ANM().get_fluctuations() to obtain the fluctuations.
            - Endmode needs to be put in the code if and when required.
        """
        EVec=self.eigen_vectors.T
        mode_bfactors=[]
        for numi,i in enumerate(self.eigen_values[6:]):
            evec_row=EVec[numi+6]
            mode_bfactors.append([ float(evec_row[j]**2 + evec_row[j+1]**2 + evec_row[j+2]**2)/i for j in range(0,len(self.eigen_values),3)])
            
        mode_bfactors=numpy.array(mode_bfactors)
        self.fluctuations=[numpy.sum(i) for i in mode_bfactors.T]
        return True
    
    
    def calculate_stiffness_compliance(self):
        """Carry out the Stiffness and Compliance analysis of the ANM model.

        Citation:
        Scaramozzino, D., Khade, P.M., Jernigan, R.L., Lacidogna, G. and Carpinteri, A.
        (2020), Structural Compliance ‐ A New Metric for Protein Flexibility. Proteins. Accepted Author Manuscript.
        doi:10.1002/prot.25968
        
        Note:
            * Obtain the following properties by using functions followed by it:
                Stiffness Map     : ANM().get_stiffness_map()      
                Compliance Map    : ANM().get_compliance_map()
                Stiffness Profile : ANM().get_stiffness_profile()
                Compliance Profile: ANM().get_compliance_profile()
                
        """
        EVec=self.eigen_vectors.T
        hessian_inv= numpy.matmul(  numpy.matmul( EVec[6:].transpose() , numpy.diag(1/self.get_eigenvalues()[6:]) )  , EVec[6:]  )
        compliance_map=numpy.zeros((len(self.coords),len(self.coords)))
        stiffness_map=numpy.zeros((len(self.coords),len(self.coords)))

        for numi,i in enumerate(self.coords):
            for numj,j in enumerate(self.coords):
                if(numi!=numj):
                    FMAT=numpy.zeros((len(self.coords),3))
                    Fij=(j-i)/numpy.linalg.norm(i-j)
                    FMAT[numi]=-Fij
                    FMAT[numj]=Fij
                    FMAT=FMAT.flatten()
                    
                    d=numpy.matmul(hessian_inv,FMAT)
                    
                    compliance= (d[(numj*3)+0]-d[(numi*3)+0])*Fij[0] + (d[(numj*3)+1]-d[(numi*3)+1])*Fij[1] + (d[(numj*3)+2]-d[(numi*3)+2])*Fij[2]
                    stiffness= 1/ compliance

                    compliance_map[numi][numj]= compliance
                    stiffness_map[numi][numj] = stiffness
        
        self.stiffness_map      = stiffness_map
        self.compliance_map     = compliance_map
        self.stiffness_profile  = [numpy.nanmean(i) for i in stiffness_map]
        self.compliance_profile = [numpy.nanmean(i) for i in compliance_map]
        return True
    
    def calculate_movie(self, mode_number, scale=1.5, n=20, ftype='cif'):
        """Get the movie of the obtained LINEAR modes. The first frame is the original structure and the projection progresses in positive (+) direction, returning to original structure and then in negative direction (-) again returning to the original structure.

        Args:
            mode_number (int)                   : Mode number. (first non-rigid mode is 6th)
            scale (float)                       : Multiplier; extent to which mode will be extrapolated.                 Defaults to 1.5
            n (int)                             : Number of frames in output (should be =>8 and ideally multiple of 4)   Defaults to 20
            ftype (string)                      : Extension of the output file (.cif / .pdb)

        Note:
            - Scale and n parameters should be redesigned.
            - direction is the variable which be allow user to explore only positive or only negative direction of the modes.

        Returns:
            True if successful; false otherwise.
        """
        x0=self.coords
        new_coords=[]
    
        movement = [ numpy.sin(k*(1.0/float(n))*2*numpy.pi) for k in range(0,n+1,1) ]
        ModelsOfTheProtein = []
        for j in movement:
            AtomsOfTheFrame = {}
            for numi,i in enumerate(x0):
                new_x=i[0]+scale*j*self.eigen_vectors[:,mode_number][(numi*3)+0]
                new_y=i[1]+scale*j*self.eigen_vectors[:,mode_number][(numi*3)+1]
                new_z=i[2]+scale*j*self.eigen_vectors[:,mode_number][(numi*3)+2]
                new_x , new_y, new_z = new_x.real , new_y.real , new_z.real
                currentatm = Atom(self.atoms[numi].get_id() , self.atoms[numi].get_name(), numpy.array([new_x,new_y,new_z]), self.atoms[numi].get_occupancy(), self.atoms[numi].get_bfactor(), self.atoms[numi].get_element(),self.atoms[numi].get_charge(), self.atoms[numi].get_parent() )
                AtomsOfTheFrame[numi] = currentatm

            ModelsOfTheProtein.append( Model(j, AtomsOfTheFrame, None, None, None, None) )
        Annotations = self.atoms[0].get_parent().get_parent().get_parent().get_parent().get_data()
        prot = Protein(str(mode_number), ModelsOfTheProtein)
        prot.set_data(Annotations)
        prot.write_structure( str(mode_number)+'.'+ftype )            
        return True