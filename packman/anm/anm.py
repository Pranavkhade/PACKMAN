'''
Author: Pranav Khade (pranavk@iastate.edu)
'''

import numpy


'''
##################################################################################################
#                                              ANM                                              #
##################################################################################################
'''

class ANM:
    def __init__(self, coords , gamma=1.0, dr=15.0, power=0, pf=None):
        """
        Author: Pranav Khade
        Part Credits: Ambuj Kumar (ambuj@iastate.edu)
        """
        self.gamma   = gamma
        self.dr      = dr
        self.power   = power
        self.pf      = pf
        self.coords  = numpy.array(coords)
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
        return self.hessian
    
    def get_eigenvalues(self):
        return self.eigen_values
    
    def get_eigenvectors(self):
        return self.eigen_vectors
    
    def get_fluctuations(self):
        return self.fluctuations
        
    def get_fluctuations(self):
        return self.fluctuations
    
    def get_stiffness_map(self):
        return self.stiffness_map
    
    def get_compliance_map(self):
        return self.compliance_map
    
    def get_stiffness_profile(self):
        return self.stiffness_profile
    
    def get_compliance_profile(self):
        return self.compliance_profile
    

    '''Calculate Functions'''
    def calculate_hessian(self):
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
        '''
        '''
        self.eigen_values,self.eigen_vectors=numpy.linalg.eigh(self.hessian)
        return True

    def calculate_fluctuations(self,endmode=None):
        '''
        '''
        self.eigen_vectors=self.eigen_vectors.T
        mode_bfactors=[]
        for numi,i in enumerate(self.eigen_values[6:]):
            evec_row=self.eigen_vectors[numi+6]
            mode_bfactors.append([ float(evec_row[j]**2 + evec_row[j+1]**2 + evec_row[j+2]**2)/i for j in range(0,len(self.eigen_values),3)])
            
        mode_bfactors=numpy.array(mode_bfactors)
        self.fluctuations=[numpy.sum(i) for i in mode_bfactors.T]
        return True
    
    
    def calculate_stiffness_compliance(self):
        '''
        '''

        hessian_inv=  numpy.matmul(  numpy.matmul( self.get_eigenvectors()[6:].transpose() , numpy.diag(1/self.get_eigenvalues()[6:]) )  , self.get_eigenvectors()[6:]  )
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
