# -*- coding: utf-8 -*-
"""The 'hdANM' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'hdANM' object documentation and the publication for details.

Citation::

    Pranav M. Khade, Domenico Scaramozzino, Ambuj Kumar, Giuseppe Lacidogna, Alberto Carpinteri, Robert L. Jernigan, hdANM: a new comprehensive dynamics model for protein hinges, Biophysical Journal, 2021, https://doi.org/10.1016/j.bpj.2021.10.017

About .hng File:
- hdANM requires the information about hinges and domains in the .hng format.
- Each column in the .hng file is tab delimited.
- Each row in the .hng file follows collowing pattern: 

Filename_ChainID    Domain/HingeId   ResidueStartPosition:ResidueEndPosition

Example of .hng file for PDBID 1EXR::

    1EXR_A  D1  1:70
    1EXR_A  H1  70:90
    1EXR_A  D2  90:148

Above format means that there are two domains (D1 and D2) separated by a hinge (H1). D1 stretches from residue 1 to 70; D2 stretches from 90 to 148 and hinge H1 is in the middle.

Example::

    from packman.anm import hdANM
    help( hdANM )

Notes:
    * Tutorial: https://py-packman.readthedocs.io/en/latest/tutorials/hdANM.html#tutorials-hdanm

Todo:
    * Finish optimizing the performance.
    * Add publication details in the Notes
    * Make sure that parameter to the ANM is changed from [float] to packman.molecule.atom

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""

import logging
from .. import molecule, Atom, Model, Protein
from ..constants import amino_acid_molecular_weight
from ..constants import atomic_weight
from ..utilities import load_hinge

import numpy
import itertools

from scipy.linalg import eig as scipy_eig
'''
##################################################################################################
#                                            hd-ANM                                              #
##################################################################################################
'''

class hdANM:
    """This class contains the functions essential to carry out the Hinge-Domain-Anisotropic Network Model and Compliance analysis.

        Notes:
        * Tutorial: 
        * For more details about the parameters for compliance, or to site this, read the following paper:

        Todo:
            * Read individual functions to know
        
        Args:
            atoms ([packman.molecule.Atom]) : Two dimentional array of atoms.
            hng_file (string)               : .hng filename and path. (Contains the information about hinge and domains on the protein)
            gamma (float, optional)         : Spring Constant Value.                                      Defaults to 1.0.
            dr (float, optional)            : Distance Cutoff.                                            Defaults to 15.0.
            power (int, optional)           : Power of distance (mainly useful in non-parametric mode).   Defaults to 0.
            pf (None, optional)             : Parameter free model?.                                      Defaults to None.
        
        Raises:
            Exception: [description]
            Exception: [description]
            Exception: [description]
        """
    
    def __init__(self, atoms , hng_file , gamma=1.0, dr=15.0, power=0, pf=None):
        self.gamma   = gamma
        self.dr      = dr
        self.power   = power
        self.pf      = pf
        self.atoms   = atoms
        self.HNGinfo = load_hinge(hng_file)
        self.RT_eigen_vectors = None

        #Coords are in the same order as the atoms
        self.coords  = numpy.array([i.get_location() for i in atoms])
        if self.pf != None and self.pf <= 0:
            raise Exception("pf value cannot be zero or negative")
        if self.gamma <= 0:
            raise Exception("gamma value cannot be zero or negative")
        if self.dr <= 0:
            raise Exception("distance cutoff value cannot be zero or negative")
        
        ##To Verify that the hessian is same as ANM (Part 1/2)
        #self.calculate_hessian()


    '''Get Functions'''
    def get_hessian(self):
        """Get the Hessian Matrix of the hd-ANM model.

        Notes:
            * Make sure that the hdANM().calculate_hessian() is called before calling this function. (will return None otherwise)
        
        Returns:
            numpy.ndarray: Hessian matrix if successful; None otherwise
        """
        return self.hessian
    
    def get_eigenvalues(self):
        """Get the Eigenvalues obtained by decomposing the Hessian Matrix of the hd-ANM model.
        
        Notes:
            * Make sure that the hdANM().calculate_hessian() and ANM().calculate_decomposition() is called before calling this function. (will return None otherwise)

        Returns:
            numpy.ndarray: Eigenvalues if successful; None otherwise
        """
        return self.eigen_values
    
    def get_eigenvectors(self):
        """Get the Eigenvectors obtained by decomposing the Hessian Matrix of the hd-ANM model.
        
        Notes:
            * Make sure that the hdANM().calculate_hessian() and hdANM().calculate_decomposition() is called before calling this function. (will return None otherwise)

        Returns:
            numpy.ndarray: Eigenvectors if successful; None otherwise
        """
        return self.eigen_vectors
    
    def get_RT_eigen_vectors(self):
        """Get the Reverse Transformed vectors from the hdANM eigenvectors. 
    
        Reverse transformed means that the hdANM eigenvector of dimension: 6D+3H x 6D+3H  (D: Number of domains; H: Number of hinge atoms) are converted to 3N x 6D+3H (N: Number of atoms)

        Note:
            - It was refered as 'exploded vector' utill version 1.3.3
    
        Returns:
            numpy.ndarray: Reverse Transformed Vector if successful; None otherwise
        """
        return self.RT_eigen_vectors
    
    def get_fluctuations(self):
        """Get the Fluctuations obtained from Eigenvectors and Eigenvalues
        
        Notes:
            * Make sure that the hdANM().calculate_hessian(), hdANM().calculate_decomposition() and hdANM().calculate_fluctuations() is called before calling this function. (will return None otherwise)

        Returns:
            numpy.ndarray: Fluctuations if successful; None otherwise
        """
        return self.fluctuations
    
    def get_hessian_pseudoinverse(self,n_modes="all"):
        """Get the Pseudoinverse of the hdANM model.

        Returns:
            numpy.ndarray: Pseudoinverse if successul; None otherwise
        """
        try:
            return self.hessian_pseudoinverse
        except:
            self.calculate_hessian_pseudoinverse(n_modes=n_modes)
            return self.hessian_pseudoinverse
    
    def get_crosscorrelation_matrix(self,n_modes="all"):
        try:
            return self.crosscorrelation_matrix
        except:
            self.calculate_cross_correlation(n_modes=n_modes)
            return self.crosscorrelation_matrix
    
    def get_hessian_block(self,Index1,Index2):
        """Calculate Hij (Hessian matrix component) using equation . ()

        Notes:
            
        Returns:
            numpy.ndarray: Hij if successful; None otherwise
        """
        i=self.coords[Index1]
        j=self.coords[Index2]
        dist_ij = numpy.linalg.norm(j-i)

        if(Index1!=Index2 and dist_ij<=self.dr):
            return numpy.array(
            [  [ (j[0]-i[0])*(j[0]-i[0]) , (j[0]-i[0])*(j[1]-i[1]) , (j[0]-i[0])*(j[2]-i[2]) ]  , 
            [ (j[1]-i[1])*(j[0]-i[0]) , (j[1]-i[1])*(j[1]-i[1]) , (j[1]-i[1])*(j[2]-i[2]) ]  ,
            [ (j[2]-i[2])*(j[0]-i[0]) , (j[2]-i[2])*(j[1]-i[1]) , (j[2]-i[2])*(j[2]-i[2]) ]  ]
            ) * ( -self.gamma/(dist_ij**2) )
        elif(Index1!=Index2 and dist_ij>self.dr):
            return numpy.zeros((3,3))
        elif(Index1==Index2):
            diag = numpy.zeros((3,3))
            for numk,k in enumerate(self.coords):
                dist_ik = numpy.linalg.norm(k-i)
                if(dist_ik<=self.dr and Index1!=numk):
                    diag=numpy.add(diag, 
                        numpy.array(
                        [  [ (k[0]-i[0])*(k[0]-i[0]) , (k[0]-i[0])*(k[1]-i[1]) , (k[0]-i[0])*(k[2]-i[2]) ]  , 
                        [ (k[1]-i[1])*(k[0]-i[0]) , (k[1]-i[1])*(k[1]-i[1]) , (k[1]-i[1])*(k[2]-i[2]) ]  ,
                        [ (k[2]-i[2])*(k[0]-i[0]) , (k[2]-i[2])*(k[1]-i[1]) , (k[2]-i[2])*(k[2]-i[2]) ]  ]
                        ) * ( -self.gamma/(dist_ik**2) )
                    )
            return -diag

        #To Verify that the hessian is same as ANM (Part 2/2)
        #print(Index1,Index2,dist_ij,'\n',self.get_hessian()[Index1*3:Index1*3+3,Index2*3:Index2*3+3],"\n####\n")
        #return self.get_hessian()[Index1*3:Index1*3+3,Index2*3:Index2*3+3]


    '''Calculate Functions'''
    #Discontinued (But kept for testing)
    '''
    def calculate_hessian(self):
        """Build the Hessian Matrix of the ANM model.

        This is an essential step for hdANM
        
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
        
        #Class Change
        self.hessian=hessian
        return True
    '''
    
    def calculate_decomposition(self, include_mass=True):
        """Decompose the Hessian Matrix of the hdANM model.

        Args:
            include_mass (bool): Amino Acid Residue mass/Atomic mass should be (True) or shouldn't be (False) included for the decomposition. Defaults to True
        
        Note:
            - Eigen values and Eigen Vectors are calculated. use hdANM().get_eigenvalues() and hdANM().get_eigenvectors() to obtain them.
            - Currently only molecular weight is included in case of Amino Acid Residue(coarse grained) mass.
            - Mass of the Amino Acid Residues/Atoms can be found in /packman/constants/Constants.py 
            ie...
            from packman.constants import amino_acid_molecular_weight
            from packman.constants import atomic_weigh
        """
        if(include_mass):
            self.eigen_values, self.eigen_vectors= scipy_eig(self.domain_hessian,self.domain_mass_matrix)

            idx                = self.eigen_values.argsort()
            self.eigen_values  = self.eigen_values[idx]
            self.eigen_vectors = self.eigen_vectors[:,idx]
        else:
            self.eigen_values, self.eigen_vectors = numpy.linalg.eigh(self.domain_hessian)
        return True
    
    
    def calculate_hessian(self,mass_type='unit'):
        """Build the Hessian Matrix of the hdANM model.

        This is the most essential step for hdANM. It picks up blocks from ANM hessian matrix and puts it in the format described in the paper.

        Args:
            mass_type (unit/atom/residue): Whether to use unit (1), atomic weights or residue mass for the mass matrix.
        
        Notes:
            * Possible argument removal for the decomposition function
            * use hdANM().domain_hessian, hdANM().domain_mass_matrix and hdANM().domain_info to access output of this function
            * Mass of the Amino Acid Residues/Atoms can be found in /packman/constants/Constants.py
            ie...
            from packman.constants import amino_acid_molecular_weight
            from packman.constants import atomic_weigh
        """
        all_mass_types=['unit','atom','residue']
        #Error Handling
        if(mass_type not in all_mass_types):
            print('Provide correct \'mass_type\' parameter. See help for more details')


        for i in self.atoms:
            for j in self.HNGinfo:
                HNGrange = self.HNGinfo[j]
                IDs  = j.split('_')
                if( HNGrange[0] <= i.get_parent().get_id() <= HNGrange[1] and IDs[1] == i.get_parent().get_parent().get_id() ):
                    i.get_parent().set_domain_id( IDs[-1].strip() )

        #Domain groups stores the atom index(of self.atoms) with Domain ID as a key (D5: [1,2] domain 5 with atom index 1 and 2 )
        DomainGroups={}
        for numi,i in enumerate(self.atoms):
            try:
                DomainGroups[i.get_domain_id()].append(numi)
            except:
                DomainGroups[i.get_domain_id()]=[]
                DomainGroups[i.get_domain_id()].append(numi)
        
        #HHH
        HingeAtoms=[]
        Domains=[]   #Used in HHD
        for i in DomainGroups.keys():
            if(i[0]=='H'):
                HingeAtoms.extend(DomainGroups[i])
            if(i[0]=='D'):
                Domains.append(i)
        HHH=numpy.zeros((len(HingeAtoms)*3,len(HingeAtoms)*3))
        for numi,i in enumerate(HingeAtoms):
            for numj,j in enumerate(HingeAtoms):
                HHH[numi*3:numi*3+3,numj*3:numj*3+3]=self.get_hessian_block(i,j)
        
        
        #HHD (Scheduled to be paralalysed)
        DomainInfo={}
        HHD=numpy.zeros((len(HingeAtoms)*3,len(Domains)*6))
        for numd,d in enumerate(Domains):
            COM=numpy.average([self.coords[i] for i in DomainGroups[d]] , 0)
            DomainInfo[d]=(numd,COM)
            for numi,i in enumerate(HingeAtoms):
                HHDijd_sum=HHDijd=numpy.zeros((3,6))
                for j in DomainGroups[d]:
                    HHDijd=numpy.zeros((3,6))
                    HHDijd[0:3,0:3]=self.get_hessian_block(i,j)

                    HHDijd[0,3:]= [  -HHDijd[0][1]*(self.coords[j][2]-COM[2]) + HHDijd[0][2]*(self.coords[j][1]-COM[1])  ,  HHDijd[0][0]*(self.coords[j][2]-COM[2]) - HHDijd[0][2]*(self.coords[j][0]-COM[0])  ,  -HHDijd[0][0]*(self.coords[j][1]-COM[1]) + HHDijd[0][1]*(self.coords[j][0]-COM[0])]
                    HHDijd[1,3:]= [  -HHDijd[1][1]*(self.coords[j][2]-COM[2]) + HHDijd[1][2]*(self.coords[j][1]-COM[1])  ,  HHDijd[1][0]*(self.coords[j][2]-COM[2]) - HHDijd[1][2]*(self.coords[j][0]-COM[0])  ,  -HHDijd[1][0]*(self.coords[j][1]-COM[1]) + HHDijd[1][1]*(self.coords[j][0]-COM[0])]
                    HHDijd[2,3:]= [  -HHDijd[2][1]*(self.coords[j][2]-COM[2]) + HHDijd[2][2]*(self.coords[j][1]-COM[1])  ,  HHDijd[2][0]*(self.coords[j][2]-COM[2]) - HHDijd[2][2]*(self.coords[j][0]-COM[0])  ,  -HHDijd[2][0]*(self.coords[j][1]-COM[1]) + HHDijd[2][1]*(self.coords[j][0]-COM[0])]
                    HHDijd_sum=numpy.add(HHDijd_sum,HHDijd)
                HHD[numi*3:(numi*3)+3,numd*6:((numd*6)+6)]=HHDijd_sum
        

        #HDH
        HDH=HHD.T

        #HDD
        HDD=numpy.zeros((6*len(Domains),6*len(Domains)))
        for numd,Dd in enumerate(Domains):
            for nume,De in enumerate(Domains):
                HDdDe=numpy.zeros((6,6))
                Dd_COM=numpy.average([self.coords[i] for i in DomainGroups[Dd]] , 0)
                De_COM=numpy.average([self.coords[i] for i in DomainGroups[De]] , 0)
                for i in DomainGroups[Dd]:
                    for j in DomainGroups[De]:
                        HDDij=self.get_hessian_block(i,j)
                        HDdDe_iter=numpy.zeros((6,6))
                        #F_Delta
                        HDdDe_iter[0:3,0:3]=HDDij
                        #F_Phi
                        HDdDe_iter[0,3:]= [  -HDDij[0][1]*(self.coords[j][2]-De_COM[2]) + HDDij[0][2]*(self.coords[j][1]-De_COM[1])  ,  HDDij[0][0]*(self.coords[j][2]-De_COM[2]) - HDDij[0][2]*(self.coords[j][0]-De_COM[0])  ,  -HDDij[0][0]*(self.coords[j][1]-De_COM[1]) + HDDij[0][1]*(self.coords[j][0]-De_COM[0])]
                        HDdDe_iter[1,3:]= [  -HDDij[1][1]*(self.coords[j][2]-De_COM[2]) + HDDij[1][2]*(self.coords[j][1]-De_COM[1])  ,  HDDij[1][0]*(self.coords[j][2]-De_COM[2]) - HDDij[1][2]*(self.coords[j][0]-De_COM[0])  ,  -HDDij[1][0]*(self.coords[j][1]-De_COM[1]) + HDDij[1][1]*(self.coords[j][0]-De_COM[0])]
                        HDdDe_iter[2,3:]= [  -HDDij[2][1]*(self.coords[j][2]-De_COM[2]) + HDDij[2][2]*(self.coords[j][1]-De_COM[1])  ,  HDDij[2][0]*(self.coords[j][2]-De_COM[2]) - HDDij[2][2]*(self.coords[j][0]-De_COM[0])  ,  -HDDij[2][0]*(self.coords[j][1]-De_COM[1]) + HDDij[2][1]*(self.coords[j][0]-De_COM[0])]
                        #Torque_Delta
                        HDdDe_iter[3,0:3]= [  -HDDij[1][0]*(self.coords[i][2]-Dd_COM[2]) + HDDij[2][0]*(self.coords[i][1]-Dd_COM[1])  ,  -HDDij[1][1]*(self.coords[i][2]-Dd_COM[2]) + HDDij[2][1]*(self.coords[i][1]-Dd_COM[1])  ,  -HDDij[1][2]*(self.coords[i][2]-Dd_COM[2]) + HDDij[2][2]*(self.coords[i][1]-Dd_COM[1])]
                        HDdDe_iter[4,0:3]= [   HDDij[0][0]*(self.coords[i][2]-Dd_COM[2]) - HDDij[2][0]*(self.coords[i][0]-Dd_COM[0])  ,   HDDij[0][1]*(self.coords[i][2]-Dd_COM[2]) - HDDij[2][1]*(self.coords[i][0]-Dd_COM[0])  ,   HDDij[0][2]*(self.coords[i][2]-Dd_COM[2]) - HDDij[2][2]*(self.coords[i][0]-Dd_COM[0])]
                        HDdDe_iter[5,0:3]= [  -HDDij[0][0]*(self.coords[i][1]-Dd_COM[1]) + HDDij[1][0]*(self.coords[i][0]-Dd_COM[0])  ,  -HDDij[0][1]*(self.coords[i][1]-Dd_COM[1]) + HDDij[1][1]*(self.coords[i][0]-Dd_COM[0])  ,  -HDDij[0][2]*(self.coords[i][1]-Dd_COM[1]) + HDDij[1][2]*(self.coords[i][0]-Dd_COM[0])]
                        #Torque_Phi
                        HDdDe_iter[3,3:]= [   HDDij[1][1]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][2]-De_COM[2]) - HDDij[1][2]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][1]-De_COM[1]) - HDDij[2][1]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][2]-De_COM[2]) + HDDij[2][2]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][1]-De_COM[1])  ,  -HDDij[1][0]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][2]-De_COM[2]) + HDDij[1][2]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][0]-De_COM[0]) + HDDij[2][0]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][2]-De_COM[2]) - HDDij[2][2]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][0]-De_COM[0])  ,   HDDij[1][0]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][1]-De_COM[1]) - HDDij[1][1]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][0]-De_COM[0]) - HDDij[2][0]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][1]-De_COM[1]) + HDDij[2][1]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][0]-De_COM[0])  ]
                        HDdDe_iter[4,3:]=[  -HDDij[0][1]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][2]-De_COM[2]) + HDDij[0][2]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][1]-De_COM[1]) + HDDij[2][1]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][2]-De_COM[2]) - HDDij[2][2]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][1]-De_COM[1])  ,   HDDij[0][0]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][2]-De_COM[2]) - HDDij[0][2]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][0]-De_COM[0]) - HDDij[2][0]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][2]-De_COM[2]) + HDDij[2][2]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][0]-De_COM[0])  ,  -HDDij[0][0]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][1]-De_COM[1]) + HDDij[0][1]*(self.coords[i][2]-Dd_COM[2])*(self.coords[j][0]-De_COM[0]) + HDDij[2][0]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][1]-De_COM[1]) - HDDij[2][1]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][0]-De_COM[0])  ]
                        HDdDe_iter[5,3:]=[   HDDij[0][1]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][2]-De_COM[2]) - HDDij[0][2]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][1]-De_COM[1]) - HDDij[1][1]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][2]-De_COM[2]) + HDDij[1][2]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][1]-De_COM[1])  ,  -HDDij[0][0]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][2]-De_COM[2]) + HDDij[0][2]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][0]-De_COM[0]) + HDDij[1][0]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][2]-De_COM[2]) - HDDij[1][2]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][0]-De_COM[0])  ,   HDDij[0][0]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][1]-De_COM[1]) - HDDij[0][1]*(self.coords[i][1]-Dd_COM[1])*(self.coords[j][0]-De_COM[0]) - HDDij[1][0]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][1]-De_COM[1]) + HDDij[1][1]*(self.coords[i][0]-Dd_COM[0])*(self.coords[j][0]-De_COM[0])  ]
                        HDdDe=numpy.add(HDdDe,HDdDe_iter)
                HDD[numd*6:(numd*6)+6 , nume*6:(nume*6)+6]=HDdDe
        
        
        #Reconstruction
        new_dim=6*len(Domains)+3*len(HingeAtoms)
        H_new=numpy.zeros((new_dim,new_dim))
        H_new[0:6*len(Domains), 0:6*len(Domains)]=HDD
        H_new[0:6*len(Domains), 6*len(Domains):]=HDH
        H_new[6*len(Domains): ,0:6*len(Domains)]=HHD
        H_new[6*len(Domains):,6*len(Domains):]=HHH
        
        #Mass Matrix
        M=numpy.zeros((H_new.shape))

        #MHH
        if(mass_type == 'unit'):
            M[6*len(Domains):,6*len(Domains):]=numpy.identity(len(HingeAtoms)*3)
        elif(mass_type == 'atom'):
            try:
                lst= [atomic_weight[self.atoms[i].get_element()] for i in HingeAtoms]
                lst= list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in lst))
                numpy.fill_diagonal( M[6*len(Domains):,6*len(Domains):] , lst )
            except:
                print('Unknown atom or the atomic weight not available in packman.constants')
        elif(mass_type == 'residue'):
            try:
                lst= [amino_acid_molecular_weight[self.atoms[i].get_parent().get_name()] for i in HingeAtoms]
                lst= list(itertools.chain.from_iterable(itertools.repeat(x, 3) for x in lst))
                numpy.fill_diagonal( M[6*len(Domains):,6*len(Domains):] , lst )
            except:
                print('Unknown Amino Acid encountered or molecular weight not available in packman.constants')

        #MDD (Scheduled to be paralalysed)
        for numd,d in enumerate(Domains):
            #F-Delta
            domain_block=numpy.zeros((6,6))
            if(mass_type == 'unit'):
                numpy.fill_diagonal( domain_block[0:3,0:3] , len(DomainGroups[d]) )
            elif(mass_type == 'atom'):
                DomainMass = numpy.sum( [atomic_weight[self.atoms[i].get_element()] for i in DomainGroups[d]] )
                numpy.fill_diagonal( domain_block[0:3,0:3] , DomainMass )
            elif(mass_type == 'residue'):
                DomainMass = numpy.sum( [amino_acid_molecular_weight[self.atoms[i].get_parent().get_name()] for i in DomainGroups[d]] )
                numpy.fill_diagonal( domain_block[0:3,0:3] , DomainMass )

            #Torque_Phi
            COM=DomainInfo[d][1]
            for numi,i in enumerate(DomainGroups[d]):
                MDdDd_local=numpy.array(
                    [   
                        [  (self.coords[i][2]-COM[2])**2+(self.coords[i][1]-COM[1])**2  ,  -(self.coords[i][1]-COM[1])*(self.coords[i][0]-COM[0])        ,  -(self.coords[i][0]-COM[0])*(self.coords[i][2]-COM[2])        ],
                        [ -(self.coords[i][1]-COM[1])*(self.coords[i][0]-COM[0])        ,   (self.coords[i][2]-COM[2])**2+(self.coords[i][0]-COM[0])**2  ,  -(self.coords[i][1]-COM[1])*(self.coords[i][2]-COM[2])        ],
                        [ -(self.coords[i][0]-COM[0])*(self.coords[i][2]-COM[2])        ,  -(self.coords[i][1]-COM[1])*(self.coords[i][2]-COM[2])        ,   (self.coords[i][0]-COM[0])**2+(self.coords[i][1]-COM[1])**2  ],
                    ])
                
                #Maybe remove if conditions later, create separate functions for improving the speed
                if(mass_type == 'unit'):
                    domain_block[3:,3:]=numpy.add(domain_block[3:,3:],MDdDd_local)
                elif(mass_type == 'atom'):
                    MDdDd_local = MDdDd_local * atomic_weight[self.atoms[i].get_element()]
                    domain_block[3:,3:]=numpy.add(domain_block[3:,3:],MDdDd_local)
                elif(mass_type == 'residue'):
                    MDdDd_local = MDdDd_local * amino_acid_molecular_weight[self.atoms[i].get_parent().get_name()]
                    domain_block[3:,3:]=numpy.add(domain_block[3:,3:],MDdDd_local)

            M[numd*6:(numd*6)+6 , numd*6:(numd*6)+6]=domain_block

        self.domain_hessian, self.domain_mass_matrix, self.domain_info = H_new, M, DomainInfo #DomainInfo is saved to check the sequence in which H_new is formed
        return True
    
    def calculate_RT_eigen_vectors(self):
        """Calculate the reverse transformed vectors from the hdANM eigenvectors. 
    
        Reverse transformed means that the hdANM eigenvector of dimension: 6D+3H x 6D+3H  (D: Number of domains; H: Number of hinge atoms) are converted to 3N x 6D+3H (N: Number of atoms)

        Note:
            - It was refered as 'exploded vector' utill version 1.3.3
        """
        x0 = numpy.array([i.get_location() for i in self.atoms])
        d0 = [i.get_domain_id() for i in self.atoms]
        
        exploded_vectors=numpy.zeros((len(x0)*3, len(self.eigen_values) ))

        for mode_number in range(0,self.eigen_vectors.shape[1],1):
            new_column=[]
            HingeResidue=0
            for numi,i in enumerate(x0):
                if(d0[numi][0]=='D'):
                    D_delta_phi= self.eigen_vectors[:,mode_number][self.domain_info[d0[numi]][0]*6 : (self.domain_info[d0[numi]][0]*6)+6]
                    D_COM=self.domain_info[d0[numi]][1]
                    new_x= ( D_delta_phi[0] + (D_delta_phi[4]*(i[2]-D_COM[2])) - (D_delta_phi[5]*(i[1]-D_COM[1])) )
                    new_y= ( D_delta_phi[1] + (D_delta_phi[3]*(i[2]-D_COM[2])) - (D_delta_phi[5]*(i[0]-D_COM[0])) )
                    new_z= ( D_delta_phi[2] + (D_delta_phi[3]*(i[1]-D_COM[1])) - (D_delta_phi[4]*(i[0]-D_COM[0])) )
                if(d0[numi][0]=='H'):
                    delta=self.eigen_vectors[:,mode_number][6*len(self.domain_info):][HingeResidue*3:(HingeResidue*3)+3]
                    new_x=delta[0]
                    new_y=delta[1]
                    new_z=delta[2]
                    HingeResidue+=1
                new_column.extend([new_x,new_y,new_z])
            exploded_vectors[:,mode_number]=numpy.array(new_column)
        
        self.RT_eigen_vectors = exploded_vectors
        return True

    def calculate_fluctuations(self):
        """Calculate the Fluctuations of the hd-ANM model.

        The fluctualtions/ theoretical b-factors are calculated using this method.
        
        Note:
            - Fluctuations are calculated. use hdANM().get_fluctuations() to obtain the fluctuations.
            - Endmode needs to be put in the code if and when required.
            - hdANM().fluctuations stores the output of this function.
        """
        x0 = numpy.array([i.get_location() for i in self.atoms])

        if(self.RT_eigen_vectors == None):
            self.calculate_RT_eigen_vectors()
        
        exploded_vectors = self.RT_eigen_vectors
        exploded_vectors = exploded_vectors.T

        mode_bfactors=[]
        for numi,i in enumerate(self.eigen_values[6:]):
            evec_row=exploded_vectors[numi+6]
            mode_bfactors.append([ float(evec_row[j]**2 + evec_row[j+1]**2 + evec_row[j+2]**2)/i for j in range(0,len(x0)*3,3)])
        
        mode_bfactors=numpy.array(mode_bfactors)
        self.fluctuations = [numpy.sum(i) for i in mode_bfactors.T]
        return True
    
    def calculate_hessian_pseudoinverse(self, n_modes="all"):
        """Calculate the hessian pseudoinverse for the hd-ANM modes

        Pseudoinverse is calculated using this methd. Please note that first 6 rigid modes are eliminated in this calculation.

        Note:
            - Can be called on demand or can be called from get method automatically (Needs thinking)
        """
        
            
        if(self.RT_eigen_vectors == None):
            self.calculate_RT_eigen_vectors()
        EVec=self.get_RT_eigen_vectors().T
        if(n_modes=="all"):
            self.hessian_pseudoinverse = numpy.matmul(  numpy.matmul( EVec[6:].transpose() , numpy.diag(1/self.get_eigenvalues()[6:]) )  , EVec[6:]  )
        else:
            self.hessian_pseudoinverse = numpy.matmul(  numpy.matmul( EVec[6:6+n_modes].transpose() , numpy.diag(1/self.get_eigenvalues()[6:6+n_modes]) )  , EVec[6:6+n_modes]  )
    
    def calculate_cross_correlation( self, n_modes = "all" ):
        """Calculate the cross correlation matrix for the hdANM modes.

        Crosscorrelation matrix is generated for all modes by default. Please change the n_modes parameter to restrict modes.

        Args:
            -   n_modes (int): Number of modes that need to be considered to calculate the cross correlation matrix. 
        """
        self.calculate_hessian_pseudoinverse(n_modes=n_modes)
        
        self.crosscorrelation_matrix = numpy.zeros( (len(self.atoms), len(self.atoms)) )
        for i in range(0,len(self.atoms)):
            for j in range(0,len(self.atoms)):
                trace_H_inv_ij = numpy.trace( self.get_hessian_pseudoinverse()[i*3:i*3+3,j*3:j*3+3] )
                trace_H_inv_ii = numpy.trace( self.get_hessian_pseudoinverse()[i*3:i*3+3,i*3:i*3+3] )
                trace_H_inv_jj = numpy.trace( self.get_hessian_pseudoinverse()[j*3:j*3+3,j*3:j*3+3] )

                self.crosscorrelation_matrix[i][j] = trace_H_inv_ij / numpy.sqrt( trace_H_inv_ii*trace_H_inv_jj )

    def calculate_movie(self, mode_number, scale=1.5, n=20, extrapolation="curvilinear", ftype='cif', ca_to_aa=False):
        """This function generates the dynamic 3D projection of the normal modes obtained using hd-ANM. The 3D projection can be linearly extrapolated or curvilinearly extrapolated depending on the choices. The first frame is the original structure and the projection progresses in positive (+) direction, returning to original structure and then in negative direction (-) again returning to the original structure.

        Args:
            mode_number (int)                   : Mode number. (first non-rigid mode is 6th)
            scale (float)                       : Multiplier; extent to which mode will be extrapolated.                 Defaults to 1.5
            n (int)                             : Number of frames in output (should be =>8 and ideally multiple of 4)   Defaults to 20
            extrapolation (linear/curvilinear)  : Extrapolation method                                                   Defaults to "curvilinear"
            ftype (string)                      : Extension of the output file (.cif / .pdb)
            ca_to_aa (boolean)                  : If only single atom type is used (often C-alpha atom only), enabling extrapolates the C-alpha motion to all the atoms. (Default: False)

        Note:
            - Scale and n parameters should be redesigned.
            - direction is the variable which be allow user to explore only positive or only negative direction of the modes.

        Returns:
            True if successful; false otherwise.
        """

        if(ca_to_aa):
            if(len(list(set([i.get_name() for i in self.atoms]))) > 1 ): logging.warning('Multiple atom types were detected, but the "ca_to_aa" option is still "True." If the model is not coarse-grained, consider making it "False." See the function description for more details.')
        else:
          if(len(list(set([i.get_name() for i in self.atoms]))) == 1 ): logging.info('A single type of atom is detected. Try enabling the "ca_to_aa" parameter. See the function description for more details.')

        x0=numpy.array([i.get_location() for i in self.atoms])
        d0=[i.get_domain_id() for i in self.atoms]
        new_coords=[]

        if(extrapolation=="linear"):
            movement = [ numpy.sin(k*(1.0/float(n))*2*numpy.pi) for k in range(0,n+1,1) ]
            ModelsOfTheProtein = []
            for j in movement:
                HingeResidue=0
                AtomsOfTheFrame = {}
                non_model_atoms = 0 #Only when ca_to_aa option is enabled.
                for numi,i in enumerate(x0):
                    if(d0[numi][0]=='D'):
                        D_delta_phi= self.eigen_vectors[:,mode_number][self.domain_info[d0[numi]][0]*6 : (self.domain_info[d0[numi]][0]*6)+6]
                        D_COM=self.domain_info[d0[numi]][1]
                        new_x=i[0]+scale*j* ( D_delta_phi[0] + (D_delta_phi[4]*(i[2]-D_COM[2])) - (D_delta_phi[5]*(i[1]-D_COM[1])) )
                        new_y=i[1]+scale*j* ( D_delta_phi[1] + (D_delta_phi[3]*(i[2]-D_COM[2])) - (D_delta_phi[5]*(i[0]-D_COM[0])) )
                        new_z=i[2]+scale*j* ( D_delta_phi[2] + (D_delta_phi[3]*(i[1]-D_COM[1])) - (D_delta_phi[4]*(i[0]-D_COM[0])) )
                        if(ca_to_aa):
                            for x in self.atoms[numi].get_parent().get_atoms():
                                temp_new_x=x.get_location()[0]+scale*j* ( D_delta_phi[0] + (D_delta_phi[4]*(x.get_location()[2]-D_COM[2])) - (D_delta_phi[5]*(x.get_location()[1]-D_COM[1])) )
                                temp_new_y=x.get_location()[1]+scale*j* ( D_delta_phi[1] + (D_delta_phi[3]*(x.get_location()[2]-D_COM[2])) - (D_delta_phi[5]*(x.get_location()[0]-D_COM[0])) )
                                temp_new_z=x.get_location()[2]+scale*j* ( D_delta_phi[2] + (D_delta_phi[3]*(x.get_location()[1]-D_COM[1])) - (D_delta_phi[4]*(x.get_location()[0]-D_COM[0])) )

                                non_model_atoms+=1
                                AtomsOfTheFrame[len(x0)+non_model_atoms] = Atom(x.get_id() , x.get_name(), numpy.array([temp_new_x,temp_new_y,temp_new_z]), x.get_occupancy(), x.get_bfactor(), x.get_element(), x.get_charge(), x.get_parent() )
                                
                    if(d0[numi][0]=='H'):
                        delta=self.eigen_vectors[:,mode_number][6*len(self.domain_info):][HingeResidue*3:(HingeResidue*3)+3]
                        new_x=i[0]+scale*j*delta[0]
                        new_y=i[1]+scale*j*delta[1]
                        new_z=i[2]+scale*j*delta[2]
                        HingeResidue+=1
                        if(ca_to_aa):
                            for x in self.atoms[numi].get_parent().get_atoms():
                                temp_new_x=x.get_location()[0]+scale*j*delta[0]
                                temp_new_y=x.get_location()[1]+scale*j*delta[1]
                                temp_new_z=x.get_location()[2]+scale*j*delta[2]

                                non_model_atoms+=1
                                AtomsOfTheFrame[len(x0)+non_model_atoms] = Atom(x.get_id() , x.get_name(), numpy.array([temp_new_x,temp_new_y,temp_new_z]), x.get_occupancy(), x.get_bfactor(), x.get_element(), x.get_charge(), x.get_parent() )
                                
                    #All the atoms in the models (CA or all atom); Atoms from ca_to_aa option are not included.
                    new_x , new_y, new_z = new_x.real , new_y.real , new_z.real
                    new_coords.append([new_x,new_y,new_z])
                    currentatm = Atom(self.atoms[numi].get_id() , self.atoms[numi].get_name(), numpy.array([new_x,new_y,new_z]), self.atoms[numi].get_occupancy(), self.atoms[numi].get_bfactor(), self.atoms[numi].get_element(),self.atoms[numi].get_charge(), self.atoms[numi].get_parent() )
                    AtomsOfTheFrame[numi] = currentatm
                
                ModelsOfTheProtein.append( Model(j, AtomsOfTheFrame, None, None, None, None) )
            
            Annotations = self.atoms[0].get_parent().get_parent().get_parent().get_parent().get_data()
            prot = Protein(str(mode_number), ModelsOfTheProtein)
            prot.set_data(Annotations)
            prot.write_structure( str(mode_number)+'.'+ftype )
            

        elif(extrapolation=="curvilinear"):
            movement = [ numpy.sin(k*(1.0/float(n))*2*numpy.pi) for k in range(0,n+1,1) ]
            ModelsOfTheProtein = []
            for j in movement:
                HingeResidue=0
                AtomsOfTheFrame = {}
                non_model_atoms = 0 #Only when ca_to_aa option is enabled.
                for numi,i in enumerate(x0):
                    if(d0[numi][0]=='D'):
                        D_delta_phi= self.eigen_vectors[:,mode_number][self.domain_info[d0[numi]][0]*6 : (self.domain_info[d0[numi]][0]*6)+6]
                        #D_mu= D_delta_phi[3:]
                        Q_D_n=numpy.linalg.norm(D_delta_phi[3:])
                        D_mu=D_delta_phi[3:] / Q_D_n
                        D_COM=self.domain_info[d0[numi]][1]
                        #D_delta_phi[0] is questionable ; D_delta_phi[3:] is also questionable
                        new_x= D_COM[0]+ scale*j*D_delta_phi[0] \
                        + ( numpy.cos(scale*j*Q_D_n)+ D_mu[0]**2 * (1-numpy.cos(scale*j*Q_D_n)) )             * (i[0] -D_COM[0]) \
                        + ( D_mu[0]*D_mu[1]*(1-numpy.cos(scale*j*Q_D_n)) - D_mu[2]*numpy.sin(scale*j*Q_D_n) ) * (i[1] -D_COM[1]) \
                        + ( D_mu[0]*D_mu[2]*(1-numpy.cos(scale*j*Q_D_n)) + D_mu[1]*numpy.sin(scale*j*Q_D_n) ) * (i[2] -D_COM[2])
                        
                        new_y= D_COM[1]+ scale*j*D_delta_phi[1] \
                        + ( D_mu[0]*D_mu[1]*(1-numpy.cos(scale*j*Q_D_n)) + D_mu[2]*numpy.sin(scale*j*Q_D_n) ) * (i[0] -D_COM[0]) \
                        + ( numpy.cos(scale*j*Q_D_n)+D_mu[1]**2 * (1-numpy.cos(scale*j*Q_D_n)) )              * (i[1] -D_COM[1]) \
                        + ( D_mu[1]*D_mu[2]*(1-numpy.cos(scale*j*Q_D_n)) - D_mu[0]*numpy.sin(scale*j*Q_D_n) ) * (i[2] -D_COM[2])

                        new_z= D_COM[2]+ scale*j*D_delta_phi[2] \
                        + ( D_mu[0]*D_mu[2]*(1-numpy.cos(scale*j*Q_D_n)) - D_mu[1]*numpy.sin(scale*j*Q_D_n) ) * (i[0] -D_COM[0]) \
                        + ( D_mu[1]*D_mu[2]*(1-numpy.cos(scale*j*Q_D_n)) + D_mu[0]*numpy.sin(scale*j*Q_D_n) ) * (i[1] -D_COM[1]) \
                        + ( numpy.cos(scale*j*Q_D_n)+ D_mu[2]**2 * (1-numpy.cos(scale*j*Q_D_n)) )             * (i[2] -D_COM[2])

                        if(ca_to_aa):
                            for x in self.atoms[numi].get_parent().get_atoms():
                                temp_new_x= D_COM[0]+ scale*j*D_delta_phi[0] \
                                + ( numpy.cos(scale*j*Q_D_n)+ D_mu[0]**2 * (1-numpy.cos(scale*j*Q_D_n)) )             * (x.get_location()[0] -D_COM[0]) \
                                + ( D_mu[0]*D_mu[1]*(1-numpy.cos(scale*j*Q_D_n)) - D_mu[2]*numpy.sin(scale*j*Q_D_n) ) * (x.get_location()[1] -D_COM[1]) \
                                + ( D_mu[0]*D_mu[2]*(1-numpy.cos(scale*j*Q_D_n)) + D_mu[1]*numpy.sin(scale*j*Q_D_n) ) * (x.get_location()[2] -D_COM[2])
                                
                                temp_new_y= D_COM[1]+ scale*j*D_delta_phi[1] \
                                + ( D_mu[0]*D_mu[1]*(1-numpy.cos(scale*j*Q_D_n)) + D_mu[2]*numpy.sin(scale*j*Q_D_n) ) * (x.get_location()[0] -D_COM[0]) \
                                + ( numpy.cos(scale*j*Q_D_n)+D_mu[1]**2 * (1-numpy.cos(scale*j*Q_D_n)) )              * (x.get_location()[1] -D_COM[1]) \
                                + ( D_mu[1]*D_mu[2]*(1-numpy.cos(scale*j*Q_D_n)) - D_mu[0]*numpy.sin(scale*j*Q_D_n) ) * (x.get_location()[2] -D_COM[2])

                                temp_new_z= D_COM[2]+ scale*j*D_delta_phi[2] \
                                + ( D_mu[0]*D_mu[2]*(1-numpy.cos(scale*j*Q_D_n)) - D_mu[1]*numpy.sin(scale*j*Q_D_n) ) * (x.get_location()[0] -D_COM[0]) \
                                + ( D_mu[1]*D_mu[2]*(1-numpy.cos(scale*j*Q_D_n)) + D_mu[0]*numpy.sin(scale*j*Q_D_n) ) * (x.get_location()[1] -D_COM[1]) \
                                + ( numpy.cos(scale*j*Q_D_n)+ D_mu[2]**2 * (1-numpy.cos(scale*j*Q_D_n)) )             * (x.get_location()[2] -D_COM[2])

                                non_model_atoms+=1
                                AtomsOfTheFrame[len(x0)+non_model_atoms] = Atom(x.get_id() , x.get_name(), numpy.array([temp_new_x,temp_new_y,temp_new_z]), x.get_occupancy(), x.get_bfactor(), x.get_element(), x.get_charge(), x.get_parent() )


                    if(d0[numi][0]=='H'):
                        delta=self.eigen_vectors[:,mode_number][6*len(self.domain_info):][HingeResidue*3:(HingeResidue*3)+3]
                        new_x=i[0]+scale*j*delta[0]
                        new_y=i[1]+scale*j*delta[1]
                        new_z=i[2]+scale*j*delta[2]
                        HingeResidue+=1

                        if(ca_to_aa):
                            for x in self.atoms[numi].get_parent().get_atoms():
                                temp_new_x=x.get_location()[0]+scale*j*delta[0]
                                temp_new_y=x.get_location()[1]+scale*j*delta[1]
                                temp_new_z=x.get_location()[2]+scale*j*delta[2]
                                
                                non_model_atoms+=1
                                AtomsOfTheFrame[len(x0)+non_model_atoms] = Atom(x.get_id() , x.get_name(), numpy.array([temp_new_x,temp_new_y,temp_new_z]), x.get_occupancy(), x.get_bfactor(), x.get_element(), x.get_charge(), x.get_parent() )

                    new_x , new_y, new_z = new_x.real , new_y.real , new_z.real
                    new_coords.append([new_x,new_y,new_z])
                    currentatm = Atom(self.atoms[numi].get_id() , self.atoms[numi].get_name(), numpy.array([new_x,new_y,new_z]), self.atoms[numi].get_occupancy(), self.atoms[numi].get_bfactor(), self.atoms[numi].get_element(),self.atoms[numi].get_charge(), self.atoms[numi].get_parent() )
                    AtomsOfTheFrame[numi] = currentatm
                    
                ModelsOfTheProtein.append( Model(j, AtomsOfTheFrame, None, None, None, None) )

            Annotations = self.atoms[0].get_parent().get_parent().get_parent().get_parent().get_data()
            prot = Protein(str(mode_number), ModelsOfTheProtein)
            prot.set_data(Annotations)
            prot.write_structure( str(mode_number)+'.'+ftype )
            
        return True