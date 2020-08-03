'''
Copyright Pranav Khade

Author: Pranav Khade(pranavk@iastate.edu)
'''

from packman import molecule
from packman.constants import amino_acid_molecular_weight
from packman.constants import atomic_weight

import numpy
import itertools

from scipy.linalg import eig as scipy_eig
'''
##################################################################################################
#                                            RD-ANM                                              #
##################################################################################################
'''

class RDANM:
    def __init__(self, atoms , gamma=1.0, dr=15.0, power=0, pf=None, HNGinfo=None):
        """
        Author: Pranav Khade
        """
        self.gamma   = gamma
        self.dr      = dr
        self.power   = power
        self.pf      = pf
        self.atoms   = atoms
        self.HNGinfo = HNGinfo
        self.coords  = numpy.array([i.get_location() for i in atoms])
        if self.pf != None and self.pf <= 0:
            raise Exception("pf value cannot be zero or negative")
        if self.gamma <= 0:
            raise Exception("gamma value cannot be zero or negative")
        if self.dr <= 0:
            raise Exception("distance cutoff value cannot be zero or negative")
        
        self.calculate_hessian()


    '''Get Functions'''
    def get_hessian(self):
        return self.hessian
    
    def get_eigenvalues(self):
        return self.eigen_values
    
    def get_eigenvectors(self):
        return self.eigen_vectors
    
    def get_fluctuations(self):
        return self.fluctuations
    
    def get_hessian_block(self,Index1,Index2):
        #To demonstrate how index selection works
        #toy_example=numpy.reshape( numpy.arange(0,81,1), (9,9))
        #print(toy_example,'\n',toy_example[Index1*3:Index1*3+3,Index2*3:Index2*3+3])
        return self.get_hessian()[Index1*3:Index1*3+3,Index2*3:Index2*3+3]


    '''Calculate Functions'''
    def calculate_hessian(self):
        """
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
    
    def calculate_decomposition(self, include_mass=True):
        """
        """
        if(include_mass):
            self.eigen_values, self.eigen_vectors= scipy_eig(self.domain_hessian,self.domain_mass_matrix)

            idx                = self.eigen_values.argsort()
            self.eigen_values  = self.eigen_values[idx]
            self.eigen_vectors = self.eigen_vectors[:,idx]
        else:
            self.eigen_values,self.eigen_vectors = numpy.linalg.eigh(self.domain_hessian)
        return True
    
    
    def calculate_coarse_grained_hessian(self,mass_type='unit'):
        """
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
                    i.get_parent().set_domain_id( IDs[-1] )
        
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
        
        
        #HHD
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

        #MDD
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

        self.domain_hessian, self.domain_mass_matrix, self.domain_info= H_new, M, DomainInfo #DomainInfo is saved to check the sequence in which H_new is formed
        return True
    

    def calculate_fluctuations(self):
        '''
        '''
        x0=numpy.array([i.get_location() for i in self.atoms])
        d0=[i.get_domain_id() for i in self.atoms]
        
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
        
        numpy.savetxt('Explooded.txt',exploded_vectors)
        exploded_vectors=exploded_vectors.T
        mode_bfactors=[]
        for numi,i in enumerate(self.eigen_values[6:]):
            evec_row=exploded_vectors[numi+6]
            mode_bfactors.append([ float(evec_row[j]**2 + evec_row[j+1]**2 + evec_row[j+2]**2)/i for j in range(0,len(x0)*3,3)])
        
        mode_bfactors=numpy.array(mode_bfactors)
        self.fluctuations = [numpy.sum(i) for i in mode_bfactors.T]
        return True
    

    def calculate_movie(self,mode_number,scale=1.5,n=10):
        x0=numpy.array([i.get_location() for i in self.atoms])
        d0=[i.get_domain_id() for i in self.atoms]
        new_coords=[]
        with open(str(mode_number)+'.pdb','w') as fh:
            for j in [k for k in range(n)]+[k for k in range(n)[::-1]]:
                HingeResidue=0
                for numi,i in enumerate(x0):
                    if(d0[numi][0]=='D'):
                        D_delta_phi= self.eigen_vectors[:,mode_number][self.domain_info[d0[numi]][0]*6 : (self.domain_info[d0[numi]][0]*6)+6]
                        D_COM=self.domain_info[d0[numi]][1]
                        new_x=i[0]+scale*j* ( D_delta_phi[0] + (D_delta_phi[4]*(i[2]-D_COM[2])) - (D_delta_phi[5]*(i[1]-D_COM[1])) )
                        new_y=i[1]+scale*j* ( D_delta_phi[1] + (D_delta_phi[3]*(i[2]-D_COM[2])) - (D_delta_phi[5]*(i[0]-D_COM[0])) )
                        new_z=i[2]+scale*j* ( D_delta_phi[2] + (D_delta_phi[3]*(i[1]-D_COM[1])) - (D_delta_phi[4]*(i[0]-D_COM[0])) )
                    if(d0[numi][0]=='H'):
                        delta=self.eigen_vectors[:,mode_number][6*len(self.domain_info):][HingeResidue*3:(HingeResidue*3)+3]
                        new_x=i[0]+scale*j*delta[0]
                        new_y=i[1]+scale*j*delta[1]
                        new_z=i[2]+scale*j*delta[2]
                        HingeResidue+=1
                    new_x , new_y, new_z = new_x.real , new_y.real , new_z.real
                    new_coords.append([new_x,new_y,new_z])
                    fh.write("ATOM  %5s %-4s %3s %1s%4s    %8s%8s%8s%6s%6s         %-4s%2s%2s\n"%(self.atoms[numi].get_id(),self.atoms[numi].get_name(),self.atoms[numi].get_parent().get_name(),self.atoms[numi].get_parent().get_parent().get_id(),self.atoms[numi].get_parent().get_id(),round(new_x,3),round(new_y,3),round(new_z,3),self.atoms[numi].get_occupancy(),self.atoms[numi].get_bfactor(),'',self.atoms[numi].get_element(),''))
                    
                fh.write('ENDMDL')
        return True

    def calculate_new_movie(self,mode_number,scale=1.5,n=10):
        x0=numpy.array([i.get_location() for i in self.atoms])
        d0=[i.get_domain_id() for i in self.atoms]
        new_coords=[]
        with open(str(mode_number)+'.pdb','w') as fh:
            for j in [k for k in range(n)]+[k for k in range(n)[::-1]]:
                HingeResidue=0
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
                        


                    if(d0[numi][0]=='H'):
                        delta=self.eigen_vectors[:,mode_number][6*len(self.domain_info):][HingeResidue*3:(HingeResidue*3)+3]
                        new_x=i[0]+scale*j*delta[0]
                        new_y=i[1]+scale*j*delta[1]
                        new_z=i[2]+scale*j*delta[2]
                        HingeResidue+=1
                    new_x , new_y, new_z = new_x.real , new_y.real , new_z.real
                    new_coords.append([new_x,new_y,new_z])
                    fh.write("ATOM  %5s %-4s %3s %1s%4s    %8s%8s%8s%6s%6s         %-4s%2s%2s\n"%(self.atoms[numi].get_id(),self.atoms[numi].get_name(),self.atoms[numi].get_parent().get_name(),self.atoms[numi].get_parent().get_parent().get_id(),self.atoms[numi].get_parent().get_id(),round(new_x,3),round(new_y,3),round(new_z,3),self.atoms[numi].get_occupancy(),self.atoms[numi].get_bfactor(),'',self.atoms[numi].get_element(),''))
                    
                fh.write('ENDMDL')
        return True