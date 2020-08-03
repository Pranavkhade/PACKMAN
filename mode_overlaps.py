from packman.molecule import load_structure

from packman.utilities import superimporse

from packman.molecule import Protein

from packman.anm import RDANM

from scipy.spatial.transform import Rotation as SciPyRotation

import numpy

import sys



def load_hinges(filename):
    """
    Note: Its important not to leave the hng file open ended. ie.. Remove 'Inf'
    """

    HNGinfo={}
    for i in open(filename):
        line=i.strip().split()
        HNGinfo[ line[0]+'_'+line[1] ]=[float(j) for j in line[2].split(':')]
    return HNGinfo


def main():

    #Input format: abc.pdb_A (filename_chain)
    filename1, chain1 = sys.argv[1].split('_')
    filename2, chain2 = sys.argv[2].split('_')
    hngfile = sys.argv[3]

    #File Change
    mol1=load_structure(filename1)
    mol2=load_structure(filename2)

    #Chain Change
    R,T=superimporse(mol1[0][chain1],mol2[0][chain2],use='calpha')

    #Change Hinge
    hng_info=load_hinges(hngfile)

    #Calpha
    calpha1=[i for i in mol1[0][chain1].get_calpha()]
    calpha2=[i for i in mol2[0][chain2].get_calpha()]

    common_residues = list(set([j.get_parent().get_id() for j in calpha1]).intersection([j.get_parent().get_id() for j in calpha2]))
    calpha1 = [i for i in calpha1 if i.get_parent().get_id() in common_residues]
    calpha2 = [i for i in calpha2 if i.get_parent().get_id() in common_residues]


    #Vector Tracing
    compare_domain=[]
    compare_hinge =[]
    for numi,i in enumerate(hng_info):
        ResPos = i.split('_')[2][0]
        ResRNG = hng_info[i]
        if(ResPos=='D'):
            #print(ResRNG,numpy.arange(ResRNG[0],ResRNG[1]+1,1))
            #Calculate Rotation and translation of the superimposed domain
            R_D,t_D = superimporse(mol1[0]['A'], mol2[0]['A'], ids=numpy.arange(ResRNG[0],ResRNG[1]+1,1), change_target=False)
            
            #Convert rotation to quarternion
            rot=SciPyRotation.from_matrix(R_D).as_rotvec().tolist()
            
            #Normalize with axis of rotation
            rot=rot/ numpy.linalg.norm(rot)
            
            compare_domain.extend(t_D.tolist()[0])
            compare_domain.extend(rot)

        elif(ResPos=='H'):
            common_residues = list(set([j.get_parent().get_id() for j in calpha1]).intersection([j.get_parent().get_id() for j in calpha2]))

            calpha1_H=numpy.array([j.get_location() for j in calpha1 if  ResRNG[0] <= j.get_parent().get_id() <= ResRNG[1] ])#and j.get_parent().get_id() in common_residues ])
            calpha2_H=numpy.array([j.get_location() for j in calpha2 if  ResRNG[0] <= j.get_parent().get_id() <= ResRNG[1] ])#and j.get_parent().get_id() in common_residues ])

            dir_vec=calpha1_H-calpha2_H
            dir_vec=dir_vec.flatten().tolist()
            compare_hinge.extend(dir_vec)
    

    #IMPORTANT VECTOR
    compare=compare_domain+compare_hinge

    MODEL=RDANM(calpha1,dr=15.0,power=0,HNGinfo=hng_info)
    MODEL.calculate_coarse_grained_hessian()
    MODEL.calculate_decomposition()
    
    EIGVEC=MODEL.get_eigenvectors()
    EIGVAL=MODEL.get_eigenvalues()
    MODEL.calculate_movie(4)

    numpy.savetxt('EVEC.txt',EIGVEC[:,5:9])
    numpy.savetxt('OVRLP.txt',compare)

    score=[]
    for numi,i in enumerate(EIGVEC.T):
        i=EIGVAL[numi]*i
        score.append( numpy.abs(numpy.dot(compare,i)) / ( numpy.sqrt(numpy.dot(compare,compare)) * numpy.sqrt(numpy.dot(i,i)) ) )
    
    for i in score:
        print(i)

    return True


if(__name__ == "__main__"):
    main()