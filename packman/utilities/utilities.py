"""This class contains all the utilities that are neccessary for carriying out structural superimposition, sequence matching etc.


"""

import numpy

def superimporse(reference,target,use='calpha',ids=[],change_target=True):
    """This function is used to superimpose the Target Chain(coordinates will be changed) on the Reference Chain(coordinates will change).

    The superimposition currently is done on the basis of matching Residue ID. If both the proteins have unequal amount of residues,
    still their matching residues will be used for the superimposition. 
    
    It is important to note that sometimes proteins (although same type and chain) have different numbering scheme.
    In such case, the superimposition will not be carried out. Change the IDs of the taget/reference in such a way that it will match each other. 
    
    For more information about how to change ID of the residue, read : packman.molecule.Residue class

    Args:
        reference (packman.molecule.Chain): Chain whose coordinates will remain same and will be used as a reference.
        target (packman.molecule.Chain): Chain whose coordinates will be changed according to the reference chain.  
        use (str): Which atoms to be used for superimposition (Options: calpha, backbone)
        ids (list): Use only particular residues to align (Provide IDs) eg... ids=[1,2,5,77] will use only 1,2,5 and 77th residues to align two chains
        change_target (bool): Change the coordinates of the target chain based on the reference

    Returns:
        R (numpy.matrix): Rotation matrix for Target Chain w.r.t Reference Chain.
        t (numpy.array): Translation vector for Target Chain w.r.t Reference Chain.
    """

    if(use=='calpha'):
        atoms1=[i for i in reference.get_calpha() ]
        atoms2=[i for i in target.get_calpha() ]
    if(use=='backbone'):
        atoms1=[j for i in reference.get_backbone() for j in i]
        atoms2=[j for i in target.get_backbone() for j in i ]    

    #Finding common residues to align
    res1=[i.get_parent().get_id() for i in atoms1]
    res2=[i.get_parent().get_id() for i in atoms2]

    if(ids==[]):
        common_residues=list(set(res1).intersection(res2))
    else:
        common_residues=list( set(res1).intersection(res2).intersection(ids) )

    atoms1=[i for i in atoms1 if i.get_parent().get_id() in common_residues]
    atoms2=[i for i in atoms2 if i.get_parent().get_id() in common_residues]

    atoms1_location=[i.get_location() for i in atoms1]
    atoms2_location=[i.get_location() for i in atoms2]

    #Subtract Mean
    Centroid1=numpy.mean(atoms1_location,axis=0)
    Centroid2=numpy.mean(atoms2_location,axis=0)
    Am= numpy.subtract(atoms1_location,Centroid1)
    Bm= numpy.subtract(atoms2_location,Centroid2)

    #Dot is matrix multiplication for array
    H= numpy.mat(Bm.T) * numpy.mat(Am)

    #Find Rotation
    U, S, Vt = numpy.linalg.svd(H)
    R = Vt.T * U.T

    #Special Reflection Case
    if numpy.linalg.det(R) < 0:
        Vt[2,:] *= -1
        R = Vt.T * U.T
    
    t = numpy.dot(-R , Centroid2) + Centroid1
    
    #Change the location of the target
    if(change_target):
        for i in target.get_atoms():
            new_location=numpy.dot(R , i.get_location()) + t
            i.set_location( numpy.array(new_location.tolist()[0]) )
    else:
        None
    
    return R, t

def RMSD(group1,group2,use='calpha',ids=[]):
    """
    """
    if(use=='calpha'):
        atoms1=[i for i in group1.get_calpha() ]
        atoms2=[i for i in group2.get_calpha() ]
    if(use=='backbone'):
        atoms1=[j for i in group1.get_backbone() for j in i]
        atoms2=[j for i in group2.get_backbone() for j in i ]    

    #Finding common residues to align
    res1=[i.get_parent().get_id() for i in atoms1]
    res2=[i.get_parent().get_id() for i in atoms2]

    if(ids==[]):
        common_residues=list(set(res1).intersection(res2))
    else:
        common_residues=list( set(res1).intersection(res2).intersection(ids) )

    atoms1=[i for i in atoms1 if i.get_parent().get_id() in common_residues]
    atoms2=[i for i in atoms2 if i.get_parent().get_id() in common_residues]

    atoms1_location=[i.get_location() for i in atoms1]
    atoms2_location=[i.get_location() for i in atoms2]

    R,t = superimporse(group1,group2,use,ids,change_target=False)

    for i in range(0,len(atoms2_location)):
        atoms2_location[i] = numpy.dot(R , atoms2_location[i]) + t
        atoms2_location[i] = numpy.array(atoms2_location[i].tolist()[0])
    
    return numpy.mean( [numpy.abs( numpy.linalg.norm(atoms1_location[i]-atoms2_location[i]) ) for i in range(0,len(atoms1_location))] )

def load_hinge(filename):
    """Load the hinge information neccessary for the hd-ANM and other methods.

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


    Args:
        filename (string) : filepath and name of the .hng file
    
    Returns:
        HNGinfo (dictionary) : residue based hinge and domain information.
    """
    HNGinfo={}
    for i in open(filename):
        line=i.strip().split()
        HNGinfo[ line[0]+'_'+line[1] ]=[float(j) for j in line[2].split(':')]
    return HNGinfo

'''
##################################################################################################
#                                    Non Algorithm Functions                                     #
##################################################################################################
'''

def WriteOBJ(atoms,faces, fh):
    """Write the .obj file to visualize the obtain alpha shape tesselations.
    
    Args:
        atoms (packman.molecule.Atom): Atoms (Just for the node records)
        faces ([float])              : SelectedTesselations (See the packman.apps.predict_hinge)
        fh (file)                    : Output file with .obj extension
    
    """
    NewIDs={i.get_id():numi+1 for numi,i in enumerate(atoms)}
    fh.write('mtllib master.mtl\ng\n'.encode())
    fh.write('usemtl atoms\n'.encode())
    for i in atoms:
        x,y,z=i.get_location()
        fh.write("v %f %f %f\n".encode()%(x,y,z))
    
    line='usemtl bonds\nl'
    for i in atoms:
        line=line+" "+str(NewIDs[i.get_id()])
    line=line+'\n'
    fh.write(line.encode())
    
    fh.write('usemtl faces\n'.encode())
    for i in faces:
        faces=[NewIDs[j.get_id()] for j in i]
        fh.write("f %i %i %i %i\n".encode()%(faces[0],faces[1],faces[2],faces[3]))
        #fh.write("l %i %i %i %i\n"%(faces[0],faces[1],faces[2],faces[3]))
    return True
