'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

import numpy


class Model():
    def __init__(self,id,AllAtoms,AllResidues,AllChains,AllHetAtoms,AllHetMols):
        self.__id=id
        self.__AllAtoms=AllAtoms
        self.__AllResidues=AllResidues
        self.__AllChains=AllChains
        self.__AllHetAtoms=AllHetAtoms
        self.__AllHetMols=AllHetMols

    def __getitem__(self,ChainID):
        return self.__AllChains[ChainID]
    
    #Get Functions
    def get_id(self):
        """
        :returns: ID of an 'Model' object
        """
        return self.__id

    def get_chains(self):
        """
        :returns: Constructer for all the chains as a 'Chain' object
        """
        for i in sorted(self.__AllChains.keys()):yield self.__AllChains[i]

    def get_residues(self):
        """
        :returns: Generator yielding all the 'Residue' objects in a Model
        """
        return [j for i in self.__AllChains.keys() for j in self.__AllChains[i].get_residues()]
    
    def get_atoms(self):
        """
        :returns: Generator yielding all the 'Atom' objects in a Model
        """
        for i in sorted(self.__AllAtoms.keys()):yield self.__AllAtoms[i]
    
    def get_chain(self,ChainID):
        """
        :returns: All chains as a 'Chain' objects in a Model
        """
        return self.__AllChains[ChainID]
    
    def get_atom_byid(self,query_atom_id):
        """
        :param query_atom_id: ID of the atom to be returned
        :returns: Atom in from of 'Atom'.
        """
        return self.__AllAtoms[query_atom_id]
    
    def get_hetmols(self):
        """
        :returns: Generator yielding all the 'HetMol' objects in a Model
        """
        for i in sorted(self.__AllHetMols.keys()):yield self.__AllHetMols[i]
     
    def get_hetatoms(self):
        """
        :returns: Generator yielding all the 'HetAtom' objects in a Model
        """
        for i in sorted(self.__AllHetAtoms.keys()):yield self.__AllHetAtoms[i]
    
    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    def set_torsion(self,Atom1,Atom2,Atom3,Atom4,theta):
        """
        NOTE: Currently we can only set torsion for the backbone atoms (Psi and Phi)
        """
        radang=self.get_torsion(Atom1,Atom2,Atom3,Atom4)
        Atom1=self.get_atom_byid(Atom1)
        Atom2=self.get_atom_byid(Atom2)
        Atom3=self.get_atom_byid(Atom3)
        Atom4=self.get_atom_byid(Atom4)
        torsion_type=Atom1.get_name()+Atom2.get_name()+Atom3.get_name()+Atom4.get_name()

        if(torsion_type=='NCACN'):
            rotang=numpy.deg2rad(theta)-radang
        elif(torsion_type=='CNCAC'):
            rotang=radang-numpy.deg2rad(theta)
        else:
            raise Exception('Not a Psi/Phi angle. Current version of PACKMAN only supports Psi/Phi rotation.')
        
        sn=numpy.sin(rotang)
        cs=numpy.cos(rotang)
        t=1-cs

        v2x=Atom2.get_location()[0]-Atom3.get_location()[0]
        v2y=Atom2.get_location()[1]-Atom3.get_location()[1]
        v2z=Atom2.get_location()[2]-Atom3.get_location()[2]
        #Normalize the rotation vector
        mag=numpy.sqrt(numpy.square(v2x)+numpy.square(v2y)+numpy.square(v2z))
        x = float(v2x)/mag
        y = float(v2y)/mag
        z = float(v2z)/mag
        #set up the rotation matrix
        m=numpy.zeros(9)
        m[0]= t*x*x + cs
        m[1]= t*x*y + sn*z
        m[2]= t*x*z - sn*y
        m[3]= t*x*y - sn*z
        m[4]= t*y*y + cs
        m[5]= t*y*z + sn*x
        m[6]= t*x*z + sn*y
        m[7]= t*y*z - sn*x
        m[8]= t*z*z + cs

        #Rotate The Atoms
        tx = Atom2.get_location()[0]
        ty = Atom2.get_location()[1]
        tz = Atom2.get_location()[2]

        import copy
        for i in self.get_atoms():
            Coordinates=i.get_location()
            #Imporve this later on, this skips some atoms
            if(torsion_type=='CNCAC' and i.get_parent().get_id()<=Atom1.get_parent().get_id()):
                Coordinates[0]-=tx
                Coordinates[1]-=ty
                Coordinates[2]-=tz
                x=Coordinates[0]*m[0]+Coordinates[1]*m[1]+Coordinates[2]*m[2]
                y=Coordinates[0]*m[3]+Coordinates[1]*m[4]+Coordinates[2]*m[5]
                z=Coordinates[0]*m[6]+Coordinates[1]*m[7]+Coordinates[2]*m[8]
                Coordinates[0]=x
                Coordinates[1]=y
                Coordinates[2]=z
                Coordinates[0]+=tx
                Coordinates[1]+=ty
                Coordinates[2]+=tz
                i.set_location(numpy.array(Coordinates))
            if(torsion_type=='NCACN' and i.get_parent().get_id()>=Atom4.get_parent().get_id()):
                Coordinates[0]-=tx
                Coordinates[1]-=ty
                Coordinates[2]-=tz
                x=Coordinates[0]*m[0]+Coordinates[1]*m[1]+Coordinates[2]*m[2]
                y=Coordinates[0]*m[3]+Coordinates[1]*m[4]+Coordinates[2]*m[5]
                z=Coordinates[0]*m[6]+Coordinates[1]*m[7]+Coordinates[2]*m[8]
                Coordinates[0]=x
                Coordinates[1]=y
                Coordinates[2]=z
                Coordinates[0]+=tx
                Coordinates[1]+=ty
                Coordinates[2]+=tz
                i.set_location(numpy.array(Coordinates))
        return True
    
    #Compute Functions
    def get_calpha(self):
        """
        :returns: All 'C-Alpha' atoms as an 'Atom' object.
        """
        return [i.get_calpha() for i in self.get_residues()]
    
    def get_backbone(self):
        """
        :returns: All 'Backbone' atoms as an 'Atom object.'
        """
        return [i.get_backbone() for i in self.get_residues()]
    
    def get_torsion(self,Atom1,Atom2,Atom3,Atom4):
        """
        """
        Atom1=self.get_atom_byid(Atom1)
        Atom2=self.get_atom_byid(Atom2)
        Atom3=self.get_atom_byid(Atom3)
        Atom4=self.get_atom_byid(Atom4)
        p0=Atom1.get_location()
        p1=Atom2.get_location()
        p2=Atom3.get_location()
        p3=Atom4.get_location()
        b0 = -1.0*(p1 - p0)
        b1 = p2 - p1
        b2 = p3 - p2
        b1 /= numpy.linalg.norm(b1)
        v = b0 - numpy.dot(b0, b1)*b1
        w = b2 - numpy.dot(b2, b1)*b1
        x = numpy.dot(v, w)
        y = numpy.dot(numpy.cross(b1, v), w)
        radang=numpy.arctan2(y, x)
        return radang
    
    #Check Function
    def check_clashes(self,distance=0.77):
        from scipy.spatial import KDTree
        all_atoms=[i for i in self.get_atoms()]
        T=KDTree([i.get_location() for i in all_atoms])
        return len(T.query_pairs(distance))
