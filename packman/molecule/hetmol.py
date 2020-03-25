'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

import numpy


class HetMol():
    def __init__(self,id,name,parent):
        self.__id=id
        self.__name=name
        self.__parent=parent
        self.__Atoms=None
    
    def __setitem__(self,id,Atom):
        try:
            self.__Atoms[id]=Atom
        except:
            self.__Atoms={}
            self.__Atoms[id]=Atom
    
    #Get Functions
    def get_id(self):
        """
        :returns: ID of an 'Residue' object
        """
        return self.__id
    
    def get_name(self):
        """
        :returns: Name of an 'Residue' object
        """
        return self.__name

    def get_parent(self):
        """
        :returns: Parent residue of an 'Residue' Object. (If available)
        """
        return self.__parent
    
    def get_hetatoms(self):
        """
        :returns: Generator yielding all the 'Atom' objects in a residue.
        """
        for i in sorted(self.__Atoms.keys()):yield self.__Atoms[i]
    
    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    def set_name(self,new_name):
        """
        :param new_name: New name to be assigned
        """
        self.__name=new_name
    
    def set_parent(self,new_parent):
        """
        :param new_parent: New parent to be assigned
        """
        self.__parent=new_parent

    def get_centerofgravity(self):
        """
        NOTE: Yet to add the atomic masses.
        :returns: Cartesian Coordinates of the centre of the gravity.
        """
        atoms=self.get_hetatoms()
        AtomicMass=1
        XYZ_M=[0,0,0]
        MassofAA=0
        for i in atoms:
            XYZ_M[0]+=i.Coordinates[0]*AtomicMass
            XYZ_M[1]+=i.Coordinates[1]*AtomicMass
            XYZ_M[2]+=i.Coordinates[2]*AtomicMass
            MassofAA=MassofAA+AtomicMass
        return numpy.array([i/MassofAA for i in XYZ_M])
