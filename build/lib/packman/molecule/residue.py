'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

import numpy


class Residue():
    '''
    Domain Identifiers: (Obtained after running PACKMAN)
    FL: Flexible Linkers (Hinges); Followed by a number
    DM: Domain; Followed by a number
    '''
    def __init__(self,id,name,parent):
        self.__id=id
        self.__name=name
        self.__parent=parent
        self.__Atoms=None
        self.__domain_id=None
    
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
    
    def get_atoms(self):
        """
        :returns: Generator yielding all the 'Atom' objects in a residue.
        """
        for i in sorted(self.__Atoms.keys()):yield self.__Atoms[i]
    
    def get_domain_id(self):
        '''
        :returns: The assigned domain ID of the residue
        '''
        return self.__domain_id
    
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

    def set_domain_id(self,new_domain_id):
        self.__domain_id=new_domain_id

    #Calculation Functions
    def get_calpha(self):
        """
        :returns: A 'C-Alpha' atom as an 'Atom' object.
        """
        return [i for i in self.get_atoms() if i.get_name()=='CA'][0]
    
    def get_centerofgravity(self):
        """
        NOTE: Yet to add the atomic masses.
        :returns: Cartesian Coordinates of the centre of the gravity.
        """
        atoms=self.get_atoms()
        AtomicMass=1
        XYZ_M=[0,0,0]
        MassofAA=0
        for i in atoms:
            XYZ_M[0]+=i.Coordinates[0]*AtomicMass
            XYZ_M[1]+=i.Coordinates[1]*AtomicMass
            XYZ_M[2]+=i.Coordinates[2]*AtomicMass
            MassofAA=MassofAA+AtomicMass
        return numpy.array([i/MassofAA for i in XYZ_M])
    
    def get_backbone(self):
        """
        NOTE: Some files like 1k20 are showing multiple backbone atoms because seperation of models is based on [Model] in this tool
        :returns: 'Backbone' atoms as an 'Atom' object.
        """
        return [i for i in self.get_atoms() if(i.get_name()=='N' or i.get_name()=='CA' or i.get_name()=='C' or i.get_name()=='O')]
    
    def get_tip(self):
        """
        :returns: Tip atom of an amino acid as an 'Atom' object.
        """
        CAlpha=self.get_calpha()
        resname=self.get_name()
        TipofAA=None
        if(resname=='ALA' or resname=='GLY'):
            TipofAA=CAlpha
        else:
            MaxDistance=0
            for i in self.get_atoms():
                tempdistance=CAlpha.calc_dist(i)
                if(tempdistance>MaxDistance):
                    MaxDistance=tempdistance
                    TipofAA=i
        return TipofAA
