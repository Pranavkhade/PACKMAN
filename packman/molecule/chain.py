'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

import numpy


class Chain():
    '''
    NOTE1: As of now, iterating over chain only fetches the residues not the hetero atoms
    NOTE2: The object fetches the atoms and residues which are not in order as they appear in the PDB. Find a way to fix this
    NOTE3: Chain has no get parent as of now | Resolved
    '''
    def __init__(self,id):
        self.__id=id
        self.__Residues=None
        self.__HetMols=None
        self.__parent=None
        #More Features
        self.__Hinges=[]
    
    def __setitem__(self,Number,Entity,Type):
        if(Type=='Residue'):
            try:
                self.__Residues[Number]=Entity
            except:
                self.__Residues={}
                self.__Residues[Number]=Entity
        elif(Type=='HetMol'):
            try:
                self.__HetMols[Number]=Entity
            except:
                self.__HetMols={}
                self.__HetMols[Number]=Entity

    def __getitem__(self,Number,Type=None):
        try:
            return self.__Residues[Number]
        except:
            return self.__HetMols[Number]
    
    #Get Functions
    def get_id(self):
        """
        :returns: ID of an 'Chain' object
        """
        return self.__id
    
    def get_residues(self):
        """
        :returns: Generator yielding all the 'Residue' objects in a Chain.
        """
        for i in sorted(self.__Residues.keys()):yield self.__Residues[i]
    
    def get_calpha(self):
        """
        :returns: Generator yielding all the 'C-Alpha' Atom objects in a Chain.
        """
        for i in self.get_residues():yield i.get_calpha()
    
    def get_hetmols(self):
        """
        :returns: Generator yielding all the 'HetMol' objects in a Chain.
        """
        for i in sorted(self.__HetMols.keys()):yield self.__HetMols[i]

    
    def get_hinges(self):
        """
        :returns: The hinge residues
        """
        return self.__Hinges
    
    def get_parent(self):
        """
        """
        return self.__parent

    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    def set_hinges(self,new_hinges):
        """
        """
        self.__Hinges=self.__Hinges+new_hinges
    
    def set_parent(self,parent):
        """
        """
        self.__parent=parent
    
    #Compute Functions
    def get_backbone(self):
        """
        :returns: All 'Backbone' atoms as an 'Atom object.'
        """
        return [i.get_backbone() for i in self.get_residues()]
