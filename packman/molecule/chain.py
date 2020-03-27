# -*- coding: utf-8 -*-
"""The 'Chain' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Chain' object documentation for details. [ help(packman.molecule.chain.Chain) ]

Example:
    >>>from packman.molecule import chain.Chain
    >>>help( Chain )
    OR

    >>>from packman import molecule
    >>>help( molecule.Chain )

Note:
    * Iterating over chain only fetches the residues not the hetero atoms
    * The object fetches the atoms and residues which are not in order as they appear in the PDB. Find a way to fix this
    
Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""


import numpy


class Chain():
    """This class contains the information about the 'Chain' object (packman.molecule.Chain).

    This class contains all the information available about the Chain and stores the corresponding 'Residue' and 'Hetmol' objects in itself. The Chain class is the third lowest in the hierarchy of the 'molecule' API classes.
    the order of hierarchy being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
    Please read the Tutorials and Documentation for more details.

    Note:
        * Please refer to the [https://web.archive.org/web/20080905024351/http://www.wwpdb.org/docs.html] for the description of the arguments.
        * Add get_atoms()
    
    Args:
        id (str): Chain ID from the PDB file as it is. Each Chain in a PDB file Model/Frame has unique ID. (essential)
        
    """
    
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
        """Get the ID of the 'Residue'

        Returns:
            str if successful, None otherwise.
        """
        return self.__id
    
    def get_parent(self):
        """Get the 'Model' the given 'Chain' belongs to.

        Returns:
            packman.molecule.Model if successful, None otherwise.
        """
        return self.__parent
    
    def get_hinges(self):
        """Get the hinges in the chain as a dictionary of 'Hinge' objects.

        Returns:
            packman.molecule.annotations.Hinge if successful, None otherwise.
        """
        return self.__Hinges
    
    #Calculation Functions
    def get_residues(self):
        """Get the generator of corresponding 'Residue' objects of the 'Chain'

        Returns:
            generator of 'Residue' objects if successful, None otherwise.
        """
        for i in sorted(self.__Residues.keys()):yield self.__Residues[i]
    
    def get_calpha(self):
        """Get the C-Alpha atoms of the 'Chain' as an 'Atom' object.

        Returns:
            generator of packman.molecule.Atom objects if successful, None otherwise.
        """
        for i in self.get_residues():yield i.get_calpha()
    
    def get_hetmols(self):
        """Get the generator of corresponding 'HetMol' objects of the 'Chain'

        Returns:
            generator of 'HetAtom' objects if successful, None otherwise.
        """
        for i in sorted(self.__HetMols.keys()):yield self.__HetMols[i]
    
    def get_backbone(self):
        """Get the Backbone atoms of the given 'Chain' as a list of 'Atom' object

        Note:
            * Backbone Atoms: CA, O, N, C

        Returns:
            list of packman.molecule.Atom if successful, None otherwise.
        """
        return [i.get_backbone() for i in self.get_residues()]


    

    #Set Functions
    def set_id(self,new_id):
        """Set the ID of the given 'Chain'

        Args:
            new_id (str): The ID User wishes to assign to the given 'Chain'
        """
        self.__id=new_id
    
    def set_parent(self,parent):
        """Set the Parent of the given 'Residue'

        Args:
            new_parent (packman.molecule.Chain): The parent 'Chain' User wishes to assign to the given 'Residue'
        """
        self.__parent=parent
    
    def set_hinges(self,new_hinges):
        """Set/Add hinge to the 'Chain' object

        Args:
            new_hinges (packman.molecule.annotations.Hinge): The 'Hinge' User wishes to assign/add to the given 'Chain'
        """
        self.__Hinges=self.__Hinges+new_hinges