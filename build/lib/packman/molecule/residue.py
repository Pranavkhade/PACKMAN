# -*- coding: utf-8 -*-
"""The 'Residue' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Residue' object documentation for details. [ help(packman.molecule.residue.Residue) ]

Example:
    >>>from packman.molecule import residue.Residue
    >>>help( Residue )
    OR

    >>>from packman import molecule
    >>>help( molecule.Residue )

Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""

import numpy


class Residue():
    """This class contains the information about the 'Residue' object (packman.molecule.Residue).

    This class contains all the information available about the Residue and stores the corresponding 'Atom' objects in itself. The Residue class is the second lowest in the hierarchy of the 'molecule' API classes.
    the order of hierarchy being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
    Please read the Tutorials and Documentation for more details.

    Note:
        - Please refer to the [https://web.archive.org/web/20080905024351/http://www.wwpdb.org/docs.html] for the description of the arguments.
    
    Args:
        id (int): Residue ID from the PDB file as it is. Each Residue in a PDB file Model/Frame has unique ID. (essential)
        name (str): Residue Name from the PDB file as it is.
        parent (packman.molecule.Chain): The Chain Object (parent) this Residue belongs to.

    """
    def __init__(self,id,name,parent):
        self.__id=id
        self.__name=name
        self.__parent=parent
        self.__Atoms=None
        self.__domain_id=None
    
    def __setitem__(self,id,Atom):
        """Simply assign new atom to the 'Residue'

        Args:
            new_id: The ID User wishes to assign to the given 'Atom'
        """
        try:
            self.__Atoms[id]=Atom
        except:
            self.__Atoms={}
            self.__Atoms[id]=Atom
    
    #Get Functions
    def get_id(self):
        """Get the ID of the 'Residue'

        Returns:
            int if successful, None otherwise.
        """
        return self.__id
    
    def get_name(self):
        """Get the Name of the 'Residue'

        Returns:
            str if successful, None otherwise.
        """
        return self.__name

    def get_parent(self):
        """Get the Parent of the 'Residue'

        Returns:
            packman.molecule.Chain if successful, None otherwise.
        """
        return self.__parent
    
    def get_atoms(self):
        """Get the generator of corresponding 'Atom' objects of the 'Residue'

        Returns:
            generator of 'Atom' objects if successful, None otherwise.
        """
        for i in sorted(self.__Atoms.keys()):yield self.__Atoms[i]
    
    def get_domain_id(self):
        """Get the Domain Identifier of the given 'Residue'. Hinge Prediction is Necessary for this option.

        Note:
            Domain Identifiers (Obtained after running PACKMAN):
            'FL': Flexible Linkers (Hinges); Followed by the domain number
            'DM': Domain; Followed by the domain number

        Returns:
            str if successful, None otherwise.
        """
        return self.__domain_id
    
    #Set Functions
    def set_id(self,new_id):
        """Set the ID of the given 'Residue'

        Args:
            new_id (int): The ID User wishes to assign to the given 'Residue'
        """
        self.__id=new_id
    
    def set_name(self,new_name):
        """Set the Name of the given 'Residue'

        Args:
            new_name (str): The Name User wishes to assign to the given 'Residue'
        """
        self.__name=new_name
    
    def set_parent(self,new_parent):
        """Set the Parent of the given 'Residue'

        Args:
            new_parent (packman.molecule.Chain): The parent 'Chain' User wishes to assign to the given 'Residue'
        """
        self.__parent=new_parent

    def set_domain_id(self,new_domain_id):
        """Set the Domain Identifier of the given 'Residue'.

        Note:
            Domain Identifiers Format:
            'FL': Flexible Linkers (Hinges); Followed by the domain number
            'DM': Domain; Followed by the domain number

        Args:
            new_domain_id (str): The Domain Identifies User wishes to assign to the given 'Residue'
        """
        self.__domain_id=new_domain_id

    #Calculation Functions
    def get_calpha(self):
        """Get the C-Alpha atom of the residue as an 'Atom' object.

        Returns:
            packman.molecule.Atom if successful, None otherwise.
        """
        return [i for i in self.get_atoms() if i.get_name()=='CA'][0]
    
    def get_centerofgravity(self):
        """Get the center of gravity of the given 'Residue'

        Note:
            Yet to add the atomic masses.

        Returns:
            Cartesian Coordinates as numpy.array of the centre of the gravity.
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
        """Get the Backbone atoms of the given 'Residue' as a list of 'Atom' object

        Note:
            Some files like 1k20 are showing multiple backbone atoms because seperation of models is based on [Model] in this tool (Maybe Solved)

        Returns:
            list of packman.molecule.Atom if successful, None otherwise.
        """
        return [i for i in self.get_atoms() if(i.get_name()=='N' or i.get_name()=='CA' or i.get_name()=='C' or i.get_name()=='O')]
    
    def get_tip(self):
        """Get the tip atom of the given 'Residue' as an 'Atom' object

        Note:
            ALA and GLY are small so their tip is C-alpha

        Returns:
            packman.molecule.Atom if successful, None otherwise.
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
