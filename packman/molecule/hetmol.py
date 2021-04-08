# -*- coding: utf-8 -*-
"""The 'HetMol' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'HetMol' object documentation for details.

Example::

    from packman.molecule import HetMol
    help( HetMol )

Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""

import numpy


class HetMol():
    """This class contains the information about the 'HetMol' object (packman.molecule.Residue).

    This class contains all the information available about the HetMol and stores the corresponding 'Atom' objects in itself. The Residue class is the second lowest in the hierarchy of the 'molecule' API classes.
    the order of hierarchy being: Protein> Model> Chain> HetMol> Atom. This class is also the component of the 'molecule' module API.
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

        #Properties are the entities that are not included in the PDB files and are obtained by calculations
        self.__properties = {}
    
    def __setitem__(self,id,Atom):
        """Simply assign new atom to the 'HetMol'

        Args:
            new_id (int)                  : The ID User wishes to assign to the given 'Atom'
            Atom  (packman.molecule.Atom) : The new 'Atom' to be assigned to the 'HetMol'
        """
        try:
            self.__Atoms[id]=Atom
        except:
            self.__Atoms={}
            self.__Atoms[id]=Atom
    
    #Get Functions
    def get_id(self):
        """Get the ID of the 'HetMol'

        Returns:
            int if successful, None otherwise.
        """
        return self.__id
    
    def get_name(self):
        """Get the Name of the 'HetMol'

        Returns:
            str if successful, None otherwise.
        """
        return self.__name

    def get_parent(self):
        """Get the Parent of the 'HetMol'

        Returns:
            packman.molecule.Chain if successful, None otherwise.
        """
        return self.__parent
    
    def get_atoms(self):
        """Get the generator of corresponding 'Atom' objects of the 'HetMol'

        Returns:
            generator of 'Atom' objects if successful, None otherwise.
        """
        for i in sorted(self.__Atoms.keys()):yield self.__Atoms[i]
    
    def get_domain_id(self):
        """Get the Domain Identifier of the given 'HetMol'. Hinge Prediction is Necessary for this option.

        Note:
            Domain Identifiers (Obtained after running PACKMAN):
            'FL': Flexible Linkers (Hinges); Followed by the domain number
            'DM': Domain; Followed by the domain number

        Returns:
            str if successful, None otherwise.
        """
        return self.__domain_id
    
    def get_property(self,property_name):
        """Get the Property of the given 'HetMol'.

        Property is any key and value combination that can be assigned to this object. This (along with the set_property) feature is mainly useful for the user customization.
        Properties are like pinboards. You can pin anything to the object with a key as a pin.

        Args:
            property_name (object): The 'Key' or a name the user wants to assign to to the property
        
        Note:
            - Users can add custom annotations; for example: If particular chain becomes disordered, it can be annotated with this feature.
        """
        try:
            return self.__properties[property_name]
        except:
            logging.warning('The Property Name provided is not assigned.')
    
    #Set Functions
    def set_id(self,new_id):
        """Set the ID of the given 'HetMol'

        Args:
            new_id (int): The ID User wishes to assign to the given 'HetMol'
        """
        self.__id=new_id
    
    def set_name(self,new_name):
        """Set the Name of the given 'HetMol'

        Args:
            new_name (str): The Name User wishes to assign to the given 'HetMol'
        """
        self.__name=new_name
    
    def set_parent(self,new_parent):
        """Set the Parent of the given 'HetMol'

        Args:
            new_parent (packman.molecule.Chain): The parent 'Chain' User wishes to assign to the given 'HetMol'
        """
        self.__parent=new_parent

    def set_domain_id(self,new_domain_id):
        """Set the Domain Identifier of the given 'HetMol'.

        Note:
            Domain Identifiers Format:
            'FL': Flexible Linkers (Hinges); Followed by the domain number
            'DM': Domain; Followed by the domain number

        Args:
            new_domain_id (str): The Domain Identifies User wishes to assign to the given 'HetMol'
        """
        self.__domain_id=new_domain_id
    
    def set_property(self,property_name,value):
        """Set the Property of the given 'HetMol'.

        Property is any key and value combination that can be assigned to this object. This (along with the get_property) feature is mainly useful for the user customization.
        Properties are like pinboards. You can pin anything to the object with a key as a pin.
        
        Args:
            property_name (object): The 'Key' or a name the user wants to assign to to the property
            value (object):         The value the user wants to assign to the property
        
        Note:
            - Users can add custom annotations; for example: If particular amino acid becomes disordered, it can be annotated with this feature.
        """
        try:
            self.__properties[property_name] = value
        except:
            logging.warning('Please check the property name. Check the allowed Python dictionary key types for more details.')

    #Calculation Functions
    def get_centerofgravity(self):
        """Get the center of gravity of the given 'HetMol'

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