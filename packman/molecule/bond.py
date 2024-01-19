# -*- coding: utf-8 -*-
"""The 'Bond' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Bond' object documentation for details.

Citation:
    Pranav M Khade, Robert L Jernigan, PACKMAN-Molecule: Python Toolbox for Structural Bioinformatics, Bioinformatics Advances, 2022;, vbac007, https://doi.org/10.1093/bioadv/vbac007

Example::

    from packman.molecule import Bond
    help( Bond )

Authors:
    * Pranav Khade (https://github.com/Pranavkhade)
"""
import logging

from typing import TYPE_CHECKING, Union, Tuple

if(TYPE_CHECKING):
    from . import Atom

class Bond():
    """This class contains the information about the 'Bond' object (packman.molecule.Bond).

    This class contains all the information available about the Bond. The Bond class is the connection between the lowest in the hierarchy of the 'molecule' API classes (Atom).
    the order of hierarchy being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
    Please read the Tutorials and Documentation for more details.

    Note:
        - Since this is an object file, user is allowed to bond any two atoms in the structure. However, it is recommended to use caution and stick to the conventions when bonding two atoms.
        - The bonds that are not in the CONNECT records of the structure file are calculated using the following algorithm:
            Zhang, Q., Zhang, W., Li, Y. et al. A rule-based algorithm for automatic bond type perception. J Cheminform 4, 26 (2012). https://doi.org/10.1186/1758-2946-4-26
    
    Args:
        id (int)                      : ID of the bond.
        atom1 (packman.molecule.Atom) : The first atom connected by the bond.  (order is not important)
        atom2 (packman.molecule.Atom) : The second atom connected by the bond. (order is not important)
        source (str)                  : Source of information for the given bond.
        type (str)                    : Bond can be either of the following type: (covalent, ionic, hydrogen or other)
        electrons (tuple)             : A python tuple containing information about shared electrons (eg.... (1,1) means the electrons shared between a single hydrocarbon bond if the atom1 is carbon and atom2 is hydrogen touple order here is (atom1, atom2).
    """
    def __init__(self, id: int, atom1: 'Atom', atom2: 'Atom', type: str, source: Union[str, None]=None):
        self.__id = id
        self.__allowed_bond_types = ['non-covalent', 'covalent', 'covalent-single', 'covalent-double', 'covalent-triple' , 'ionic', 'hydrogen', 'salt-bridge', 'other']
        self.__atom1 = atom1
        self.__atom2 = atom2
        self.__type = type
        self.__source = source

        #Properties are the entities that are not included in the PDB files and are obtained by calculations
        self.__properties = {}
    
    def __repr__(self) -> str:
        return '<Bond ('+str(self.__type)+') '+str(self.__id)+' between Atom '+str(self.__atom1.get_id())+' and '+str(self.__atom2.get_id())+'>'
    
    #Get Functions
    def get_id(self) -> int:
        """Get the ID of the 'Bond'

        Returns:
            id (int) : ID of the bond.
        """
        return self.__id

    def get_atoms(self) -> 'Atom':
        """Get the atoms involved in the bond.

        Returns:
            (atom1, atom2) if successful; None otherwise
        """
        return ( self.__atom1, self.__atom2 )
    
    def get_type(self) -> str:
        """Get the type of the bond

        Returns:
            type of the bond (str) if successful; None otherwise
        """
        return self.__type
    
    def get_source(self) -> str:
        """Get the source of the information of the created bond.

        Returns:
            source (str): Source of information for the given bond.
        """
        return self.__source
    
    def get_property(self, property_name):
        """Get the Property of the given 'Bond'.

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
    
    def get_torsion(self, neighbor1: Union[int, 'Atom'], neighbor2: Union[int, 'Atom'], radians: bool=True) -> float:
        """Calculate the torsion angle of the given covalent bond with the corresponding selected neighbors.

        Note:
            At least four atoms are needed to form two planes that can measure the torsional angles; therefore, along with the two bond atoms, the user needs to provide the additional two atoms that are ideally non-mutual neighbors of the atoms in the bond.

        Args:
            neighbor1 (int/packman.molecule.Atom) : Neighbour of the Atom1 as an 'Atom' object or Atom ID.
            neighbor1 (int/packman.molecule.Atom) : Neighbour of the Atom2 as an 'Atom' object or Atom ID.
            radians   (bool)                      : Return value of the angle in radians (returns value in degrees if False; Default : True)
        
        Returns:
            The torsion angle in radians/degrees if sucessful, None otherwise.
        """
        Model = self.__atom1.get_parent().get_parent().get_parent()
        return Model.get_torsion(self, neighbor1=neighbor1, neighbor2=neighbor2, radians=radians)
    
    #Set Functions
    def set_atoms(self, new_pair: Tuple['Atom', 'Atom']) -> bool:
        """Set the atoms involved in the bond.
        
        Args:
            atom_pair (tuple): Tuple of :ref:`packman.molecule.Atom` objects.
        
        Returns:
            True if successful, False otherwise.
        """
        if(len(new_pair)!=2):
            logging.warning('The input provided for Bond.set_atoms() does not have valid length. Please make sure that exactly two atoms are provided. eg.. (Atom1, Atom2)')
            return False
        else:
            self.__atom1 = new_pair[0]
            self.__atom2 = new_pair[1]
            return True
    
    def set_type(self, new_type: str) -> bool:
        """Set the bond type.

        Bond can be either of the following type: (covalent, ionic, hydrogen or other)

        Args:
            new_type (str): New bond type
        
        Returns:
            True if successful, False otherwise.
        """
        if(new_type in self.__allowed_bond_types):
            self.__type = new_type
            return True
        else:
            logging.warning( 'Invalid bond type. Allowed types: '+', '.join(self.__allowed_bond_types) )
            return False

    def set_property(self, property_name, value) -> bool:
        """Set the Property of the given 'Bond'.

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
            return True
        except:
            logging.warning('Please check the property name. Check the allowed Python dictionary key types for more details.')
            return False
    
    def set_torsion(self, theta: float, neighbor1: Union[int, 'Atom'], neighbor2: Union[int, 'Atom'], radians: bool=True):
        """Set the torsion for the given covalent bond with the corresponding selected neighbors.

        Note:
            At least four atoms are needed to form two planes that change the torsional angles; therefore, along with the two bond atoms, the user needs to provide the additional two atoms that are ideally non-mutual neighbors of the atoms in the bond.

        Args:
            theta     (float)                     : Set the torsional angle (see the 'radians' parameter description)
            neighbor1 (int/packman.molecule.Atom) : Neighbour of the Atom1 as an 'Atom' object or Atom ID.
            neighbor1 (int/packman.molecule.Atom) : Neighbour of the Atom2 as an 'Atom' object or Atom ID.
            radians   (bool)                : Parameter 'theta' will be assuned to be in Radians if True, Degrees will be assumed when False. ( Default : True)
        
        Returns:
            True if successful, None otherwise
        """
        Model = self.__atom1.get_parent().get_parent().get_parent()
        return Model.set_torsion(self, theta, neighbor1=neighbor1, neighbor2=neighbor2, radians=radians)

    #Calculate Functions