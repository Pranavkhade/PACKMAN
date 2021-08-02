# -*- coding: utf-8 -*-
"""The 'Bond' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Bond' object documentation for details.

Example::

    from packman.molecule import Bond
    help( Bond )


Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""
import logging

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
    def __init__(self, id, atom1, atom2, type, source=None):
        self.__id = id
        self.__allowed_bond_types = ['covalent', 'covalent-single', 'covalent-double', 'covalent-triple' , 'ionic', 'hydrogen', 'other']
        self.__atom1 = atom1
        self.__atom2 = atom2
        self.__type = type
        self.__source = source

        #Properties are the entities that are not included in the PDB files and are obtained by calculations
        self.__properties = {}
    
    #Get Functions
    def get_id(self):
        """Get the ID of the 'Bond'

        Returns:
            id (int) : ID of the bond.
        """
        return self.__id

    def get_atoms(self):
        """Get the atoms involved in the bond.

        Returns:
            (atom1, atom2) if successful; None otherwise
        """
        return ( self.__atom1, self.__atom2 )
    
    def get_type(self):
        """Get the type of the bond

        Returns:
            type of the bond (str) if successful; None otherwise
        """
        return self.__type
    
    def get_source(self):
        """Get the source of the information of the created bond.

        Returns:
            source (str): Source of information for the given bond.
        """
        return self.__source
    
    def get_property(self,property_name):
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
    
    #Set Functions
    def set_atoms(self, new_pair):
        """Set the atoms involved in the bond.
        
        Args:
            atom_pair (tuple): Tuple of :ref:`packman.molecule.Atom` objects.
        """
        if(len(new_pair)!=2):
            logging.warning('The input provided for Bond.set_atoms() does not have valid length. Please make sure that exactly two atoms are provided. eg.. (Atom1, Atom2)')
        else:
            self.__atom1 = new_pair[0]
            self.__atom2 = new_pair[1]
    
    def set_type(self,new_type):
        """Set the bond type.

        Bond can be either of the following type: (covalent, ionic, hydrogen or other)

        Args:
            new_type (str): 
        """
        if(new_type in self.__allowed_bond_types):
            self.__type = new_type
        else:
            logging.warning( 'Invalid bond type. Allowed types: '+', '.join(self.__allowed_bond_types) )

    def set_property(self,property_name,value):
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
        except:
            logging.warning('Please check the property name. Check the allowed Python dictionary key types for more details.')


    #Calculate Functions