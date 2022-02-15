# -*- coding: utf-8 -*-
"""The 'Atom' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Atom' object documentation for details.

Citation:
    Pranav M Khade, Robert L Jernigan, PACKMAN-Molecule: Python Toolbox for Structural Bioinformatics, Bioinformatics Advances, 2022;, vbac007, https://doi.org/10.1093/bioadv/vbac007

Example::

    from packman.molecule import Atom
    help( Atom )


Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""

import numpy

import logging

from packman import molecule


class Atom():
    """This class contains the information about the 'Atom' object (packman.molecule.Atom).

    This class contains all the information available about the Atom. The Atom class is the lowest in the hierarchy of the 'molecule' API classes.
    the order of hierarchy being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
    Please read the Tutorials and Documentation for more details.

    Note:
        - Most essential class of all as it stores the actual data of the atoms.
        - Please refer to the [https://web.archive.org/web/20080905024351/http://www.wwpdb.org/docs.html] for the description of the arguments.
    
    Args:
        id (int): Atom ID from the PDB file as it is. Each Atom in a PDB file Model/Frame has unique ID. (essential)
        AtomName (str): Atom Name from the PDB file as it is.
        Coordinates (numpy.array): The 3*1 array of coordinates/location of the atoms. (essential)
        Occupancy (float): The value of the occupancy.
        bfactor (float): The B-factor/Temperature Factor/ Debye Waller Factor of the atom. (For NMR structures, it is RMSF accross the frames)
        Element (str): Element information of the given 'Atom' object. (Element from the Periodic Table)
        Charge (float): Charnge of the given 'Atom' object.
        parent (packman.molecule.Residue): The Residue Object (parent) this atom belongs to.
        AlternateLocationIndicator (bool): If the alternate location available for the atom (To be removed in future)

    """
    def __init__(self,id,AtomName,Coordinates,Occupancy,bfactor,Element,Charge,parent):
        self.__id=id
        self.__AtomName=AtomName
        self.__AlternateLocationIndicator=None #remove it later
        self.__parent=parent
        self.__Coordinates=Coordinates
        self.__Occupancy=Occupancy
        self.__bfactor=bfactor
        self.__SegmentIdentifier=None
        self.__Element=Element
        self.__Charge=Charge

        #Properties are the entities that are not included in the PDB files and are obtained by calculations
        self.__properties = {}
        self.__Bonds = []

    #Get Functions
    def get_id(self):
        """Get the ID of the 'Atom'

        Returns:
            int if successful, None otherwise.
        """
        return self.__id
    
    def get_name(self):
        """Get the Name of the 'Atom'

        Returns:
            str if successful, None otherwise.
        """
        return self.__AtomName
    
    def get_alternatelocationindicator(self):
        """Know if the given 'Atom' has Alternate Location Information

        Note:
            This feature might be removed in the future.

        Returns:
            True if successful, False otherwise.
        """
        return self.__AlternateLocationIndicator

    def get_parent(self):
        """Get the 'Residue' the given 'Atom' belongs to.

        Returns:
            packman.molecule.Residue if successful, None otherwise.
        """
        return self.__parent
    
    def get_location(self):
        """Get the Coordinates/Location of the given 'Atom'

        Returns:
            numpy.array if successful, None otherwise.
        """
        return self.__Coordinates
    
    def get_occupancy(self):
        """Get the Occupancy of the given 'Atom'

        Returns:
            float if successful, None otherwise.
        """
        return self.__Occupancy
    
    def get_bfactor(self):
        """Get the B-factor/Temperature Factor/ Debye Waller Factor of the given 'Atom'

        Returns:
            float if successful, None otherwise.
        """
        return self.__bfactor
    
    def get_segmentidentifier(self):
        """Get the Segment Identifier of the given 'Atom'

        Returns:
            str if successful, None otherwise.
        """
        return self.__SegmentIdentifier
    
    def get_element(self):
        """Get the Element of the given 'Atom'

        Returns:
            str if successful, None otherwise.
        """
        return self.__Element
    
    def get_charge(self):
        """Get the Charge of the given 'Atom'

        Returns:
            str if successful, None otherwise.
        """
        return self.__Charge
    
    def get_domain_id(self):
        """Get the Domain Identifier of the given 'Atom'. Hinge Prediction is Necessary for this option.

        Note:
            Domain Identifiers (Obtained after running PACKMAN):
            'FL': Flexible Linkers (Hinges); Followed by the domain number
            'DM': Domain; Followed by the domain number

        Returns:
            str if successful, None otherwise.
        """
        return self.get_parent().get_domain_id()
    
    def get_property(self,property_name):
        """Get the Property of the given 'Atom'.

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
    
    def get_bonds(self):
        """Get the bonds for the given 'Atom'

        Returns:
            list of packman.molecule.Bond if successful; [] otherwise.
        """
        return self.__Bonds
    
    def get_bond(self,atom2):
        """Get the specific bond with the specific atom.

        Args:
            atom2 (packman.molecule.Atom/ int) : 'Atom' object or the atom id.

        Returns:
            Bond (packman.molecule.Bond) if successful, None if bond does not exist.
        """
        if(type(atom2) == type(self)):
            None
        elif(type(atom2) ==int):
            Model = self.get_parent().get_parent().get_parent()
            atom2 = Model.get_atom(atom2)
        else:
            logging.warning('The atom2 parameter should either be an integer (atom id) or a packman.molecule.Atom object.')

        return list(set(self.get_bonds()).intersection( atom2.get_bonds() ))[0]

    #Set Functions
    def set_id(self,new_id):
        """Set the ID of the given 'Atom'

        Args:
            new_id (int): The ID User wishes to assign to the given 'Atom'
        """
        self.__id=new_id
    
    def set_name(self,new_name):
        """Set the Name of the given 'Atom'

        Args:
            new_name (str): The Name User wishes to assign to the given 'Atom'
        """
        self.__AtomName=new_name    
    
    def set_alternatelocationindicator(self,new_alternatelocationindicator):
        """Set the Alternate Location Indicator of the given 'Atom'

        Args:
            new_alternatelocationindicator: The Alternate Location Indicator User wishes to assign to the given 'Atom'
        """
        self.__AlternateLocationIndicator=new_alternatelocationindicator

    def set_parent(self,new_parent):
        """Set the Parent 'Residue' of the given 'Atom'

        Args:
            new_parent (packman.molecule.residue.Residue): The parent 'Residue' User wishes to assign to the given 'Atom'
        """
        self.__parent=new_parent
    
    def set_location(self,new_location):
        """Set the Coordinates/Location of the given 'Atom'

        Args:
            new_location (numpy.array): The Coordinates/Location User wishes to assign to the given 'Atom'
        """
        self.__Coordinates=new_location
    
    def set_occupancy(self,new_occupancy):
        """Set the Occupancy of the given 'Atom'

        Args:
            new_occupancy: The Occupancy User wishes to assign to the given 'Atom'
        """
        self.__Occupancy=new_occupancy
    
    def set_bfactor(self,new_bfactor):
        """Set the B-factor/Temperature Factor/ Debye Waller Factor of the given 'Atom'

        Args:
            new_bfactor (float): The B-factor/Temperature Factor/ Debye Waller Factor User wishes to assign to the given 'Atom'
        """
        self.__bfactor=new_bfactor
    
    def set_segmentidentifier(self,new_segmentidentifier):
        """Set the Segment Identifier of the given 'Atom'

        Args:
            new_segmentidentifier (str): The Segment Identifier User wishes to assign to the given 'Atom'
        """
        self.__SegmentIdentifier=new_segmentidentifier
    
    def set_elment(self,new_element):
        """Set the Element of the given 'Atom'

        Args:
            new_element (str): The Element User wishes to assign to the given 'Atom'
        """
        self.__Element=new_element
    
    def set_charge(self,new_charge):
        """Set the Charge of the given 'Atom'

        Args:
            new_charge (str): The Charge User wishes to assign to the given 'Atom'
        """
        self.__Charge=new_charge
    
    def set_property(self,property_name,value):
        """Set the Property of the given 'Atom'.

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

    def set_bond(self,new_bond):
        """Set the atom to the given 'Atom'

        Args:
            new_bond (packman.molecule.Bond) : Add new bond to the given 'Atom'
        
        Note:
            - Yet to add the functionality to delete the specific bonds.
        """
        self.__Bonds.append(new_bond)

    #Calculation Functions
    def calculate_distance(self,another_atom):
        """Calculate the Distance between the given 'Atom' and another 'Atom'

        Args:
            another_atom (packman.molecule.atom.Atom): The 'Atom' User wishes to calculate distance from.
        
        Returns:
            float if successful, None otherwise.
        """
        return numpy.linalg.norm(self.__Coordinates-another_atom.get_location())
