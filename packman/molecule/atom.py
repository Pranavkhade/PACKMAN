'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

import numpy


class Atom():
    """This class contains the information about the Atom.

    This class contains all the information available about the Atom. The Atom class is the lowest in the hierarchy of the 'molecule' API classes.
    the order of hierarchy being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
    Please read the Tutorials and Documentation for more details.

    Properties created with the ``@property`` decorator should be documented
    in the property's getter method.

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

    #Get Functions
    @property
    def get_id(self):
        """
        :returns: ID of an 'packman.molecule.atom.Atom' object
        """
        return self.__id
    
    @property
    def get_name(self):
        """
        :returns: Name of an 'Atom' object
        """
        return self.__AtomName
    
    @property
    def get_alternatelocationindicator(self):
        """
        :returns: Alternate Location Indicator
        """
        return self.__AlternateLocationIndicator

    @property
    def get_parent(self):
        """
        :returns: Parent residue of an 'Atom' Object. (If available)
        """
        return self.__parent
    
    @property
    def get_location(self):
        """
        :returns: Cartesian Coordinates of an 'Atom'
        """
        return self.__Coordinates
    
    @property
    def get_occupancy(self):
        """
        :returns: Occupancy of an 'Atom'
        """
        return self.__Occupancy
    
    @property
    def get_bfactor(self):
        """
        :returns: B-Factor value of an 'Atom'
        """
        return self.__bfactor
    
    @property
    def get_segmentidentifier(self):
        """
        :returns: Segment Identifier of an 'Atom'
        """
        return self.__SegmentIdentifier
    
    @property
    def get_element(self):
        """
        :returns: Element of an 'Atom'
        """
        return self.__Element
    
    @property
    def get_domain_id(self):
        """
        :returns: Domain ID of the parent Residue
        """
        return self.get_parent().get_domain_id()

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
        self.__AtomName=new_name
    
    def set_alternatelocationindicator(self,new_alternatelocationindicator):
        """
        :param new_alternatelocationindicator: New Alternate Location Indicator to be assigned
        """
        self.__AlternateLocationIndicator=new_alternatelocationindicator

    def set_parent(self,new_parent):
        """
        :param new_parent: New parent to be assigned
        """
        self.__parent=new_parent
    
    def set_location(self,new_location):
        """
        :param new_location: New Cartesian Coordinates to be assigned
        """
        self.__Coordinates=new_location
    
    def set_occupancy(self,new_occupancy):
        """
        :param new_occupancy: New occupancy value to be assigned
        """
        self.__Occupancy=new_occupancy
    
    def set_bfactor(self,new_bfactor):
        """
        :param new_bfactor: New B-Factor to be assigned
        """
        self.__bfactor=new_bfactor
    
    def set_segmentidentifier(self,new_segmentidentifier):
        """
        :param new_segmentidentifier: New Segment Identifier to be assigned
        """
        self.__SegmentIdentifier=new_segmentidentifier
    
    def set_elment(self,new_element):
        """
        :param new_element: New element to be assigned
        """
        self.__Element=new_element

    #Calculation Functions
    def calc_dist(self,another_atom):
        """
        :param another_atom: 'Atom' object to which the distance should be calculated.
        """
        return numpy.linalg.norm(self.__Coordinates-another_atom.get_location())
