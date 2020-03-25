'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

import numpy

class HetAtom():
    """
    This class contains the information about the HetAtom.
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
    def get_id(self):
        """
        :returns: ID of an 'Atom' object
        """
        return self.__id
    
    def get_name(self):
        """
        :returns: Name of an 'Atom' object
        """
        return self.__AtomName
    
    def get_alternatelocationindicator(self):
        """
        :returns: Alternate Location Indicator
        """
        return self.__AlternateLocationIndicator

    def get_parent(self):
        """
        :returns: Parent residue of an 'Atom' Object. (If available)
        """
        return self.__parent
    
    def get_location(self):
        """
        :returns: Cartesian Coordinates of an 'Atom'
        """
        return self.__Coordinates
    
    def get_occupancy(self):
        """
        :returns: Occupancy of an 'Atom'
        """
        return self.__Occupancy
    
    def get_bfactor(self):
        """
        :returns: B-Factor value of an 'Atom'
        """
        return self.__bfactor
    
    def get_segmentidentifier(self):
        """
        :returns: Segment Identifier of an 'Atom'
        """
        return self.__SegmentIdentifier
    
    def get_elment(self):
        """
        :returns: Element of an 'Atom'
        """
        return self.__Element

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
