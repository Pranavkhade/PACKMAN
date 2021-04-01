# -*- coding: utf-8 -*-
"""The various annotation objects host file.

This is file information, not the class information. This information is only for the API developers.
Please read the corrosponding object documentation for details.

Example::

    from packman.molecule import Hinge
    help( Hinge )

Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""


'''
##################################################################################################
#                                 Structural Super Elements                                      #
##################################################################################################
'''

class Hinge():
    """This class contains the information about the 'Hinge' object (packman.molecule.Hinge).

    This class contains all the information available about the Hinge. The Hinge class is the not in the hierarchy of the 'molecule' API classes.
    However, it resides in the packman.molecule.Chain object because there are unique hinges per chain.
    the order of hierarchy of Chain being being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
    Please read the Tutorials and Documentation for more details.

    Args:
        hid (int)                          : Unique Hinge ID
        alpha_value (float)                : The hinge prediction algorithm parameter (Alpha Value) of the Hinge.
        elements (packman.molecule.Residue): The elements defining the current hinge (Currently residue objects)
        stats ([float])                    : Everything about the statistics (mean/median/mode of B-factors)
        p (float)                          : p-value obtained from the permutation test (Please read the paper for more details)
    """
    def __init__(self,hid,alpha_value,elements,stats,p):
        self.__id=hid
        self.__alpha_value=alpha_value
        self.__elements=elements
        self.__stats=stats
        self.__p=p

    #Get functions
    def get_id(self):
        """Get the ID of the 'Hinge'

        Returns:
            int if successful, None otherwise.
        """
        return self.__id
    
    def get_alpha_value(self):
        """Get the hinge prediction algorithm parameter (Alpha Value) of the 'Hinge'.

        Returns:
            float if successful, None otherwise
        """
        return self.__alpha_value

    def get_elements(self):
        """Get the elements (Currently residue objects) of the 'Hinge'

        Returns:
            list of elements (Currently residue objects) if successful, None otherwise.
        """
        return self.__elements

    def get_stats(self):
        """Get the statistics of the 'Hinge'

        Returns:
            Statistics of the Hinge such as mean/median/mode of B-factors if successful, None otherwise.
        """
        return self.__stats

    def get_pvalue(self):
        """Get the statistics of the 'Hinge'

        Returns:
            p-value (float) : p-value obtained from the permutation test (Please read the paper for more details)
        """
        return self.__p
