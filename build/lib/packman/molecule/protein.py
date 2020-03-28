# -*- coding: utf-8 -*-
"""The 'Protein' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Protein' object documentation for details. [ help(packman.molecule.protein.Protein) ]

Example:
    >>>from packman.molecule import protein.Protein
    >>>help( Protein )
    OR

    >>>from packman import molecule
    >>>help( molecule.Protein )

Note:
    * Top in the hierarchy
    
Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""

from numpy import around

class Protein():
    """This class contains the information about the 'Protein' object (packman.molecule.Protein).

    This class contains all the information available about the Protein and stores everything in itself. The Protein class is the highest in the hierarchy of the 'molecule' API classes.
    the order of hierarchy being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
    Please read the Tutorials and Documentation for more details.

    Note:
        * Please refer to the [https://web.archive.org/web/20080905024351/http://www.wwpdb.org/docs.html] for the description of the arguments.
    

    Args:
        id (int)                         : Protein ID
        name (str)                       : Protein Name
        Models ([packman.molecule.Model]): Protein models/frames of the structure. (NMR files usually have multiple conformers of the same protein)
    """
    
    def __init__(self,id,name,Models):
        self.id=id
        self.name=name
        self.Models=Models
    
    def __getitem__(self,ModelNumber):
        return self.Models[ModelNumber]
    
    def write_pdb(self,filename):
        """[summary]
        
        Args:
            filename (str): Name of the output file user wishes to assign.
        """
        open(filename,'w').write('')
        fh=open(filename,'a')
        for num_,_ in enumerate(self):
            fh.write("Model\t"+str(num_)+'\n')
            for i in _.get_atoms():
                fh.write("ATOM  %5s %-4s %3s %1s%4s    %8s%8s%8s%6s%6s         %-4s%2s%2s\n"%(i.get_id(),i.get_name(),i.get_parent().get_name(),i.get_parent().get_parent().get_id(),i.get_parent().get_id(),around(i.get_location()[0],decimals=3),around(i.get_location()[1],decimals=3),around(i.get_location()[2],decimals=3),i.get_occupancy(),i.get_bfactor(),'',i.get_element(),''))
            fh.write("ENDMDL\n")
        return True