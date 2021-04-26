# -*- coding: utf-8 -*-
"""The 'Protein' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Protein' object documentation for details.

Example::

    from packman.molecule import Protein
    help( Protein )

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
    
    def __init__(self,id,Models):
        self.__id=id
        self.__Models=Models
        self.__Data = None
    
    #Get Functions
    def __getitem__(self,ModelNumber):
        return self.__Models[ModelNumber]
    
    def get_id(self):
        """Get the ID for the Protein object.

        Returns:
            String if successful, None otherwise.
        """
        return self.__id

    def get_data(self):
        """Get the misc data (other than coordinates) from the file.

        Returns:
            Array of Strings
        """
        return self.__Data

    #Set functions
    def set_data(self,data):
        """Set the misc data (other than coordiantes) to the Protein object.

        Args:
            data (array): Array of String 
        
        Note:
            - All the properties are planned to be put in specific format to achieve complete interformat conversion.
        """
        self.__Data = data


    #Wite Functions
    def write_pdb(self,filename):
        """Write a PDB (.pdb) file from the Protein object.
        
        Args:
            filename (str): Name of the output file user wishes to assign.
        """
        open(filename,'w').write('')
        fh=open(filename,'a')
        #Annotations
        try:
            for i in self.get_data():
                fh.write(i+'\n')
        except:
            None
        #Atoms part
        for num_,_ in enumerate(self):
            fh.write("Model\t"+str(num_)+'\n')
            for i in _.get_atoms():
                fh.write("ATOM  %5s %-4s %3s %1s%4s    %8s%8s%8s%6s%6s         %-4s%2s%2s\n"%(i.get_id(), i.get_name(), i.get_parent().get_name(), i.get_parent().get_parent().get_id(), i.get_parent().get_id(), around(i.get_location()[0],decimals=3).real, around(i.get_location()[1],decimals=3).real, around(i.get_location()[2],decimals=3).real, i.get_occupancy(), round(i.get_bfactor(),2), '', i.get_element(), ''))
            fh.write("ENDMDL\n")
        fh.flush()
        fh.close()
        return True
    
    def write_cif(self,filename):
        """Write a PDBx/mmCIF (.cif) file from the Protein object.
        
        Args:
            filename (str): Name of the output file user wishes to assign.
        """
        #Annotations
        open(filename,'w').write('')
        fh=open(filename,'a')
        if(self.__Data != None):
            for i in self.__Data:
                try:
                    fh.write(i+'\n')
                except:
                    None
        fh.write('#\nloop_\n')
        #'_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol', '_atom_site.label_atom_id', '_atom_site.label_alt_id', '_atom_site.label_comp_id', '_atom_site.label_asym_id', '_atom_site.label_entity_id', '_atom_site.label_seq_id', '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy', '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge', '_atom_site.auth_seq_id', '_atom_site.auth_comp_id', '_atom_site.auth_asym_id', '_atom_site.auth_atom_id', '_atom_site.pdbx_PDB_model_num'
        cols = ['_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol', '_atom_site.label_atom_id', '_atom_site.label_alt_id', '_atom_site.label_comp_id', '_atom_site.label_asym_id', '_atom_site.label_entity_id', '_atom_site.label_seq_id', '_atom_site.pdbx_PDB_ins_code', '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z', '_atom_site.occupancy', '_atom_site.B_iso_or_equiv', '_atom_site.pdbx_formal_charge', '_atom_site.auth_seq_id', '_atom_site.auth_comp_id', '_atom_site.auth_asym_id', '_atom_site.auth_atom_id', '_atom_site.pdbx_PDB_model_num']

        for i in cols:
            fh.write(i+'\n')
        
        for num_,_ in enumerate(self):
            for i in _.get_atoms():
                fh.write("ATOM\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%( i.get_id(), i.get_element(), i.get_name(), '.' , i.get_parent().get_name(), i.get_parent().get_parent().get_id() , '1', i.get_parent().get_id(), '?', around(i.get_location()[0],decimals=3).real, around(i.get_location()[1],decimals=3).real, around(i.get_location()[2],decimals=3).real, i.get_occupancy(), round(i.get_bfactor(),2), i.get_charge(), i.get_parent().get_id(), i.get_parent().get_name(), i.get_parent().get_parent().get_id(), i.get_name(), num_+1  ) )
        
        fh.flush()
        fh.close()

        return True
    
    def write_structure(self,filename,ftype='cif'):
        """Write the 'Protein' object to the file.
        
        CIF file format is default because it has more advantages over PDB format and PDB format is 'frozen'. Please read following for more information::
            https://www.wwpdb.org/documentation/file-formats-and-the-pdb

        Args:
            filename (str): Name of the output file user wishes to assign.
            ftype    (str): Format for the file (pdb / cif)
        """
        try:
            ftype = filename.split('.')[1]
        except:
            None
        
        if(ftype == 'cif'):
            self.write_cif(filename)
        elif(ftype == 'pdb'):
            self.write_pdb(filename)
        else:
            raise Exception('Please provide appropriate "ftype" argument. (cif/pdb).')

        return True