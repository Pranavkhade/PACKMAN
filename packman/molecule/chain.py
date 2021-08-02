# -*- coding: utf-8 -*-
"""The 'Chain' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Chain' object documentation for details.

Example::

    from packman.molecule import Chain
    help( Chain )

Note:
    * Iterating over chain only fetches the residues not the hetero atoms
    * The object fetches the atoms and residues which are not in order as they appear in the PDB. Find a way to fix this
    
Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""
#Bond information
aa_connectivity= {
'ALA': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'HB1', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('OXT', 'HXT', 'SING')],
'ARG': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'CD', 'SING'), ('CG', 'HG2', 'SING'), ('CG', 'HG3', 'SING'), ('CD', 'NE', 'SING'), ('CD', 'HD2', 'SING'), ('CD', 'HD3', 'SING'), ('NE', 'CZ', 'SING'), ('NE', 'HE', 'SING'), ('CZ', 'NH1', 'SING'), ('CZ', 'NH2', 'DOUB'), ('NH1', 'HH11', 'SING'), ('NH1', 'HH12', 'SING'), ('NH2', 'HH21', 'SING'), ('NH2', 'HH22', 'SING'), ('OXT', 'HXT', 'SING')],
'ASN': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'OD1', 'DOUB'), ('CG', 'ND2', 'SING'), ('ND2', 'HD21', 'SING'), ('ND2', 'HD22', 'SING'), ('OXT', 'HXT', 'SING')],
'ASP': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'OD1', 'DOUB'), ('CG', 'OD2', 'SING'), ('OD2', 'HD2', 'SING'), ('OXT', 'HXT', 'SING'), ('"SYSTEMATIC', 'NAME"', 'ACDLabs')],
'CYS': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'SG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('SG', 'HG', 'SING'), ('OXT', 'HXT', 'SING')],
'GLN': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'CD', 'SING'), ('CG', 'HG2', 'SING'), ('CG', 'HG3', 'SING'), ('CD', 'OE1', 'DOUB'), ('CD', 'NE2', 'SING'), ('NE2', 'HE21', 'SING'), ('NE2', 'HE22', 'SING'), ('OXT', 'HXT', 'SING')],
'GLU': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'CD', 'SING'), ('CG', 'HG2', 'SING'), ('CG', 'HG3', 'SING'), ('CD', 'OE1', 'DOUB'), ('CD', 'OE2', 'SING'), ('OE2', 'HE2', 'SING'), ('OXT', 'HXT', 'SING'), ('"SYSTEMATIC', 'NAME"', 'ACDLabs')],
'GLY': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'HA2', 'SING'), ('CA', 'HA3', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('OXT', 'HXT', 'SING')],
'HIS': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'ND1', 'SING'), ('CG', 'CD2', 'DOUB'), ('ND1', 'CE1', 'DOUB'), ('ND1', 'HD1', 'SING'), ('CD2', 'NE2', 'SING'), ('CD2', 'HD2', 'SING'), ('CE1', 'NE2', 'SING'), ('CE1', 'HE1', 'SING'), ('NE2', 'HE2', 'SING'), ('OXT', 'HXT', 'SING')],
'ILE': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG1', 'SING'), ('CB', 'CG2', 'SING'), ('CB', 'HB', 'SING'), ('CG1', 'CD1', 'SING'), ('CG1', 'HG12', 'SING'), ('CG1', 'HG13', 'SING'), ('CG2', 'HG21', 'SING'), ('CG2', 'HG22', 'SING'), ('CG2', 'HG23', 'SING'), ('CD1', 'HD11', 'SING'), ('CD1', 'HD12', 'SING'), ('CD1', 'HD13', 'SING'), ('OXT', 'HXT', 'SING')],
'LEU': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'CD1', 'SING'), ('CG', 'CD2', 'SING'), ('CG', 'HG', 'SING'), ('CD1', 'HD11', 'SING'), ('CD1', 'HD12', 'SING'), ('CD1', 'HD13', 'SING'), ('CD2', 'HD21', 'SING'), ('CD2', 'HD22', 'SING'), ('CD2', 'HD23', 'SING'), ('OXT', 'HXT', 'SING')],
'LYS': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'CD', 'SING'), ('CG', 'HG2', 'SING'), ('CG', 'HG3', 'SING'), ('CD', 'CE', 'SING'), ('CD', 'HD2', 'SING'), ('CD', 'HD3', 'SING'), ('CE', 'NZ', 'SING'), ('CE', 'HE2', 'SING'), ('CE', 'HE3', 'SING'), ('NZ', 'HZ1', 'SING'), ('NZ', 'HZ2', 'SING'), ('NZ', 'HZ3', 'SING'), ('OXT', 'HXT', 'SING')],
'MET': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'SD', 'SING'), ('CG', 'HG2', 'SING'), ('CG', 'HG3', 'SING'), ('SD', 'CE', 'SING'), ('CE', 'HE1', 'SING'), ('CE', 'HE2', 'SING'), ('CE', 'HE3', 'SING'), ('OXT', 'HXT', 'SING')],
'PHE': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'CD1', 'DOUB'), ('CG', 'CD2', 'SING'), ('CD1', 'CE1', 'SING'), ('CD1', 'HD1', 'SING'), ('CD2', 'CE2', 'DOUB'), ('CD2', 'HD2', 'SING'), ('CE1', 'CZ', 'DOUB'), ('CE1', 'HE1', 'SING'), ('CE2', 'CZ', 'SING'), ('CE2', 'HE2', 'SING'), ('CZ', 'HZ', 'SING'), ('OXT', 'HXT', 'SING')],
'PRO': [('N', 'CA', 'SING'), ('N', 'CD', 'SING'), ('N', 'H', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'CD', 'SING'), ('CG', 'HG2', 'SING'), ('CG', 'HG3', 'SING'), ('CD', 'HD2', 'SING'), ('CD', 'HD3', 'SING'), ('OXT', 'HXT', 'SING')],
'SER': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'OG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('OG', 'HG', 'SING'), ('OXT', 'HXT', 'SING')],
'THR': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'OG1', 'SING'), ('CB', 'CG2', 'SING'), ('CB', 'HB', 'SING'), ('OG1', 'HG1', 'SING'), ('CG2', 'HG21', 'SING'), ('CG2', 'HG22', 'SING'), ('CG2', 'HG23', 'SING'), ('OXT', 'HXT', 'SING')],
'TRP': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'CD1', 'DOUB'), ('CG', 'CD2', 'SING'), ('CD1', 'NE1', 'SING'), ('CD1', 'HD1', 'SING'), ('CD2', 'CE2', 'DOUB'), ('CD2', 'CE3', 'SING'), ('NE1', 'CE2', 'SING'), ('NE1', 'HE1', 'SING'), ('CE2', 'CZ2', 'SING'), ('CE3', 'CZ3', 'DOUB'), ('CE3', 'HE3', 'SING'), ('CZ2', 'CH2', 'DOUB'), ('CZ2', 'HZ2', 'SING'), ('CZ3', 'CH2', 'SING'), ('CZ3', 'HZ3', 'SING'), ('CH2', 'HH2', 'SING'), ('OXT', 'HXT', 'SING')],
'TYR': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG', 'SING'), ('CB', 'HB2', 'SING'), ('CB', 'HB3', 'SING'), ('CG', 'CD1', 'DOUB'), ('CG', 'CD2', 'SING'), ('CD1', 'CE1', 'SING'), ('CD1', 'HD1', 'SING'), ('CD2', 'CE2', 'DOUB'), ('CD2', 'HD2', 'SING'), ('CE1', 'CZ', 'DOUB'), ('CE1', 'HE1', 'SING'), ('CE2', 'CZ', 'SING'), ('CE2', 'HE2', 'SING'), ('CZ', 'OH', 'SING'), ('OH', 'HH', 'SING'), ('OXT', 'HXT', 'SING')],
'VAL': [('N', 'CA', 'SING'), ('N', 'H', 'SING'), ('N', 'H2', 'SING'), ('CA', 'C', 'SING'), ('CA', 'CB', 'SING'), ('CA', 'HA', 'SING'), ('C', 'O', 'DOUB'), ('C', 'OXT', 'SING'), ('CB', 'CG1', 'SING'), ('CB', 'CG2', 'SING'), ('CB', 'HB', 'SING'), ('CG1', 'HG11', 'SING'), ('CG1', 'HG12', 'SING'), ('CG1', 'HG13', 'SING'), ('CG2', 'HG21', 'SING'), ('CG2', 'HG22', 'SING'), ('CG2', 'HG23', 'SING'), ('OXT', 'HXT', 'SING')],
}

import numpy
from networkx import Graph
import logging

from .bond import Bond


class Chain():
    """This class contains the information about the 'Chain' object (packman.molecule.Chain).

    This class contains all the information available about the Chain and stores the corresponding 'Residue' and 'Hetmol' objects in itself. The Chain class is the third lowest in the hierarchy of the 'molecule' API classes.
    the order of hierarchy being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
    Please read the Tutorials and Documentation for more details.

    Note:
        * Please refer to the [https://web.archive.org/web/20080905024351/http://www.wwpdb.org/docs.html] for the description of the arguments.
        * Add get_atoms()
    
    Args:
        id (str): Chain ID from the PDB file as it is. Each Chain in a PDB file Model/Frame has unique ID. (essential)
        
    """
    
    def __init__(self,id):
        self.__id = id
        self.__Residues = None
        self.__HetMols = None
        self.__parent = None
        #More Features
        self.__Hinges = []

        #Properties are the entities that are not included in the PDB files and are obtained by calculations
        self.__properties = {}
    
    def __setitem__(self,Number,Entity,Type):
        if(Type=='Residue'):
            try:
                self.__Residues[Number]=Entity
            except:
                self.__Residues={}
                self.__Residues[Number]=Entity
        elif(Type=='HetMol'):
            try:
                self.__HetMols[Number]=Entity
            except:
                self.__HetMols={}
                self.__HetMols[Number]=Entity

    def __getitem__(self,Number,Type=None):
        try:
            return self.__Residues[Number]
        except:
            return self.__HetMols[Number]
    
    #Get Functions
    def get_id(self):
        """Get the ID of the 'Residue'

        Returns:
            str if successful, None otherwise.
        """
        return self.__id
    
    def get_parent(self):
        """Get the 'Model' the given 'Chain' belongs to.

        Returns:
            packman.molecule.Model if successful, None otherwise.
        """
        return self.__parent
    
    def get_hinges(self):
        """Get the hinges in the chain as a dictionary of 'Hinge' objects.

        Returns:
            packman.molecule.annotations.Hinge if successful, None otherwise.
        """
        return self.__Hinges
    
    def get_entropy(self,entropy_type):
        """Get the Packing Entropy of the given 'Chain'.

        Please note that if the Entropy is calculated using specific atoms, this option might not give results for the amino acids that are not included because of the specific selection. Use the get_total_chain_entropy() function. Please see the documentation for more details.

        Args:
            type (str): Type of entropy (Allowed Values: 1. PackingEntropy)

        Note:
            - More type of Entropies might be added in the future.
        """
        EntropyTypes = ['PackingEntropy']
        try:
            return numpy.sum( [self.__Residues[i].get_entropy(entropy_type) for i in self.__Residues] )
        except:
            if(entropy_type in EntropyTypes):
                logging.warning('This Entropy type is not yet calculated for parameters provided. Please check the "calculate_entropy" function in the documentation for the details.')
            else:
                logging.warning('The Entropy type provided is invalid. Please check the documentation for the details.')

    def get_property(self,property_name):
        """Get the Property of the given 'Chain'.

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
        """Get the Bonds in the given 'Chain'

        Please check out calculate_bonds() for more information on the bond calculations.
        """
        try:
            yield self.__AllBonds
        except:
            logging.warn('Bonds are not calculated yet. Use calculate_bond()')
    
    #Set Functions
    def set_id(self,new_id):
        """Set the ID of the given 'Chain'

        Args:
            new_id (str): The ID User wishes to assign to the given 'Chain'
        """
        self.__id=new_id
    
    def set_parent(self,parent):
        """Set the Parent of the given 'Residue'

        Args:
            new_parent (packman.molecule.Chain): The parent 'Chain' User wishes to assign to the given 'Residue'
        """
        self.__parent=parent
    
    def set_hinges(self,new_hinges):
        """Set/Add hinge to the 'Chain' object

        Args:
            new_hinges (packman.molecule.annotations.Hinge): The 'Hinge' User wishes to assign/add to the given 'Chain'
        """
        self.__Hinges=self.__Hinges+new_hinges

    def set_property(self,property_name,value):
        """Set the Property of the given 'Chain'.

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
    def get_atoms(self):
        """Get the generator of corresponding 'atom' objects of the 'Chain'

        Returns:
            generator of 'atom' objects if successful, None otherwise.
        
        Note:
            find a way to deal with hetatoms
        """
        for i in sorted(self.__Residues.keys()):
            for j in self.__Residues[i].get_atoms():
                yield j

    def get_residues(self):
        """Get the generator of corresponding 'Residue' objects of the 'Chain'

        Returns:
            generator of 'Residue' objects if successful, None otherwise.
        """
        for i in sorted(self.__Residues.keys()):yield self.__Residues[i]
    
    def get_calpha(self):
        """Get the C-Alpha atoms of the 'Chain' as an 'Atom' object.

        Returns:
            generator of packman.molecule.Atom objects if successful, None otherwise.
        """
        for i in self.get_residues():yield i.get_calpha()
    
    def get_hetmols(self):
        """Get the generator of corresponding 'HetMol' objects of the 'Chain'

        Returns:
            generator of 'HetAtom' objects if successful, None otherwise.
        """
        for i in sorted(self.__HetMols.keys()):yield self.__HetMols[i]
    
    def get_backbone(self):
        """Get the Backbone atoms of the given 'Chain' as a list of 'Atom' object

        Note:
            * Backbone Atoms: CA, O, N, C

        Returns:
            list of packman.molecule.Atom if successful, None otherwise.
        """
        return [i.get_backbone() for i in self.get_residues()]
    
    def calculate_entropy(self,entropy_type,chains=None, probe_size=1.4, onspherepoints=30):
        """Calculate the entropy for the each amino acid will be returned.
    
        The 'chains' argument should be used when the user wants to restrict the analysis to a chain or group of chains rather than the whole structure.

        Args:
            entropy_type (str)              : Type of entropy to be calculated (Options: PackingEntropy)
            chains ([str]/str)              : Chain IDs for the Entropy calculation (None means all the chains are included; single string means only one chain ID; multiple chains should be an array of strings).
            probe_size (float)              : Radius of the probe to generate the surface points (This value should not be less than 1;Read the Publication for more details)
            onspherepoints (int)            : Number of points to be generated around each point for the surface (Read the Publication for more details)
        """
        self.get_parent().calculate_entropy(entropy_type,chains=chains, probe_size=probe_size, onspherepoints=onspherepoints)
    
    def calculate_bonds(self):
        """Calculate the bonds in the given 'Chain'.

        Currently, bonds are only calculated based on the following RCSB PDB resource file:
        http://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif.gz
        """
        resi = [i for i in self.get_residues()]
        counter = 0
        self.__AllBonds = []
        self.__ChainGraph = Graph()

        for i in resi:
            try:
                for j in aa_connectivity[i.get_name()]:
                    try:
                        atom1, atom2 =  i.get_atom(j[0]), i.get_atom(j[1])
                        if(atom1!=None and atom2!=None and atom1.get_parent().get_parent().get_id() == atom2.get_parent().get_parent().get_id() ):
                            counter += 1
                            if(j[2]=='SING'):
                                bond = Bond(counter, atom1, atom2, 'covalent-single', source='RCSB/aa-variants-v1.cif')
                            elif(j[2]=='DOUB'):
                                bond = Bond(counter, atom1, atom2, 'covalent-double', source='RCSB/aa-variants-v1.cif')
                            else:
                                bond = Bond(counter, atom1, atom2, 'covalent', source='RCSB/aa-variants-v1.cif')
                            self.__AllBonds.append(bond)
                            atom1.set_bond(bond)
                            atom2.set_bond(bond)
                            self.__ChainGraph.add_node( atom1.get_id() )
                            self.__ChainGraph.add_node( atom2.get_id() )
                            self.__ChainGraph.add_edge( atom1.get_id(), atom2.get_id() , id = counter )
                    except Exception as e:
                        #Bond pair not found
                        None
            except:
                logging.warn('Residue Number|Name|Chain '+str(i.get_id())+'|'+str(i.get_name())+'|'+str(i.get_parent().get_id())+' is not a standard amino acid.')

        #Peptide bonds (Assumption: All the reisudes are in the incremental/decremental order)
        for i in range(0,len(resi)-1):
            if( abs(resi[i].get_id() - resi[i+1].get_id()) == 1 ):
                counter+=1
                C, N = resi[i].get_atom('C'), resi[i+1].get_atom('N')
                bond = Bond(counter, C, N, 'covalent-single', source='PACKMAN Peptide Bond calculation' )
                self.__AllBonds.append(bond)
                C.set_bond(bond)
                N.set_bond(bond)
                self.__ChainGraph.add_edge( C.get_id(), N.get_id() , id = counter )
            else:
                logging.info('The peptide bond following residue is missing: '+str(resi[i].get_id())+' Chain: '+resi[i].get_parent().get_id())
