# -*- coding: utf-8 -*-
"""The 'Model' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Model' object documentation for details.

Citation:
    Pranav M Khade, Robert L Jernigan, PACKMAN-Molecule: Python Toolbox for Structural Bioinformatics, Bioinformatics Advances, 2022;, vbac007, https://doi.org/10.1093/bioadv/vbac007

Example::

    from packman.molecule import Model
    help( Model )

Note:
    * The models are frames of the PDB file.

Authors:
    * Pranav Khade (https://github.com/Pranavkhade)
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

from ..entropy import PackingEntropy
from .bond import Bond

import numpy
import logging

from networkx import Graph, connected_components

from typing import TYPE_CHECKING, Union, List, Iterable, Dict

if(TYPE_CHECKING):
    from .. import Atom, Residue, HetMol, Chain, Protein


class Model():
    """This class contains the information about the 'Model' object (packman.molecule.Model).

        This class contains all the information available about the Model and stores the corresponding 'Chain' objects. The Model class is the second lowest in the hierarchy of the 'molecule' API classes.
        the order of hierarchy being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
        Please read the Tutorials and Documentation for more details.

        Notes:
            * Model index starts from 0 while retriving but identifiers (ids) start from 1.
            * Important: get_hetmols and get_residues must be used based on your interest in protein/ non-protein components.
            * Please refer to the [https://web.archive.org/web/20080905024351/http://www.wwpdb.org/docs.html] for the description of the arguments.
        
        Args:
            id (int)                                      : Model ID from the PDB file ordered from first to the last. Each Model in a PDB file has unique ID. (essential)
            AllAtoms ({ref:`packman.molecule.Atom`})      : Dictionary of all the 'Atom' in the given model.
            AllResidues ({ref:`packman.molecule.Residue`}): Dictionary of all the 'Residue' in the given model.
            AllChains ({ref:`packman.molecule.Chain`})    : Dictionary of all the 'Chain' in the given model.
            AllHetAtoms ({ref:`packman.molecule.HetAtom`}): Dictionary of all the 'HetAtom' in the given model.
            AllHetMols ({ref:`packman.molecule.HetMol`})  : Dictionary of all the 'HetMol' in the given model.

        """
        
    def __init__(self, id: int, AllAtoms: Dict[int, 'Atom'], AllResidues: Dict[str, 'Residue'], AllChains: Dict[str, 'Chain'], AllHetAtoms: Dict[int, 'Atom'], AllHetMols: Dict[str, 'HetMol']):                
        self.__id               = id
        self.__AllAtoms         = AllAtoms
        self.__AllAtoms_inverse = { self.__AllAtoms[i]:i for i in self.__AllAtoms }
        self.__AllResidues      = AllResidues
        self.__AllChains        = AllChains
        self.__AllHetAtoms      = AllHetAtoms
        self.__AllHetMols       = AllHetMols
        self.__parent           = None
        
        #Properties are the entities that are not included in the PDB files and are obtained by calculations
        self.__properties = {}

    def __getitem__(self, ChainID: str) -> 'Chain':
        return self.__AllChains[ChainID]
    
    def __repr__(self) -> str:
        return '<Model '+str(self.__id)+' Index: '+str(self.__id-1)+'>'
    
    #Get Functions
    def get_id(self) -> int:
        """Get the ID of the 'Model'

        Returns:
            int if successful, None otherwise.
        """
        return self.__id

    def get_chains(self) -> Iterable['Chain']:
        """Get the list of corresponding 'Chain' objects of the 'Model'

        Returns:
            Generator of ref:`packman.molecule.Chain` if successful, None otherwise.
        """
        for i in sorted(self.__AllChains.keys()): yield self.__AllChains[i]

    def get_residues(self) -> List['Chain']:
        """Get the generator of corresponding 'Residue' objects of the 'Model'

        Returns:
            array of 'Residue' objects if successful, None otherwise.
        """
        residues = []
        for i in self.__AllChains.keys():
            try:
                residues.extend( self.__AllChains[i].get_residues() )
            except:
                logging.warning("Chain "+str(i)+" either doesn't have residues or an error occurred; Model.get_residues() may have loaded other chains.")
        return residues
    
    def get_atoms(self) -> Iterable['Atom']:
        """Get the generator of corresponding 'Atom' objects of the 'Model'

        Returns:
            generator of 'Atom' objects if successful, None otherwise.
        """
        for i in sorted(self.__AllAtoms.keys()): yield self.__AllAtoms[i]
    
    def get_atom(self, idx: int) -> 'Atom':
        """Get the atom of the given ID in the current model.

        Note: - This function is different from :py:func:`packman.molecule.Chain.get_atoms` and also :py:func:`packman.molecule.Chain.get_atom`
              - If the PDB file is constructed manually/ has multiple atoms of the same ID, the first instance of the atom with that id is returned. Please avoid saving two atoms with same ID in a same structure file in a given frame/model.

        Args:
            idx (int): Get atom by the id
        
        Returns:
            atom (:py:class:`packman.molecule.Atom`): Atom of the given ID if successful; None otherwise.
        """
        the_atom = None
        for i in self.__AllResidues:
            the_atom =  self.__AllResidues[i].get_atom(idx)
            if(the_atom is not None):
                break
            
        if(the_atom is None):
            for i in self.__AllHetMols:
                the_atom =  self.__AllHetMols[i].get_atom(idx)
                if(the_atom is not None):
                    break
        if(the_atom is not None):
            return the_atom
        else:
            logging.info('The atom with the given ID is not found in this Model')
            return None
    
    def get_chain(self, ChainID: str) -> 'Chain':
        """Get the corresponding 'Chain' object of the current Model.

        Returns:
            'Chain' object if successful, None otherwise.
        """
        return self.__AllChains[ChainID]
    
    def get_hetmols(self) -> Iterable['HetMol']:
        """Get the generator of corresponding 'HetMol' objects of the 'Model'

        Returns:
            generator of ref:`packman.molecule.HetMol` objects if successful, None otherwise.
        """
        for i in sorted(self.__AllHetMols.keys()): yield self.__AllHetMols[i]
     
    def get_hetatoms(self) -> Iterable['Atom']:
        """Get the generator of corresponding 'HetAtom' objects of the 'Model'

        Returns:
            generator of 'packman.molecule.HetAtom' objects if successful, None otherwise.
        """
        for i in sorted(self.__AllHetAtoms.keys()): yield self.__AllHetAtoms[i]
    
    def get_parent(self) -> 'Protein':
        """Get the 'Protein' parent of the 'Model' object.

        Returns:
            ref:`packman.molecule.Protein` object if successful, None otherwise.
        """
        return self.__parent
    
    def get_entropy(self, entropy_type: str) -> float:
        """Get the Packing Entropy of the given 'Chain' for the Current Model.

        Args:
            type (str): Type of entropy (Allowed Values: 1. PackingEntropy)

        Note:
            * More type of Entropies might be added in the future.
        
        Returns:
            Total Entropy (float) of the complex if successful, None otherwise.
        """
        EntropyTypes = ['PackingEntropy']
        try:
            return numpy.sum( [self.__AllResidues[i].get_entropy(entropy_type) for i in self.__AllResidues] )
        except:
            if(entropy_type in EntropyTypes):
                logging.warning('This Entropy type might not be calculated for all the chains/residues. Please check the "calculate_entropy" function in the documentation for the details.')
            else:
                logging.warning('The Entropy type provided is invalid. Please check the documentation for the details.')

    def get_property(self, property_name):
        """Get the Property of the given 'Model'.

        Property is any key and value combination that can be assigned to this object. This (along with the set_property) feature is mainly useful for the user customization.
        Properties are like pinboards. You can pin anything to the object with a key as a pin.

        Args:
            property_name (object): The 'Key' or a name the user wants to assign to to the property
        
        Note:
            - Users can add custom annotations; for example: If particular model is disordered, it can be annotated with this feature.
        """
        try:
            return self.__properties[property_name]
        except:
            logging.warning('The Property Name provided is not assigned.')
    
    def get_bonds(self) -> Iterable['Bond']:
        """Return all the bonds in the given 'Model'.

        Returns:
            generator of ref:`packman.molecule.Bond` objects if successful, None otherwise.
        """
        try:
            for i in self.__AllBonds: yield self.__AllBonds[i]
        except:
            logging.warning('Failed to return the bonds.')
            
    def get_bond(self, idx: int) -> 'Bond':
        """Return the specific bond with ID in the current Model.

        Args:
            idx (int): Get the 'Bond' by the id
        
        Returns:
            ref:`packman.molecule.Bond` object if successful, None otherwise.
        """
        try:
            return self.__AllBonds[ idx ]
        except:
            logging.warning('Failed to return the bond. Please check the ID.')

    #Compute Functions
    def get_calpha(self) -> List['Atom']:
        """Get the C-Alpha atom of the 'Model' as an 'Atom' object from the current Model.

        Returns:
            list of packman.molecule.Atom if successful, None otherwise.
        """
        return [i.get_calpha() for i in self.get_residues()]
    
    def get_backbone(self) -> List[List['Atom']]:
        """Get the Backbone atoms of the given 'Model' as a list of 'Atom' object

        Note:
            * Backbone Atoms: CA, O, N, C

        Returns:
            list of packman.molecule.Atom if successful, None otherwise.
        """
        return [i.get_backbone() for i in self.get_residues()]
    
    def get_torsion(self, bond: 'Bond', neighbor1: Union[int, 'Atom'], neighbor2: Union[int, 'Atom'], radians: bool=True) -> float:
        """Calculate the torsion angle of the given covalent bond with the corresponding selected neighbors in the current Model.

        Note:
            At least four atoms are needed to form two planes that can measure the torsional angles; therefore, along with the two bond atoms, the user needs to provide the additional two atoms that are ideally non-mutual neighbors of the atoms in the bond.

        Args:
            bond      (packman.molecule.Bond)     : The bond user wishes to calculate torsion angle to.
            neighbor1 (int/packman.molecule.Atom) : Neighbour of the Atom1 as an 'Atom' object or Atom ID.
            neighbor1 (int/packman.molecule.Atom) : Neighbour of the Atom2 as an 'Atom' object or Atom ID.
            radians   (bool)                : Return value of the angle in radians (returns value in degrees if False; Default : True)
        
        Returns:
            The torsion angle in radians/degrees if sucessful, None otherwise.
        """
        assert type(bond) == Bond, 'The bond varible should be a packman.molecule.Bond object'

        try:
            if(self.__AllBonds == {}):
                logging.warning('Please check if the Model.calculate_bonds() was executed successfully.')
                return None
        except:
            logging.warning('Please check if the Model.calculate_bonds() was executed successfully.')
            return None
        
        atom1, atom2 = bond.get_atoms()

        if( neighbor1 is None ):
            logging.error('neighbour1 not selected. Options available for neighbour1: ' +  ', '.join( [str(i) for i in self.__ModelGraph[atom1.get_id()] if i != atom2.get_id() ] ) )
            return None

        if( neighbor2 is None ):
            logging.error('neighbour2 not selected. Options available for neighbour2: ' + ', '.join( [str(i) for i in self.__ModelGraph[atom2.get_id()] if i != atom1.get_id() ]  ) )
            return None
        
        if(isinstance(neighbor1, int)):
            neighbor1 = self.__AllAtoms[neighbor1]
        elif(type(neighbor1) == type(atom1)):
            None
        else:
            logging.error('neighbour1 should either be an integer or a packman.molecule.Atom object.')
            return None
        
        if(isinstance(neighbor2, int)):
            neighbor2 = self.__AllAtoms[neighbor2]
        elif(type(neighbor2) == type(atom2)):
            None
        else:
            logging.error('neighbour2 should either be an integer or a packman.molecule.Atom object.')
            return None

        #Actual code
        b0 = -1.0*(atom1.get_location() - neighbor1.get_location())
        b1 = atom2.get_location() - atom1.get_location()
        b2 = neighbor2.get_location() - atom2.get_location()
        b1 /= numpy.linalg.norm(b1)
        v = b0 - numpy.dot(b0, b1)*b1
        w = b2 - numpy.dot(b2, b1)*b1
        x = numpy.dot(v, w)
        y = numpy.dot(numpy.cross(b1, v), w)
        radang=numpy.arctan2(y, x)
        
        if(radians):
            return radang
        else:
            return numpy.rad2deg(radang)
    
    def get_sequence(self) -> str:
        """Get the Amino acid sequence of the chain. (Protein chains only) of the current Model.
        
        Returns:
            FASTA format string of the chain sequence.
        """
        return '\n'.join( [chain.get_sequence() for chain in self.get_chains() if chain.get_sequence()!=''] )
    
        #Set Functions
    def set_id(self, new_id: int):
        """Set the ID of the given 'Model'

        Args:
            new_id (int): The ID User wishes to assign to the given 'Model'
        """
        self.__id = new_id
    
    def set_parent(self, new_parent: 'Protein'):
        """Set the 'Protein' object as a parent to the 'Model' object.

        Args:
            new_parent (ref:`packman.molecule.Protein`): The 'Protein' object as a parent to the given 'Model'
        """
        self.__parent = new_parent
    
    def set_property(self, property_name, value):
        """Set the Property of the given 'Model'.

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
    def calculate_entropy(self, entropy_type: str,chains: Union[List[str], str, None]=None, probe_size: float=1.4, onspherepoints: int=30):
        """Calculate the entropy for the each amino acid.
    
        The 'chains' argument should be used when the user wants to restrict the analysis to a chain or group of chains rather than the whole structure.

        Args:
            entropy_type (str)              : Type of entropy to be calculated (Options: PackingEntropy)
            chains ([str]/str)              : Chain IDs for the Entropy calculation (None means all the chains are included; single string means only one chain ID; multiple chains should be an array of strings).
            probe_size (float)              : Radius of the probe to generate the surface points (This value should not be less than 1;Read the Publication for more details)
            onspherepoints (int)            : Number of points to be generated around each point for the surface (Read the Publication for more details)
        """
        if(entropy_type=='PackingEntropy'):
            PackingEntropy(self.get_atoms(), chains=chains, probe_size=probe_size, onspherepoints=onspherepoints)
        else:
            logging.warning("Please provide the valid type for the Entropy Calculation.")
        
    def set_torsion(self, bond: 'Bond', theta: float, neighbor1: Union[int, 'Atom'], neighbor2: Union[int, 'Atom'], radians: bool=True) -> Union[bool, None]:
        """Set the torsion for the given covalent bond with the corresponding selected neighbors.

        Note:
            At least four atoms are needed to form two planes that change the torsional angles; therefore, along with the two bond atoms, the user needs to provide the additional two atoms that are ideally non-mutual neighbors of the atoms in the bond.

        Args:
            bond      (packman.molecule.Bond)     : The bond user wishes to rotate.
            theta     (float)                     : Set the torsional angle (see the 'radians' parameter description)
            neighbor1 (int/packman.molecule.Atom) : Neighbour of the Atom1 as an 'Atom' object or Atom ID.
            neighbor1 (int/packman.molecule.Atom) : Neighbour of the Atom2 as an 'Atom' object or Atom ID.
            radians   (bool)                      : Parameter 'theta' will be assuned to be in Radians if True, Degrees will be assumed when False. ( Default : True)
        
        Returns:
            True if successful; None otherwise
        """
        assert type(bond) == Bond, 'The bond varible should be a packman.molecule.Bond object'

        try:
            if(self.__AllBonds=={}):
                logging.warning('Please check if the Model.calculate_bonds() was executed successfully.')
                return None
        except:
            logging.warning('Please check if the Model.calculate_bonds() was executed successfully.')
            return None
        
        atom1, atom2 = bond.get_atoms()

        if( neighbor1 is None ):
            logging.error('neighbour1 not selected. Options available for neighbour1: ' +  ', '.join( [str(i) for i in self.__ModelGraph[atom1.get_id()] if i != atom2.get_id() ] ) )
            return None

        if( neighbor2 is None ):
            logging.error('neighbour2 not selected. Options available for neighbour2: ' + ', '.join( [str(i) for i in self.__ModelGraph[atom2.get_id()] if i != atom1.get_id() ]  ) )
            return None
        
        if(isinstance(neighbor1, int)):
            neighbor1 = self.__AllAtoms[neighbor1]
        elif(type(neighbor1)==type(atom1)):
            None
        else:
            logging.error('neighbour1 should either be an integer or a packman.molecule.Atom object.')
            return None
        
        if(isinstance(neighbor2, int)):
            neighbor2 = self.__AllAtoms[neighbor2]
        elif(type(neighbor2)==type(atom2)):
            None
        else:
            logging.error('neighbour2 should either be an integer or a packman.molecule.Atom object.')
            return None

        if(radians):
            None
        else:
            theta = numpy.deg2rad(theta)
        
        #Actual code
        current_torsion = self.get_torsion(bond, neighbor1=neighbor1, neighbor2=neighbor2, radians=True)

        rotang = current_torsion - theta
        
        #Find the section to rotate
        test = self.__ModelGraph.copy()
        test.remove_edge( atom1.get_id(), atom2.get_id() )
        before_components = [j for j in connected_components(self.__ModelGraph)]
        after_componets = [j for j in connected_components(test)]
        #Changed components tell you what parts have been changed
        changed_components = []
        for i in after_componets:
            if(i in before_components):
                #All the unchanged componets after breaking the bridge
                None
            else:
                changed_components.append(i)

        # Whatever component neighbour 1 is in, rotate that. This keeps $\vec{e_a}$ condition satisfied
        if( self.__AllAtoms_inverse[neighbor1] in changed_components[0] ):
            to_rotate = list(changed_components[0])
        elif( self.__AllAtoms_inverse[neighbor1] in changed_components[1] ):
            to_rotate = list(changed_components[1])
        else:
            logging.error('Atom :'+str(neighbor1.get_id())+' is not in the connected to the components')
            return None

        #Setting up rotation matrix
        sn=numpy.sin(rotang)
        cs=numpy.cos(rotang)
        t=1-cs

        v2x = atom1.get_location()[0] - atom2.get_location()[0]
        v2y = atom1.get_location()[1] - atom2.get_location()[1]
        v2z = atom1.get_location()[2] - atom2.get_location()[2]
        #Normalize the rotation vector
        mag=numpy.sqrt(numpy.square(v2x)+numpy.square(v2y)+numpy.square(v2z))
        x = float(v2x)/mag
        y = float(v2y)/mag
        z = float(v2z)/mag
        #set up the rotation matrix
        m=numpy.zeros(9)
        m[0]= t*x*x + cs
        m[1]= t*x*y + sn*z
        m[2]= t*x*z - sn*y
        m[3]= t*x*y - sn*z
        m[4]= t*y*y + cs
        m[5]= t*y*z + sn*x
        m[6]= t*x*z + sn*y
        m[7]= t*y*z - sn*x
        m[8]= t*z*z + cs

        #Rotate The Atoms
        tx = atom1.get_location()[0]
        ty = atom1.get_location()[1]
        tz = atom1.get_location()[2]

        #Set the angles
        for i in to_rotate:
            try:
                COORDS = self.__AllAtoms[i].get_location()
                COORDS[0]-=tx
                COORDS[1]-=ty
                COORDS[2]-=tz
                x=COORDS[0]*m[0]+COORDS[1]*m[1]+COORDS[2]*m[2]
                y=COORDS[0]*m[3]+COORDS[1]*m[4]+COORDS[2]*m[5]
                z=COORDS[0]*m[6]+COORDS[1]*m[7]+COORDS[2]*m[8]
                COORDS[0]=x
                COORDS[1]=y
                COORDS[2]=z
                COORDS[0]+=tx
                COORDS[1]+=ty
                COORDS[2]+=tz
                self.__AllAtoms[i].set_location( numpy.array(COORDS) )
            except:
                logging.warning('Failed to rotate Atom: '+str(self.__AllAtoms[i].get_id()))
        return True


    def calculate_bonds(self):
        """Calculate the bonds in the given 'Model'.

        Currently, bonds are only calculated based on the following RCSB PDB resource file:
        http://ftp.wwpdb.org/pub/pdb/data/monomers/aa-variants-v1.cif.gz
        """
        counter = -1
        self.__AllBonds = {}
        self.__ModelGraph = Graph()

        #PDB File CONECT Records from annotations
        for i in self.get_parent().get_data():
            if(i[0:6]=='CONECT'):
                #Information on the following line is obtained from: https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html and http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_69.html and https://mmcif.wwpdb.org/docs/pdb_to_pdbx_correspondences.html#CONECT
                other_atoms = [ (11,16,'covalent'), (16,21,'covalent'), (21,26,'covalent'), (26,31,'covalent'), (31,36,'hydrogen'), (36,41,'hydrogen'), (41,46,'salt-bridge'), (46,51,'hydrogen'), (51,56,'hydrogen'), (56,61,'salt-bridge') ]
                main_atom = self.get_atom( int(i[6:11]) )
                for j in other_atoms:
                    try:
                        temp_atom = self.get_atom(int(i[j[0]:j[1]]))
                        counter+=1
                        bond = Bond(counter, main_atom, temp_atom, j[2], source='CONECT-section')
                        self.__AllBonds[counter] = bond
                        main_atom.set_bond(bond)
                        temp_atom.set_bond(bond)
                        self.__ModelGraph.add_node( main_atom.get_id() )
                        self.__ModelGraph.add_node( temp_atom.get_id() )
                        if(bond.get_type().split('-')[0]=='covalent'):
                            self.__ModelGraph.add_edge( main_atom.get_id(), temp_atom.get_id() , id = counter )
                    except:
                        None
        
        
        #CIF FILE (https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/struct_conn.html)
        if(counter==-1):
            data = '\n'.join(self.get_parent().get_data())

            for i in data.split('#\nloop_\n'):
                if(i[:12] == '_struct_conn' ):
                    col_seq = []
                    for j in i.split('#'):
                        if(j[:12] == '_struct_conn' ):
                            for k in j.split('\n'):
                                try:
                                    if(k[0]=='_'):
                                        col_seq.append( k )
                                    else:
                                        temp_dict = {col_seq[numl]:l for numl,l in enumerate(k.split())}
                                        #Actual approach: find chains -> residue/hetatm -> get the atom with the same name
                                        chain1, chain2 = self[temp_dict['_struct_conn.ptnr1_label_asym_id']], self[temp_dict['_struct_conn.ptnr2_label_asym_id']]
                                                                                
                                        ptnr1 = chain1.get_residue( int(temp_dict['_struct_conn.ptnr1_auth_seq_id']) )
                                        if(ptnr1 is None):
                                            ptnr1 = chain1.get_hetmol( int(temp_dict['_struct_conn.ptnr1_auth_seq_id']) )
                                            if(ptnr1 is None):
                                                ptnr1 = chain1.get_hetmol( temp_dict['_struct_conn.ptnr1_label_comp_id']+temp_dict['_struct_conn.ptnr1_auth_seq_id'] )

                                        ptnr2 = chain2.get_residue( int(temp_dict['_struct_conn.ptnr2_auth_seq_id']) )
                                        if(ptnr2 is None):
                                            ptnr2 = chain2.get_hetmol( int(temp_dict['_struct_conn.ptnr2_auth_seq_id']) )
                                            if(ptnr2 is None):
                                                ptnr2 = chain2.get_hetmol( temp_dict['_struct_conn.ptnr2_label_comp_id']+temp_dict['_struct_conn.ptnr2_auth_seq_id'] )
                                        
                                        if(ptnr1 is not None and ptnr2 is not None):
                                            atm1 = ptnr1.get_atom(temp_dict['_struct_conn.ptnr1_label_atom_id'])
                                            atm2 = ptnr2.get_atom(temp_dict['_struct_conn.ptnr2_label_atom_id'])

                                            cif_vocab2pacman= {'covale':'covalent' , 'disulf':'covalent-single', 'hydrog':'hydrogen', 'metalc':'non-covalent'}

                                            bond_type = None
                                            try:        
                                                bond_type = cif_vocab2pacman[temp_dict['_struct_conn.conn_type_id']]
                                            except:
                                                bond_type = 'other'
                                            
                                            counter+=1
                                            bond = Bond(counter, atm1, atm2, bond_type, source='CONECT-section')
                                            self.__AllBonds[counter] = bond
                                            atm1.set_bond(bond)
                                            atm2.set_bond(bond)
                                            self.__ModelGraph.add_node( atm1.get_id() )
                                            self.__ModelGraph.add_node( atm2.get_id() )
                                            if(bond.get_type().split('-')[0]=='covalent'):
                                                self.__ModelGraph.add_edge( atm1.get_id(), atm2.get_id() , id = counter )
                                        else:
                                            try:
                                                logging.debug('Some atom(s) are incorrectly identified and not added in the mmCIF file while adding bond information from the _struct_conn field. Atoms'+str(atm1.get_id())+' '+str(atm2.get_id()) )
                                            except:
                                                logging.debug('Some atom(s) are incorrectly identified and not added in the mmCIF file while adding bond information from the _struct_conn field.')

                                except Exception as e:
                                    #Check whats up with string error (very minor)
                                    None

        #Generate bonds from the default settings (No connect records)
        for chain in self.get_chains():
            try:
                resi = [i for i in chain.get_residues()]
            except:
                logging.info('Residues are not found for bonds without CONNECT records (default bonds) for the chain: '+str(chain.get_id()))
                continue
            #This section is for the default bonds (Not in the connect record)
            for i in resi:
                try:
                    for j in aa_connectivity[i.get_name()]:
                        try:
                            atom1, atom2 =  i.get_atom(j[0]), i.get_atom(j[1])
                            if(atom1 is not None and atom2 is not None and atom1.get_parent().get_parent().get_id() == atom2.get_parent().get_parent().get_id() ):
                                counter += 1
                                if(j[2]=='SING'):
                                    bond = Bond(counter, atom1, atom2, 'covalent-single', source='RCSB/aa-variants-v1.cif')
                                elif(j[2]=='DOUB'):
                                    bond = Bond(counter, atom1, atom2, 'covalent-double', source='RCSB/aa-variants-v1.cif')
                                else:
                                    bond = Bond(counter, atom1, atom2, 'covalent', source='RCSB/aa-variants-v1.cif')
                                self.__AllBonds[counter] = bond
                                atom1.set_bond(bond)
                                atom2.set_bond(bond)

                                self.__ModelGraph.add_node( atom1.get_id() )
                                self.__ModelGraph.add_node( atom2.get_id() )
                                if(bond.get_type().split('-')[0]=='covalent'):
                                    self.__ModelGraph.add_edge( atom1.get_id(), atom2.get_id() , id = counter )
                        except Exception as e:
                            #Bond pair not found
                            None
                except:
                    logging.debug('Residue Number|Name|Chain '+str(i.get_id())+'|'+str(i.get_name())+'|'+str(i.get_parent().get_id())+' is not a standard amino acid; sidechain bonds are not calculated.')

            #Peptide bonds (Assumption: All the reisudes are in the incremental/decremental order)
            for i in range(0,len(resi)-1):
                if( abs(resi[i].get_id() - resi[i+1].get_id()) == 1 ):
                    counter+=1
                    C, N = resi[i].get_atom('C'), resi[i+1].get_atom('N')
                    bond = Bond(counter, C, N, 'covalent-single', source='PACKMAN Peptide Bond calculation' )
                    self.__AllBonds[counter] = bond
                    C.set_bond(bond)
                    N.set_bond(bond)
                    self.__ModelGraph.add_edge( C.get_id(), N.get_id() , id = counter )
                else:
                    logging.info('The peptide bond following residue is missing: '+str(resi[i].get_id())+' Chain: '+resi[i].get_parent().get_id())
        
        try:
            None
            #Successful run
            '''
            set_to = -90
            print('Before:', self.calculate_torsion(common_bond, neighbor1=227, neighbor2=231, radians=False) )
            self.calculate_rotation(common_bond, set_to, neighbor1=227, neighbor2=231, radians=False)
            print('Set:',set_to)
            print('After:', self.calculate_torsion(common_bond, neighbor1=227, neighbor2=231, radians=False) )
            '''
        except Exception as e:
            print(str(e))


    #Check Function
    def check_clashes(self, distance: float=0.77) -> int:
        """Check if any atoms are too close to each other. This is important since too close atoms in the elastic network models can be very bad for the results.

        Notes:
            * This function will be moved to the molecule manipulation package later

        Args:
            distance (float): The distance cutoff user wishes to defined as a clash radius (default:0.77; max bond length)
        
        Returns:
            clashes (int)   : Number of clashes present according to the set cutoff.
        """
        from scipy.spatial import KDTree
        all_atoms=[i for i in self.get_atoms()]
        T=KDTree([i.get_location() for i in all_atoms])
        return len(T.query_pairs(distance))
