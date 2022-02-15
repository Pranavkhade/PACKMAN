# -*- coding: utf-8 -*-
"""The 'Chain' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Chain' object documentation for details.

Citation:
    Pranav M Khade, Robert L Jernigan, PACKMAN-Molecule: Python Toolbox for Structural Bioinformatics, Bioinformatics Advances, 2022;, vbac007, https://doi.org/10.1093/bioadv/vbac007

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

import numpy
import logging



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
        self.__Residues = {}
        self.__HetMols = {}
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
                self.__Residues[Number]=Entity
        elif(Type=='HetMol'):
            try:
                self.__HetMols[Number]=Entity
            except:
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
        """Get the generator of corresponding 'atom' objects of the residues of the 'Chain'

        Returns:
            generator of 'atom' objects if successful, None otherwise.
        """
        for i in sorted(self.__Residues.keys()):
            for j in self.__Residues[i].get_atoms():
                yield j
    
    def get_hetatoms(self):
        """Get the generator of corresponding 'atom' objects of the hetmols of the 'Chain'

        Returns:
            generator of 'atom' objects if successful, None otherwise.
        """
        for i in sorted(self.__HetMols.keys()):
            for j in self.__HetMols[i].get_atoms():
                yield j
    
    def get_atom(self, idx):
        """Get the atom of the given ID.

        Note: - This function is different from :py:func:`packman.molecule.chain.get_atoms` and also :py:func:`packman.molecule.residue.get_atom`
              - If the PDB file is constructed manually/ has multiple atoms of the same ID, the first instance of the atom with that id is returned. Please avoid saving two atoms with same ID in a same structure file in a given frame/model.

        Args:
            idx (int): Get atom by the id
        
        Returns:
            atom (:py:class:`packman.molecule.Atom`): Atom of the given ID if successful; None otherwise.
        """
        the_atom = None
        for i in self.__Residues:
            the_atom =  self.__Residues[i].get_atom(idx)
            if(the_atom!=None):
                break
            
        if(the_atom==None):
            for i in self.__HetMols:
                the_atom =  self.__HetMols[i].get_atom(idx)
                if(the_atom!=None):
                    break
        if(the_atom!=None):
            return the_atom
        else:
            logging.info('The atom with the given ID is not found in the Residues/HetMols')
            return None
    
    def get_residue(self,idx):
        """Get the residue of the given ID.

        Note: - This function is different from :py:func:`packman.molecule.chain.get_residues`
              - If the PDB file is constructed manually/ has multiple residues of the same ID, the first instance of the residue with that id is returned. Please avoid saving two residues with same ID in a same structure file in a given frame/model.

        Args:
            idx (int): Get residue by the id
        
        Returns:
            residue (:py:class:`packman.molecule.Residue`): Residue of the given ID if successful; None otherwise.
        """
        for i in self.get_residues():
            if(i.get_id() == idx):
                return i
                break

    def get_hetmol(self,idx):
        """Get the hetmol of the given ID.

        Note: - This function is different from :py:func:`packman.molecule.chain.get_hetmols`
              - If the PDB file is constructed manually/ has multiple hetmol of the same ID, the first instance of the hetmol with that id is returned. Please avoid saving two hetmol with same ID in a same structure file in a given frame/model.

        Args:
            idx (int): Get hetmol by the id
        
        Returns:
            residue (:py:class:`packman.molecule.HetMol`): HetMol of the given ID if successful; None otherwise.
        """
        for i in self.get_hetmols():
            if(i.get_id() == idx):
                return i
                break
    
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