# -*- coding: utf-8 -*-
# Author: Pranav Khade, Iowa State University
# Please read the project licence file for the Copyrights.

"""The 'packman.molecule' module is used to read, write, manipulate and analyze the molecule.

This module is the base of the tool packman. It is used as a base module for all the packman utilities
such as HingePrediction, Compliance and Right Domain ANM. The molecule module can also be an API to utilize the objects such
as Atom, Residue, Chain, Model and Protein. Please read the documentation and tutorials for more details.

Example:
    To Load the 
    
    >>>from packman import molecule
    >>>molecule.download_structure('1prw','1prw.pdb')
    >>>molecule.load_structure('1prw.pdb')
    <packman.molecule.protein.Protein object at Memory>

Todo:
    * Add new features
    * Use ``sphinx.ext.todo`` extension

"""


#Parser Functions
from .molecule import download_structure
from .molecule import load_structure

#Building Functions
from .protein import Protein
from .model import Model

from .chain import Chain
from .residue import Residue
from .atom import Atom

from .hetmol import HetMol
from .hetatom import HetAtom

from .annotations import Hinge