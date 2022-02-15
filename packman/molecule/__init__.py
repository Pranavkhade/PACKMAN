# -*- coding: utf-8 -*-
# Author: Pranav Khade, Iowa State University
# Please read the project licence file for the Copyrights.

"""The 'packman.molecule' module is used to read, write, manipulate and analyze the molecule.

This module is the base of the tool packman. It is used as a base module for all the packman utilities
such as HingePrediction, Compliance and Right Domain ANM. The molecule module can also be an API to utilize the objects such
as Atom, Residue, Chain, Model and Protein. Please read the documentation and tutorials for more details.

Citation:
    Pranav M Khade, Robert L Jernigan, PACKMAN-Molecule: Python Toolbox for Structural Bioinformatics, Bioinformatics Advances, 2022;, vbac007, https://doi.org/10.1093/bioadv/vbac007

Notes:
    * Tutorial link: https://py-packman.readthedocs.io/en/latest/tutorials/molecule.html#tutorials-molecule

Example:
    To Load the molecule::
        from packman import molecule
        molecule.download_structure('1prw','1prw.pdb')
        mol = molecule.load_structure('1prw.pdb')

Todo:
    * Add new features
    * Use ``sphinx.ext.todo`` extension
    * Add Tutorial Link
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
from .bond import Bond

from .hetmol import HetMol

from .annotations import Hinge