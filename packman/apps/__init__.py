# -*- coding: utf-8 -*-
# Author: Pranav Khade, Iowa State University
# Please read the project licence file for the Copyrights.

"""The 'packman.apps' module is a collection of applications built on packman.molecule API.

Notes:
    * Current apps list: - predict_hinge : A program to predict the hinge on the molecule given the atoms and relevent parameters

Example::
    * Review the packman.bin.PACKMAN.py file for the app use.

Todo:
    * Add new features
    * Use ``sphinx.ext.todo`` extension
    * Clean the code
    * Redesign print statements
    * Add Tutorial Link

"""

from .predict_hinge import predict_hinge, hinge_cli
from .hdanm import hdanm_cli
from .calculate_entropy import entropy_cli
from .dci import DCI, dci_cli