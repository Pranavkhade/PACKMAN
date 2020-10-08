# -*- coding: utf-8 -*-
# Author: Pranav Khade, Iowa State University
# Please read the project licence file for the Copyrights.

"""The 'packman' is a collection of subpackages built on the packman.molecule API.

Please check the corresponding packages and tutorials for more information about the package use.
"""
from .molecule import *
from .apps import *
from .bin import *
from .anm import *
from .constants import *
from .utilities import *


#VERSION CHANGE HERE CHANGES IT IN docs AND setup.py; STILL NEED TO CHANGE IT IN requirement.txt
__version__='1.1.4'