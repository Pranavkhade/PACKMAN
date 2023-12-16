# -*- coding: utf-8 -*-
# Author: Pranav Khade
"""The 'packman.apps' module is a collection of applications built on packman.molecule API.

Note:
    Current apps list: 
    * predict_hinge : Predict Protein Hinges (https://doi.org/10.1016/j.jmb.2019.11.018)
    * hdanm: Hinge-Domain ANM (https://doi.org/10.1016/j.bpj.2021.10.017)
    * calculate_entropy: Calculate Protein Packing Entropy (https://doi.org/10.1021/acsomega.2c00999)
    * dci : Predict Protein Dynamics Communities (https://doi.org/10.1093/bioinformatics/btac159)

Example::
    * Review the packman.bin.PACKMAN.py file for the app use.
"""

from .predict_hinge import predict_hinge, hinge_cli
from .hdanm import hdanm_cli
from .calculate_entropy import entropy_cli
from .dci import DCI, dci_cli