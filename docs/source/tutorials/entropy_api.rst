.. _tutorials_entropy_api:


Entropy Analysis API
====================
This tutorial familiarises users with the Application Programming Interface (API) of the Entropy Analysis. The same functionality is also available in the PACKMAN CLI ( :ref:`tutorials_entropy_cli` ) and GUI (:ref:`tutorials_gui`).

PACKMAN Entropy Analysis Web Server: https://packing-entropy.bb.iastate.edu/

How to cite::

    Pranav M. Khade and Robert L. Jernigan. Entropies Derived from the Packing Geometries within a Single Protein Structure. 
    ACS Omega 2022 7 (24), 20719-20730 DOI: 10.1021/acsomega.2c00999

There are two ways the entropy can be calculated. The second way is when the user explicitly wants to leave specific atoms/residues out of the entropy calculation. Chain exclusion can be done via the first method as well. Also, the user has access to the extra get, set, and calculate methods of the Entropy objects using the second way.

The First Way (Through packman.molecule Objects)
------------------------------------------------

Step 1
~~~~~~
1. Import the 'molecule' submodule from PACKMAN. (If PACKMAN is not installed, please follow the link: https://github.com/Pranavkhade/PACKMAN)

2. Download the structure from PDB and save it with the appropriate extension.

3. Load the structure using the submodule::
    
    #Step 1.1
    from packman import molecule

    #Step 1.2
    molecule.download_structure('1LF7')
    #OR (Default is CIF format; to change it to PDB)
    #molecule.download_structure('1LF7',ftype = 'pdb')

    #Step 1.3
    mol=molecule.load_structure('1LF7.cif')
    #OR
    #mol=molecule.load_structure('1LF7.pdb')


Step 2
~~~~~~
Please note that selecting --chain option is used to calculate the entropy with or without the other chains. For example, if a protein has three chains A, B & C, and the user wants to calculate the entropy of chains A & B in the absence of chain C, the user can use --chain A, B parameter-option do so. However, the presence of the chain C is not an issue, but the user wants to calculate the entropy of chain A & B (with C present), --chain option can be ignored, and chain column in the output should be used to select only A & B chains. The use of API is recommended to control these types of situations efficiently that require more control. Also, we recommend reading the publication for more details on the other parameters.

Do any of the following to calculate the Entropy::

    #Via a Model (frame) object
    mol[0].calculate_entropy('PackingEntropy',chains=None, probe_size=1.4, onspherepoints=30)

    #OR

    #Via a Chain object
    mol[0]['A'].calculate_entropy('PackingEntropy',chains=None, probe_size=1.4, onspherepoints=30)

    #OR

    #Via a Residue object (Other parameters are not defined to show they are optional)
    mol[0].get_residues()[0].calculate_entropy('PackingEntropy')

Step 3
~~~~~~
Please note that the get_entropy (This step) is just a retrieval method. The way entropy is calculated defined in Step 2.

Retrieve the calculated entropy in the exactly same way except it will be Residue/Chain/Model (Frame) specific depending on the object it is being retrieved from.::

    #For the Entropy of the 0th frame
    mol[0].get_entropy('PackingEntropy')

    #For the Entropy of the Chain A
    mol[0]['A'].get_entropy('PackingEntropy')

    #For the Entropy of the First residue in the sequence
    mol[0].get_residues()[0].get_entropy('PackingEntropy')


The Second Way (Via 'Entropy' Objects)
--------------------------------------

For the example, we are going to use 'PackingEntropy' object. However, other entropies can be calculated in a similar way. The user can use this way to leave out specific atoms/residues and even chains (chains are also possible to leave out in the first way) in the entropy calculation. Also, the user has access to the extra get, set, and calculate methods of the Entropy objects using this way.

Step 1
~~~~~~
1. Import the 'molecule' submodule from PACKMAN. (If PACKMAN is not installed, please follow the link: https://github.com/Pranavkhade/PACKMAN)

2. Download the structure from PDB and save it with the appropriate extension.

3. Load the structure using the submodule::

    #Step 1.1
    from packman import molecule

    #Step 1.2
    molecule.download_structure('1LF7')
    #OR (Default is CIF format; to change it to PDB)
    #molecule.download_structure('1LF7',ftype = 'pdb')

    #Step 1.3
    mol=molecule.load_structure('1LF7.cif')
    #OR
    #mol=molecule.load_structure('1LF7.pdb')

Step 2
~~~~~~

Please note that selecting --chain option is used to calculate the entropy with or without the other chains. For example, if a protein has three chains A, B & C, and the user wants to calculate the entropy of chains A & B in the absence of chain C, the user can use --chain A, B parameter-option do so. However, the presence of the chain C is not an issue, but the user wants to calculate the entropy of chain A & B (with C present), --chain option can be ignored, and chain column in the output should be used to select only A & B chains. The use of API is recommended to control these types of situations efficiently that require more control. Also, we recommend reading the publication for more details on the other parameters.

1. Import the 'PackingEntropy' (for example)

2. Use the 'PackingEntropy' object with specific 'Atoms' (that user can select or filter based on choice).::

    #Step 1
    from packman.entropy import PackingEntropy

    #Step 2
    result = PackingEntropy(mol[0].get_atoms(),chains='A,B',probe_size=1.4,onspherepoints=30)

Step 3
~~~~~~
Please note that the get_entropy (This step) is just a retrieval method. The way entropy is calculated defined in Step 2.

The entropy can be retrieved using the same procedure explained in Step 3 of the first way. However, the PackingEntropy also has get, set, and calculate methods that can be used. Please check the :mod:`packman.entropy.PackingEntropy` for more details.
