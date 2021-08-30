.. _tutorials_entropy_cli:


Entropy Analysis CLI
====================

This tutorial familiarises users with the Command-line Interface (CLI) of the Entropy Analysis. The same functionality is also available in the PACKMAN API ( :ref:`tutorials_entropy_api` ) and GUI (:ref:`tutorials_gui`).

PACKMAN Entropy Analysis Web Server: * Packing Entropy: https://packing-entropy.bb.iastate.edu/

How to cite::

    Paper Under Review.


CLI INSTRUCTIONS
----------------

Please note that selecting --chain option is used to calculate the entropy with or without the other chains. For example, if a protein has three chains A, B & C, and the user wants to calculate the entropy of chains A & B in the absence of chain C, the user can use --chain A, B parameter-option do so. However, the presence of the chain C is not an issue, but the user wants to calculate the entropy of chain A & B (with C present), --chain option can be ignored, and chain column in the output should be used to select only A & B chains. The use of API is recommended to control these types of situations efficiently that require more control. Also, we recommend reading the publication for more details on the other parameters.

Following is the PACKMAN hinge prediction interface description::

    usage: packman entropy [-h] [-type entropy_type] [-pdbid PDB_ID]
                            [--chains Chains to be used for the entropy calculation]
                            [--probe_size Size surface probe radius]
                            [--onspherepoints Number of points on a sphere]
                            FILENAME

    positional arguments:
    FILENAME              Path and filename of the PDB file.

    optional arguments:
    -h, --help            show this help message and exit
    -type entropy_type, --type entropy_type
                            Provide the Entropy type (Options: 1. PackingEntropy)
    -pdbid PDB_ID, --pdbid PDB_ID
                            If provided, the PBD with this ID will be downloaded
                            and saved to FILENAME.
    --chains Chains to be used for the entropy calculation
                            Recommended: None. Chain IDs for the Entropy
                            calculation (None means all the chains are included;
                            single string means only one chain ID; multiple chains
                            should be comma separated).
    --probe_size Size surface probe radius
                            Recommended: 1.4 (radius of a water molecule), Please
                            refer to the paper for more details
    --onspherepoints Number of points on a sphere
                            Recommended: 30. Number of points to be generated
                            around each point for the surface (Read the
                            Publication for more details)

EXAMPLES
--------

1. `python.exe -m packman entropy 1prw.pdb -type PackingEntropy --chains A,B` The options are self-explanatory. This is the example where PDB file is already present.
2. `python.exe -m packman entropy 1prw.pdb -type PackingEntropy` Chains are not defined, and therefore entropy is calculated for all the chains.