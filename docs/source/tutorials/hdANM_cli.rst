.. _tutorials_hdanm_cli:


hdANM CLI
=========

This tutorial familiarises users with the Command-line Interface (CLI) of the hdanm Analysis. The same functionality is also available in the PACKMAN API ( :ref:`tutorials_hdANM:` ) and GUI (:ref:`tutorials_gui`).

PACKMAN hdANM Analysis Web Server: https://hdanm.bb.iastate.edu/

How to cite::

    Paper Under Review.


CLI INSTRUCTIONS
----------------

Please note that selecting --chain option is used to calculate the hdanm with or without the other chains. For example, if a protein has three chains A, B & C, and the user wants to calculate the hdanm of chains A & B in the absence of chain C, the user can use --chain A, B parameter-option do so. However, the presence of the chain C is not an issue, but the user wants to calculate the hdanm of chain A & B (with C present), --chain option can be ignored, and chain column in the output should be used to select only A & B chains. The use of API is recommended to control these types of situations efficiently that require more control. Also, we recommend reading the publication for more details on the other parameters.

Following is the PACKMAN hinge prediction interface description::

    usage: packman hdanm [-h] [-pdbid PDB_ID] [--chain CHAIN] [--dr DR]
                            [--power POWER] [--mass MASS] [--scale SCALE]
                            [--frames FRAMES] [--modes MODES] [--make_tar]
                            FILENAME HNG

    positional arguments:
    FILENAME              Path and filename of the PDB file.
    HNG                   Path and filename of the corresponding HNG file.

    optional arguments:
    -h, --help            show this help message and exit
    -pdbid PDB_ID, --pdbid PDB_ID
                            If provided, the PBD with this ID will be downloaded
                            and saved to FILENAME.
    --chain CHAIN         Enter The Chain ID
    --dr DR               Distance cutoff for the ANM.
    --power POWER         Power of the distance in non-parametric ANM.
    --mass MASS           Mass of the residue; unit or molecular weight
    --scale SCALE         movie scale
    --frames FRAMES       number of frames
    --modes MODES         how many modes
    --make_tar            package output files into a tar.gz file