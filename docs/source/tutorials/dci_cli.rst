.. _tutorials_dci_cli:


DCI Analysis CLI
================

This tutorial familiarises users with the Command-line Interface (CLI) of the DCI Analysis. The same functionality is also available in the PACKMAN API ( :ref:`tutorials_dci_api` )

PACKMAN DCI Analysis Web Server:

How to cite::

    Kumar, A., Khade, P. M., Dorman, K. S., & Jernigan, R. L. (2022). Coarse-graining protein structures into their 
    dynamic communities with DCI, a dynamic community identifier. 
    Bioinformatics, 38(10). https://doi.org/10.1093/bioinformatics/btac159


CLI INSTRUCTIONS
----------------


Following is the PACKMAN hinge prediction interface description::

    usage: packman dci [-h] [-pdbid PDB_ID] [-chain CHAIN] [-cutoff CUTOFF]
                        [-n_com N_COM]
                        FILENAME

    positional arguments:
    FILENAME              Path and filename of the PDB file.

    optional arguments:
    -h, --help            show this help message and exit
    -pdbid PDB_ID, --pdbid PDB_ID
                            If provided, the PBD with this ID will be downloaded
                            and saved to FILENAME.
    -chain CHAIN, --chain CHAIN
                            Enter The Chain ID
    -cutoff CUTOFF, --cutoff CUTOFF
                            Enter the cutoff for DCI. (Read the Publication for
                            more details)
    -n_com N_COM, --n_com N_COM
                            Enter the number of communities. (Read the Publication
                            for more details)

EXAMPLES
--------
1. `python -m packman dci 1prw.pdb --chain A` The options are self-explanatory. This is the example where PDB file is already present.
2. `python -m packman dci 1prw.pdb` Chains are not defined, and therefore DCI is calculated for all the chains.