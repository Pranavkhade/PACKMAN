PACKMAN: PACKing and Motion ANalysis
------------------------------------
<img src="https://github.com/Pranavkhade/PACKMAN/blob/master/logo.gif" width="250">

PACKMAN is a multiutility tool to study protein packing and its effect on protein dynamics. It currently has the functionality to identify the protein hinges (separating the domains). It can also be used to read, write, manipulate and analyze protein molecules and it's properties through its API.

DOCUMENTATION
-------------
* [Web Server](https://packman.bb.iastate.edu/)
* [Reference Hinge Prediction](https://doi.org/10.1016/j.jmb.2019.11.018)

INSTALLATION
------------

1. Installing from source
```
git clone https://github.com/Pranavkhade/PACKMAN
cd PACKMAN
python setup.py install
```

2. Installing with pip
```
pip install py-packman
```
OR
```
pip install git+git://github.com/Pranavkhade/PACKMAN
```

PREREQUISITES
-------------

* [numpy](http://www.numpy.org/)

* [scipy](https://www.scipy.org/)

* [networkx](https://networkx.github.io/)

* [mlxtend](http://rasbt.github.io/mlxtend/)

* [sklearn](https://scikit-learn.org/stable/)

OPTIONAL
---------

If you wish to visualize the plane of the hinge, please refer to the following article: 

* [Plane Wizard (PyMOL wiki)](https://pymolwiki.org/index.php/Plane_Wizard)

INSTRUCTIONS
------------

```
usage: packman [-h] [-pdbid PDB_ID]
                  [--e_clusters NumberOfEccentricityClusters]
                  [--minhnglen MinimumHingeLength] [--chain CHAIN]
                  [--generateobj GENERATEOBJ] [--outputfile OUTPUTFILE]
                  [--logfile LOGFILE] [--callbackurl CALLBACKURL]
                  [--nodeid NODEID]
                  AlphaValue FILENAME

PACKMAN: PACKing and Motion ANalysis. (https://github.com/Pranavkhade/PACKMAN)

positional arguments:
  AlphaValue            Recommended: 2.8 for closed; 4.5 for open form, Please
                        refer to the paper for more details
  FILENAME              Path and filename of the PDB file.

optional arguments:
  -h, --help            show this help message and exit
  -pdbid PDB_ID, --pdbid PDB_ID
                        If provided, the PBD with this ID will be downloaded
                        and saved to FILENAME.
  --e_clusters NumberOfEccentricityClusters
                        Recommended: 4, Please refer to the paper for more
                        details
  --minhnglen MinimumHingeLength
                        Recommended: 5, Please refer to the paper for more
                        details
  --chain CHAIN         Enter The Chain ID
  --generateobj GENERATEOBJ
                        Path and filename to save the .obj file at. Ignored
                        unless --chain is provided.

Web server parameters:
  Used by the web form

  --outputfile OUTPUTFILE
                        Path and filename write output to
  --logfile LOGFILE     Path and filename write log messages to
  --callbackurl CALLBACKURL
                        Optional callback url if this script was called from
                        Drupal.
  --nodeid NODEID       Optional node id if this script was called from
                        Drupal.

```

EXAMPLES
--------

1. `python -m packman --pdbid 1prw 2.8 1prw.pdb`
OR
`packman --pdbid 1prw 2.8 1prw.pdb`
for --pdbid, the parameter is the PDB ID the user submits to download a corresponding PDB file. First positional parameter 2.8 is the alpha value parameter which can be changed. and second is the name of the file you wish to save the downloaded PDB file.


2. `python -m packman 2.8 1prw.pdb`
OR
`python packman 2.8 1prw.pdb`
First parameter 2.8 is the alpha value parameter which can be changed. second parameter is the parameter is location and name of the PDB file.
