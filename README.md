# PACKMAN: PACKing and Motion ANalysis.
Description

### Prerequisites:
#### Required modules. 

Modules are available in most GNU/Linux distributions, or from their respective websites.

* [numpy](http://www.numpy.org/)

* [scipy](https://www.scipy.org/)

* [networkx](https://networkx.github.io/)

* [mlxtend](http://rasbt.github.io/mlxtend/)

* [sklearn](https://scikit-learn.org/stable/)

#### Hinge Plane Visualization

If you wish to visualize the plane of the hinge, please refer to the following article: 

* [Plane Wizard (PyMOL wiki)](https://pymolwiki.org/index.php/Plane_Wizard)


### Installation
Description
<!--
Installing from source
```
git clone https://github.com/Pranavkhade/PACKMAN
cd PACKMAN
python setup.py install
```

Installing with pip
```
pip install PACKMAN
```
OR
```
pip install git+git://github.com/Pranavkhade/PACKMAN
```
-->

### Files and instructions
Description

```
usage: PACKMAN.py [-h] (-pdbid PDB_ID PDB_ID | -filename FILENAME)
                  [--e_clusters NumberOfEccentricityClusters]
                  [--minhnglen MinimumHingeLength] [--chain CHAIN]
                  [--generateobj {yes,no}]
                  AlphaValue

PACKMAN: PACKing and Motion ANalysis. (https://github.com/Pranavkhade/PACKMAN)

positional arguments:
  AlphaValue            Recommended: Start from 2 and keep increasing the
                        parameter value till the hinges become redundant
                        compared to the previous alpha values (Typically
                        around 5), Please refer to the paper for more details

optional arguments:
  -h, --help            show this help message and exit
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
### Examples

1. `python PPACKMAN.py --pdbid 1prw 2.8 1prw.pdb`
for --pdbid, the parameter is the PDB ID the user submits to download a corresponding PDB file. First positional parameter 2.8 is the alpha value parameter which can be changed. and second is the name of the file you wish to save the downloaded PDB file.


2. `python PPACKMAN.py 2.8 1prw.pdb`
First parameter 2.8 is the alpha value parameter which can be changed. second parameter is the parameter is location and name of the PDB file.