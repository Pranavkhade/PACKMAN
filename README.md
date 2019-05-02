# PPACKMAN: Protein PACKing and Motion ANalysis.
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
git clone https://github.com/Pranavkhade/PPACKMAN
cd PPACKMAN
python setup.py install
```

Installing with pip
```
pip install PPACKMAN
```
OR
```
pip install git+git://github.com/Pranavkhade/PPACKMAN
```
-->

### Files and instructions
Description

```
usage: PPACKMAN.py [-h] [-pdbid PDB_ID] [--chain CHAIN]
                   [--generateobj GENERATEOBJ] [--outputfile OUTPUTFILE]
                   [--logfile LOGFILE] [--callbackurl CALLBACKURL]
                   [--nodeid NODEID]
                   AlphaValue FILENAME

PPACKMAN: Protein PACKing and Motion ANalysis.
(https://github.com/Pranavkhade/PPACKMAN)

positional arguments:
  AlphaValue            Recommended: 2.8 for closed; 4.5 for open form, Please
                        refer to the paper for more details
  FILENAME              Path and filename of the PDB file.

optional arguments:
  -h, --help            show this help message and exit
  -pdbid PDB_ID, --pdbid PDB_ID
                        If provided, the PBD with this ID will be downloaded
                        and saved to FILENAME.
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

1. `python PPACKMAN.py 2.8 --pdbid 1prw 1prw.pdb`
First parameter 2.8 is the alpha value parameter which can be changed. for --pdbid, the first parameter is the PDB ID the user submits to download a corresponding PDB file and second is the name of the file you wish to save the downloaded PDB file.


2. `python PPACKMAN.py 2.8 --filename 1prw.pdb`
First parameter 2.8 is the alpha value parameter which can be changed. for --filename, the parameter is location and name of the PDB file.