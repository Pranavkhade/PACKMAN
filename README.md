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
  -pdbid PDB_ID PDB_ID, --pdbid PDB_ID PDB_ID
                        (1) PDB ID of the input file (2) Location and Name by
                        which you wish to save the downloaded file
  -filename FILENAME, --filename FILENAME
                        Path and filename of the PDB file
  --e_clusters NumberOfEccentricityClusters
                        Recommended: 4, Please refer to the paper for more
                        details
  --minhnglen MinimumHingeLength
                        Recommended: 5, Please refer to the paper for more
                        details
  --chain CHAIN         Enter The Chain ID
  --generateobj {yes,no}
                        Select yes if you wish to generate the .obj file
```
### Examples

1. `python PACKMAN.py 2.8 --pdbid 1prw 1prw.pdb`
First parameter 2.8 is the alpha value parameter which can be changed. for --pdbid, the first parameter is the PDB ID the user submits to download a corresponding PDB file and second is the name of the file you wish to save the downloaded PDB file.


2. `python PACKMAN.py 2.8 --filename 1prw.pdb`
First parameter 2.8 is the alpha value parameter which can be changed. for --filename, the parameter is location and name of the PDB file.
