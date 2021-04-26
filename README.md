[![Build Status](https://travis-ci.com/Pranavkhade/PACKMAN.svg?branch=master)](https://travis-ci.com/Pranavkhade/PACKMAN) [![Documentation Status](https://readthedocs.org/projects/py-packman/badge/?version=latest)](https://py-packman.readthedocs.io/en/latest/?badge=latest) [![Downloads](https://pepy.tech/badge/py-packman)](https://pepy.tech/project/py-packman) [![PyPI version](https://badge.fury.io/py/py-packman.svg)](https://badge.fury.io/py/py-packman)


PACKMAN: PACKing and Motion ANalysis
------------------------------------
<img src="https://github.com/Pranavkhade/PACKMAN/blob/master/docs/_static/gallary/logo.gif" width="250">

This package focuses on studying molecular structures and their dynamics using a simple yet informative property known as Protein Packing. Over the last few years, we have worked on several techniques to capture and quantify the protein packing, resulting in a few publications. This package has all the code to repeat and further develop these techniques.


DOCUMENTATION
-------------
* [IMPORTANT : Documentation, Tutorials & More](https://py-packman.readthedocs.io)
* [Hinge Prediction Web Server](https://packman.bb.iastate.edu/)
* [hd-ANM Web Server](https://hdanm.bb.iastate.edu/)
* [Reference for Hinge Prediction](https://doi.org/10.1016/j.jmb.2019.11.018)
* [Reference for Compliance](https://doi.org/10.1002/prot.25968)
* [Reference for hd-ANM](coming_soon)

INSTALLATION
------------

1. Installing with pip (Recommended)
```
pip install py-packman
```

2. Installing from source
```
git clone https://github.com/Pranavkhade/PACKMAN
cd PACKMAN
python setup.py install
```

HOW TO USE
----------

PACKMAN, along with its components, can be accessed via Graphical User Interface (GUI), Command-line Interface (CLI), and Application programming interface (API).

For the GUI, please run the following command.
```
python -m packman gui
```
OR
```
python3 -m packman gui
```

For the CLI and API, please read the tutorials & documentation.

PREREQUISITES
-------------

* [numpy](http://www.numpy.org/)

* [scipy](https://www.scipy.org/)

* [networkx](https://networkx.github.io/)

* [mlxtend](http://rasbt.github.io/mlxtend/)

* [sklearn](https://scikit-learn.org/stable/)