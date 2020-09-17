PACKMAN: PACKing and Motion ANalysis
------------------------------------
<img src="https://github.com/Pranavkhade/PACKMAN/blob/master/docs/_static/gallary/logo.gif" width="250">

PACKMAN is a multiutility tool to study protein packing and its effect on protein dynamics. It currently has the functionality to identify the protein hinges (separating the domains). It can also be used to read, write, manipulate and analyze protein molecules and it's properties through its API.

DOCUMENTATION
-------------
* [IMPORTANT : Documentation, Tutorials & More](https://py-packman.readthedocs.io)
* [Hinge Prediction Web Server](https://packman.bb.iastate.edu/)
* [hd-ANM Web Server](coming_soon)
* [Reference for Hinge Prediction](https://doi.org/10.1016/j.jmb.2019.11.018)
* [Reference for Compliance](https://doi.org/10.1002/prot.25968)
* [Reference for hd-ANM](coming_soon)

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