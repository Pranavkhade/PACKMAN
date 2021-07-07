# -*- coding: utf-8 -*-
"""The Entropy method(s) host file.

This is file information, not the class information. This information is only for the API developers.
Please read the corrosponding object documentation for details.

Example::


Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""
from .. import molecule
from ..constants import vdw_surface_bondi

import numpy

import logging

from scipy.spatial import Voronoi, ConvexHull, KDTree
from scipy.constants import R

from multiprocessing.dummy import Pool as ThreadPool

class PackingEntropy():
    """This class contains all the methods required to obtain a protein complex's entropy.
    Given a group of atoms in the :mod:'packman.molecule.Atom' objects, the entropy for the each amino acid will be returned.
    The 'chains' argument should be used when the user wants to restrict the analysis to a chain or group of chains rather than the whole structure.
    Args:
        atoms ([packman.molecule.Atom]) : The group of atoms user wisher to calculate Packing Entropy with.
        chains ([str]/str)              : Chain IDs for the Entropy calculation (None means all the chains are included; single string means only one chain ID; multiple chains should be an array of strings).
        probe_size (float)              : Radius of the probe to generate the surface points (This value should not be less than 1;Read the Publication for more details)
        onspherepoints (int)            : Number of points to be generated around each point for the surface (Read the Publication for more details)
    """
    def __init__(self, atoms, chains=None, probe_size=1.4, onspherepoints=30):
        if(chains==None):
            self.atoms = [i for i in atoms]
        else:
            #Single or multiple chain IDs provided. If there are multiple chain ids, they should be in an array form.
            try:
                self.atoms = [i for i in atoms if i.get_parent().get_parent().get_id() in chains]
            except:
                self.atoms = [i for i in atoms if i.get_parent().get_parent().get_id() == chains]

        self.coordinates=[i.get_location() for i in self.atoms]
        self.probe_size = probe_size
        self.kd_tree = KDTree(self.coordinates)
        self.onspherepoints = onspherepoints

        self.residues = list(set([i.get_parent() for i in self.atoms]))

        #Need to be calculated
        self.surface_points = None

        #Generate Surface Points
        self.calculate_surafacepoints()

        #Initiate
        self.calculate_entropy()

    #Get functions
    def get_total_packing_fraction(self):
        """The sum of Packing Fraction for the Residues in the provided atoms.

        Returns:
            The sum of the Residue Packing Fraction (float)
        """
        return numpy.sum( [i.get_property('PackingFraction') for i in self.residues] )

    def get_total_entropy(self):
        """The sum of Packing Entropies for the Residues in the provided atoms.
        Returns:
            The sum of Residue Entropies (float)
        """
        return numpy.sum( [i.get_entropy('PackingEntropy') for i in self.residues] )
    
    def get_total_chain_entropy(self,chain):
        """The sum of Packing Entropies for the Residues in the provided atoms and chain.
        
        Please note that the entropy for the chain you are selecting might not exist if you have not calculated it properly. Please make sure to run calculate_entropy() function properly.
        Args:
            chain (str) : chain ID of for the chain you wish to get the entropy.
        """
        return numpy.sum( [i.get_entropy('PackingEntropy') for i in self.residues if i.get_parent().get_id()==chain] )
        
    def get_surafacepoints(self):
        """Get the surface points around the given set of atoms in the protein.
        Returns:
            [[float]]: Array of 3D points around the given set of atoms in the protein.
        """
        return self.surface_points

    #Calculate Functions
    def calculate_spherepoints(self,atom):
        """Given a single point, this function generates point cloud around the given points.
        Args:
            point ([float]) : 3D coordinate around which the sphere of point cloud to be generated.
        """
        point = atom.get_location()
        indices = numpy.arange(0, self.onspherepoints, dtype=float) + 0.5
        phi = numpy.arccos(1 - 2*indices/self.onspherepoints)
        theta = numpy.pi * (1 + 5**0.5) * indices
        x, y, z = numpy.cos(theta) * numpy.sin(phi), numpy.sin(theta) * numpy.sin(phi), numpy.cos(phi)
        sphere_multiplier = ( self.probe_size + vdw_surface_bondi[atom.get_element()] ) * 2
        x, y, z = sphere_multiplier*x, sphere_multiplier*y, sphere_multiplier*z
        x=x+point[0]
        y=y+point[1]
        z=z+point[2]
        data = zip(x.ravel(), y.ravel(), z.ravel())
        selected=[]
        for _ in data:
            if(self.kd_tree.query(_,k=2)[0][1] > sphere_multiplier):
                selected.append(_)
        return selected

    def calculate_surafacepoints(self):
        """Calculate the surface points with the current setup.
        """
        pool = ThreadPool()
        surface=pool.map( self.calculate_spherepoints, self.atoms )
        surface=[item for sublist in surface for item in sublist]
        self.surface_points = surface

    def calculate_entropy(self):
        """Calculate the Packing Entropy with the current setup.
        """
        logging.info("Packing Entropy calculation started.")

        #All the points (surface+protein)
        points=numpy.concatenate((self.surface_points,[i.get_location() for i in self.atoms]))
        voronoi=Voronoi(points)

        AvailableVolume = {j:[] for j in set([i.get_parent() for i in self.atoms])}
        for numi, i in enumerate(self.atoms):
            atom_cell=[]
            for j in voronoi.regions[  voronoi.point_region[ numi+len(self.surface_points) ]  ]:
                atom_cell.append(voronoi.vertices[j])
            AvailableVolume[i.get_parent()].extend(atom_cell)
        
        #Each value= Total volume of voronoi cells of the atoms of the residue / Total volume of the residues (Both convex hull volumes)
        PackingFraction={}
        for i in AvailableVolume:
            PackingFraction[i] = float(ConvexHull([j.get_location() for j in i.get_atoms()]).volume) / ConvexHull(AvailableVolume[i]).volume
            
        #R = Gas constant
        for i in PackingFraction:
            #Old Model:
            #ResiduePackingEntropy = R * -PackingFraction[i] * numpy.log2(PackingFraction[i])
            #New Model:
            ResiduePackingEntropy = - R * numpy.log10( PackingFraction[i] )
            
            i.set_entropy('PackingEntropy', ResiduePackingEntropy)
            i.set_property('PackingFraction',PackingFraction[i])
        
        logging.info("Packing Entropy calculated and assigned to the respective 'Residues' successfully.")