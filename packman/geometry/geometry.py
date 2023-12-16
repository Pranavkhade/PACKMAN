# -*- coding: utf-8 -*-
"""The gemotry objects host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'geometry' object documentation for details.

Example::

    from packman import geometry
    help( geometry )

Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Annotate 'error', 'cross'
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""
import numpy

from typing import TYPE_CHECKING, List, Iterable, Union, Tuple

from scipy.spatial import Delaunay

from networkx import Graph

from .. import Atom

def Circumsphere(Tetrahydron: List[Atom]) -> (float, float):
    """Get the Circumsphere of the set of four points.
    
    Given any four three dimentional points, this function calculates the features of the circumsphere having the given four points on it's surface.
    
    Args:
        Tetrahydron ([packman.molecule.Atom] or [[X,Y,Z]]): Either packman.molecule.Atom objects or 3D corrdinates. 
    
    Returns:
        [Centre, Radius] (float): The 3D coordinates of the geometrical center of the given four points, Radius of the circumsphere made up of given four points in that order.
    """
    alpha_mat, gamma_mat, Dx_mat, Dy_mat, Dz_mat=[], [], [], [], []
    for i in Tetrahydron:
        temp_coords = i.get_location()
        alpha_mat.append( [temp_coords[0],temp_coords[1],temp_coords[2],1] )
        gamma_mat.append( [temp_coords[0]**2+temp_coords[1]**2+temp_coords[2]**2,temp_coords[0],temp_coords[1],temp_coords[2]] )
        Dx_mat.append( [temp_coords[0]**2+temp_coords[1]**2+temp_coords[2]**2,temp_coords[1],temp_coords[2],1] )
        Dy_mat.append( [temp_coords[0]**2+temp_coords[1]**2+temp_coords[2]**2,temp_coords[0],temp_coords[2],1] )
        Dz_mat.append( [temp_coords[0]**2+temp_coords[1]**2+temp_coords[2]**2,temp_coords[0],temp_coords[1],1] )
    alpha = numpy.linalg.det(alpha_mat)
    gamma = numpy.linalg.det(gamma_mat)
    Dx = (+numpy.linalg.det(Dx_mat))
    Dy = (-numpy.linalg.det(Dy_mat))
    Dz = (+numpy.linalg.det(Dz_mat))
    Centre = numpy.array([Dx/2*alpha,Dy/2*alpha,Dz/2*alpha])
    Radius = numpy.sqrt(Dx**2 + Dy**2 + Dz**2 - 4*alpha*gamma)/(2*numpy.absolute(alpha))
    return Centre, Radius


def AlphaShape( atoms: Iterable[Atom], alpha: float, get_graph: bool=True ) -> Union[ List[Atom], Tuple[List[Atom], Graph] ]:
    """Get the Alpha Shape of the atoms.

    Given set of atoms as 'Atom' (packman.molecule.Atom) object(s) and the alpha value, this class will calculate and return the alpha shape tessellations. It will also write the .obj file if the filename is provided as an input.

    Note:
        * Tip: If you do not want to use the function multiple times to save computation, calculate it once with alpha = float('Inf') and then use the tessellations to calculate radius and save it as a dictionary to retrieve. Tessellations with any cutoff.
        * For more information on the alpha shape, read the following paper:
            EdelsbrunnerandE. P. M ̈ucke.Three-dimensional alpha shapes.
            Manuscript UIUCDCS-R-92-1734, Dept.Comput.Sci. ,Univ.Illinois, Urbana-Champaign, IL, 1992.
    
    Args:
        atoms (packman.molecule.Atom) : Iterable that returns Atom type.
        alpha (float)                 : Alpha cutoff value.
        get_graph (networkx.Graph)    : Return graph of the Atoms for the given alpha parameter.
    
    Returns:
        - Alpha Shape Tessellations                ; if 'get_graph' = False
        - Alpha Shape Tessellations, Protein Graph ; if 'get_graph' = True
    """
    def calculate_alphafitness(alpha: float, circumradius: float) -> bool:
        """Alpha Test as per the paper.
        
        Note:
            * Resides inside 'get_alphashape' function.

        Args:
            alpha (float)                : Alpha.                              (Read parent method description)
            circumradius (float)         : Circumradius of the circumsphere.   (Read parent method description)
        
        Returns:
            bool: True if alpha test is passed. False otherwise.
        """
        if( circumradius < alpha ):
            return True
        else:
            return False

    #get_alphashape function commands.
    DelaunayTesssellations = Delaunay( [i.get_location() for i in atoms] )
    AlphaShape   = []
    ProteinGraph = Graph()

    for a, b, c, d in DelaunayTesssellations.simplices:
        Tetrahydron = [atoms[a],atoms[b],atoms[c],atoms[d]]
        Centre, Radius = Circumsphere(Tetrahydron)
        #Alpha Test
        if( calculate_alphafitness(alpha,Radius) ):
            ProteinGraph.add_nodes_from( [a,b,c,d] )
            ProteinGraph.add_edges_from( [(a,b),(a,c),(a,d),(b,c),(b,d),(c,d)] )
            AlphaShape.append( [atoms[a],atoms[b],atoms[c],atoms[d]] )

    if(get_graph):
        return AlphaShape, ProteinGraph
    else:
        return AlphaShape
