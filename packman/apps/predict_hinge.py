# -*- coding: utf-8 -*-
"""The 'predict_hinge' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'predict_hinge' object documentation for details.

Example::

    from packman.apps import predict_hinge
    help( predict_hinge )

Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Move the marked functions to the geometry module
    * Clear and annotate the subclasses properly
    * Annotate 'error', 'cross'
    * Finish optimizing the performance.
    * Convert GetExample() functions to get_example() to make it uniform with rest of the API.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""

import numpy
import functools
import logging

from networkx import eccentricity

from scipy.stats import mode
from scipy.optimize import minimize

from sklearn.cluster import KMeans
from mlxtend.evaluate import permutation_test
from itertools import groupby, count

from ..molecule import Hinge
from ..utilities import WriteOBJ
from ..geometry import AlphaShape


def predict_hinge(atoms, outputfile, Alpha=float('Inf'),method='alpha_shape',filename='Output.pdb',MinimumHingeLength=5,nclusters=4):
    """This function is used to carry out hinge prediction given the parameters.

    Notes:
        * The packman.bin.PACKMAN uses this function

        * Ideally, the alpha values should be scanned from 0 to 10 to obtain conclusively repetative hinges (Use '' for this purpose) 

        * Please refer to the following paper for the details on the algorithm and citation:
            Pranav M. Khade, Ambuj Kumar, Robert L. Jernigan, Characterizing and Predicting Protein Hinges for Mechanistic Insight,
            Journal of Molecular Biology, Volume 432, Issue 2, 2020, Pages 508-522, ISSN 0022-2836, https://doi.org/10.1016/j.jmb.2019.11.018.
            (http://www.sciencedirect.com/science/article/pii/S0022283619306837)
        
        * Tutorial Link: https://py-packman.readthedocs.io/en/latest/tutorials/hinge_predict.html#tutorials-predict-hinge

        * The predicted hinges are stored in packman.molecule.Chain object as an packman.molecule.Hinge object.

    Todo:
        * Add reference to the function in the description which scans multiple alpha values for conclusive hinges.
    
    Args:
        atoms ([packman.molecule.Atom])   : PACKMAN uses backbone atoms of the protein. However, any number and type of atoms can be used (Alpha value range will change)
        outputfile (file)                 : Output File.
        Alpha (float, optional)           : Please refer to the paper for this parameter. Defaults to float('Inf').
        method (str, optional)            : Please refer to the paper for this parameter. Defaults to 'alpha_shape'.
        filename (str, optional)          : Please refer to the paper for this parameter. Defaults to 'Output.pdb'.
        MinimumHingeLength (int, optional): Please refer to the paper for this parameter. Defaults to 5.
        nclusters (int, optional)         : Please refer to the paper for this parameter. Defaults to 4.
    
    Returns:
        alpha_shape              : The Alpha Shape (Subset of Delaunay Tesselations) 
    """
    if(method=='alpha_shape'):

        def get_statistics(atoms,SelectedHingeResidues,filename='Output'):
            """This sub-method is used to get the statistical data on the hinges and print it into a file.
            
            Notes:
                * * Function level: 1 (1 being top)
                * Do something about the output file

            Args:
                atoms ([packman.molecule.Atom])                   : Set of atoms. (Read parent method description)
                SelectedHingeResidues ([packman.molecule.Residue]): Predicted hinge residues. 
                filename (str, optional)                          : Output file name. Defaults to 'Output'.
            
            Returns:
                [p-value, stats] (float): p-value of the predicted hinge, statistics of the hinge (in that order)
            """
            hinge_atoms=[i.get_backbone() for i in SelectedHingeResidues]
            hinge_atoms=[item for sublist in hinge_atoms for item in sublist]
            non_hinge_atoms=list(set([i for i in atoms])-set(hinge_atoms))
            all_atoms_bfactor=[i.get_bfactor() for i in atoms]
            hinge_atoms_bfactor=[i.get_bfactor() for i in hinge_atoms]
            non_hinge_atoms_bfactor=[i.get_bfactor() for i in non_hinge_atoms]

            return_stats=[]

            outputfile.write('\nSTATISTICS\n\t\tN\tMin\tMax\tMean\tMode\tMedian\tSTDDev\n')
            return_stats.append(['','N','Min','Max','Mean','Mode','Median','STDDev'])
            outputfile.write('Total   '+'\t'+str(len(all_atoms_bfactor))+'\t'+str(numpy.min(all_atoms_bfactor))+'\t'+str(numpy.max(all_atoms_bfactor))+'\t'+str(numpy.mean(all_atoms_bfactor))+'\t'+str(mode(all_atoms_bfactor)[0][0])+'\t'+str(numpy.median(all_atoms_bfactor))+'\t'+str(numpy.std(all_atoms_bfactor))+'\n')
            return_stats.append(['Total',len(all_atoms_bfactor),numpy.min(all_atoms_bfactor),numpy.max(all_atoms_bfactor),numpy.mean(all_atoms_bfactor),mode(all_atoms_bfactor)[0][0],numpy.median(all_atoms_bfactor),numpy.std(all_atoms_bfactor)])
            outputfile.write('Hinge   '+'\t'+str(len(hinge_atoms_bfactor))+'\t'+str(numpy.min(hinge_atoms_bfactor))+'\t'+str(numpy.max(hinge_atoms_bfactor))+'\t'+str(numpy.mean(hinge_atoms_bfactor))+'\t'+str(mode(hinge_atoms_bfactor)[0][0])+'\t'+str(numpy.median(hinge_atoms_bfactor))+'\t'+str(numpy.std(hinge_atoms_bfactor))+'\n')
            return_stats.append(['Hinge',len(hinge_atoms_bfactor),numpy.min(hinge_atoms_bfactor),numpy.max(hinge_atoms_bfactor),numpy.mean(hinge_atoms_bfactor),mode(hinge_atoms_bfactor)[0][0],numpy.median(hinge_atoms_bfactor),numpy.std(hinge_atoms_bfactor)])
            outputfile.write('NonHinge'+'\t'+str(len(non_hinge_atoms_bfactor))+'\t'+str(numpy.min(non_hinge_atoms_bfactor))+'\t'+str(numpy.max(non_hinge_atoms_bfactor))+'\t'+str(numpy.mean(non_hinge_atoms_bfactor))+'\t'+str(mode(non_hinge_atoms_bfactor)[0][0])+'\t'+str(numpy.median(non_hinge_atoms_bfactor))+'\t'+str(numpy.std(non_hinge_atoms_bfactor))+'\n')
            return_stats.append(['NonHinge',len(non_hinge_atoms_bfactor),numpy.min(non_hinge_atoms_bfactor),numpy.max(non_hinge_atoms_bfactor),numpy.mean(non_hinge_atoms_bfactor),mode(non_hinge_atoms_bfactor)[0][0],numpy.median(non_hinge_atoms_bfactor),numpy.std(non_hinge_atoms_bfactor)])
            
            p_value = permutation_test(hinge_atoms_bfactor, non_hinge_atoms_bfactor,method='approximate',num_rounds=10000,seed=0)
            outputfile.write('\np-value:\t'+str(p_value)+'\n')
            return p_value,return_stats
        
        def get_leastsquareplane(atoms,HingeResidues,HingePlane):
            """This sub-function gives the Least Square Plane equation of the given atoms.

            This plane is currently gives us an idea about the possible direction of the movement of the hinge residues (See figures in the publication)
            
            Todo:
                * Change the name of the parameter HingePlane (It is misleading)

            Notes:
                * Function level: 1 (1 being top)
                * Possibly will be moved to geometry module later.
                * Check the output file for this equation. The instructions to visualize this plane are given at the end of the file.

            Args:
                atoms ([packman.molecule.Atom])           : Set of atoms. (Read parent method description)
                HingeResidues ([packman.molecule.Residue]): Predicted hinge residues. 
                HingePlane (int)                          : Index of the HingePlane
            
            Returns:
                [a,b,c,d]: Four coefficients essential to define the plane in 3D space.
            """
            hinge_points=[item.get_location() for sublist in [i.get_backbone() for i in HingeResidues] for item in sublist]
            HingeAtoms=[item for sublist in [i.get_backbone() for i in HingeResidues] for item in sublist]
            NonHingeAtoms=[i.get_location() for i in atoms if i not in HingeAtoms]
            
            #Zipping the values
            hinge_xs,hinge_ys,hinge_zs = zip(*hinge_points)
            XS,YS,ZS=zip(*NonHingeAtoms)

            def project_points(x, y, z, a, b, c):
                """Project the points on a given plane

                Projects the points with coordinates x, y, z onto the plane
                defined by a*x + b*y + c*z = 1

                Todo:
                    * Explain the returning variables properly

                Notes:
                    * Function level: 2 (1 being top)
                    * Possibly will be moved to the geometry module
                
                Args:
                    x (float): X-coordinate of the point to be projected on the plane.
                    y (float): Y-coordinate of the point to be projected on the plane.
                    z (float): Z-coordinate of the point to be projected on the plane.
                    a (float): X-component of the plane.
                    b (float): Y-component of the plane.
                    c (float): Z-component of the plane.
                
                Returns:
                    Neccesary information to project the points on the given plane.
                """
                vector_norm = a*a + b*b + c*c
                normal_vector = numpy.array([a, b, c]) / numpy.sqrt(vector_norm)
                point_in_plane = numpy.array([a, b, c]) / vector_norm
                points = numpy.column_stack((x, y, z))
                points_from_point_in_plane = points - point_in_plane
                proj_onto_normal_vector = numpy.dot(points_from_point_in_plane,normal_vector)
                proj_onto_plane = (points_from_point_in_plane - proj_onto_normal_vector[:, None]*normal_vector)
                return point_in_plane + proj_onto_plane
            
            def FitPlane(points):
                """Plane fitting algorithm given the points.

                Notes:
                    * Function level: 2 (1 being top)
                
                Args:
                    points ([float]): Two dimentional array of 3D points
                
                Returns:
                    [a,b,c,d] (float) : Numbers essential to define the Least Square Plane equation
                """
                xs,ys,zs = zip(*points)
                def plane(x, y, params):
                    """Plane Object
                    
                    Args:
                        x ([type]): [description]
                        y ([type]): [description]
                        params ([type]): [description]
                    
                    Returns:
                        [float]: Equation of the plane
                    """
                    a = params[0]
                    b = params[1]
                    c = params[2]
                    z = a*x + b*y + c
                    return z

                def error(params, points):
                    """Get the error
                    
                    Args:
                        params ([type]): [description]
                        points ([type]): [description]
                    
                    Returns:
                        [type]: [description]
                    """
                    result = 0
                    for (x,y,z) in points:
                        plane_z = plane(x, y, params)
                        diff = abs(plane_z - z)
                        result += diff**2
                    return result

                def cross(a, b):
                    """Get cross
                    
                    Args:
                        a ([type]): [description]
                        b ([type]): [description]
                    
                    Returns:
                        [type]: [description]
                    """
                    return [a[1]*b[2] - a[2]*b[1],
                            a[2]*b[0] - a[0]*b[2],
                            a[0]*b[1] - a[1]*b[0]]

                fun = functools.partial(error, points=points)
                params0 = [0, 0, 0]
                res = minimize(fun, params0)

                a = res.x[0]
                b = res.x[1]
                c = res.x[2]

                point  = numpy.array([0.0, 0.0, c])
                normal = numpy.array(cross([1,0,a], [0,1,b]))
                d = -point.dot(normal)

                xx, yy = numpy.meshgrid([min(xs),max(xs)], [min(ys),max(ys)])
                z = (-normal[0] * xx - normal[1] * yy - d) * 1. /normal[2]
                
                #Three points on the plane
                p1=numpy.array([numpy.mean(xx[0]),numpy.mean(yy[:,0]),numpy.mean(z)])
                p2=numpy.array([xx[0][0], yy[0][1], z[0][0]])
                p3=numpy.array([xx[0][1], yy[0][0], z[0][1]])
                p4=numpy.array([xx[1][0], yy[1][1], z[1][0]])
                
                #Two vectors in the plane
                v1=p3-p1
                v2=p2-p1

                #vector normal to the plane
                cp = numpy.cross(v1, v2)
                a, b, c = cp

                # This evaluates a * x3 + b * y3 + c * z3 which equals d
                d = numpy.dot(cp, p3)
                a,b,c,d=a/d,b/d,c/d,d/d
                #outputfile.write('The equation is {0}x + {1}y + {2}z = {3}'.format(a, b, c, d))

                #Plane in PyMol
                outputfile.write('import plane\n')
                outputfile.write('plane.make_plane_points(name=\'HingePlane'+ str(HingePlane) +'\', l1='+str(list(p4))+', l2='+str(list(p2))+', l3='+str(list(p3))+', center=False, makepseudo=False)'+'\n')
                outputfile.write('set cgo_transparency, 0.35, HingePlane'+str(HingePlane)+'\n')
                return a,b,c,d
            
            a,b,c,d=FitPlane(hinge_points)
            projected=project_points(hinge_xs,hinge_ys,hinge_zs,a,b,c)
            return projected
        
        #Alpha shape bit
        alpha_shape, ProteinGraph = AlphaShape( atoms, Alpha, get_graph = True )

        centrality = eccentricity(ProteinGraph)
        centrality_sorted_with_keys = numpy.array([float(centrality[j]) for j in sorted([i for i in centrality.keys()])]).reshape(-1, 1)
        
        #Cluster (4 is like a resolution here)
        km=KMeans(n_clusters=nclusters)
        km.fit(centrality_sorted_with_keys)
        central_nodes=numpy.argwhere(km.labels_==numpy.argmin(km.cluster_centers_)).T[0]
        HingeResidues=list(set([atoms[i].get_parent() for i in central_nodes]))
        HingeResiduesID=[i.get_id() for i in HingeResidues]
        SortedHingeResidues=[x for _,x in sorted(zip(HingeResiduesID,HingeResidues))]
        HingeResiduesID=[i.get_id() for i in SortedHingeResidues]
        
        PredictedHinges=[]
        for i,j in groupby(HingeResiduesID,lambda n, c=count(): n-next(c)):
            Local_Hinge=[k for k in j]
            if(len(Local_Hinge)>MinimumHingeLength):
                PredictedHinges.append(SortedHingeResidues[HingeResiduesID.index(Local_Hinge[0]):HingeResiduesID.index(Local_Hinge[-1])+1])

        #Print part
        Hinges=[]
        Chains=','.join(list(set([i.get_parent().get_id() for i in SortedHingeResidues])))
        outputfile.write('Filename= '+str(filename)+'\t| Chain(s)= '+Chains+'\t| AlphaValue= '+str(Alpha)+'\t| MinimumHingeLength= '+str(MinimumHingeLength)+'\t| EccentricityClusters= '+str(nclusters)+ '\n')
        outputfile.write('Hindge Residues(Predicted):\n')
        for numi,i in enumerate(PredictedHinges):
            #Assigning domain IDs
            for _ in i:_.set_domain_id('FL'+str(numi))
            #Molecule.Hinge(numi,elements,stats,p)
            outputfile.write('\nHinge #'+str(numi+1)+'\nResidues: '+i[0].get_name()+'-'+str(i[0].get_id()) +' to '+i[-1].get_name()+'-'+str(i[-1].get_id())+'\n')
            p_value,hstats=get_statistics(atoms,i,filename=filename)
            outputfile.write('\nPymol Terminal Commands for Visualizing:\ncolor blue, resi '+str(i[0].get_id())+':'+str(i[-1].get_id())+'\n')
            get_leastsquareplane(atoms,i,numi+1)
            outputfile.write("#--------------------------------------------------#\n")
            #HingeObject
            Hinges.append( Hinge(numi,Alpha,i,hstats,p_value))
        
        #Assigning domain IDs
        AllChainResidues = [i for i in HingeResidues[0].get_parent().get_residues()]
        #If sorting is needed
        DomainNumber, flag = 0, True
        for i in AllChainResidues:
            if(i.get_domain_id() == None):
                i.set_domain_id('DM'+str(DomainNumber))
                flag=True
            elif(i.get_domain_id()[:2] == 'FL' and flag):
                DomainNumber+=1
                flag=False
        
        AllChainResidues[0].get_parent().set_hinges(Hinges)
        return alpha_shape

def hinge_cli(args,mol):
    """Command-line Interface for the 'hinge' command. Please check the packman.bin.PACKMAN file for more details.

    This function is for the CLI and not an integral function for the API.

    Args:
        args (parser.parse_args())     : The arguments that were passed by the user to the PACKMAN-hinge app.
        mol (packman.molecule.Protein) : The 'Protein' object for the anaylsis.
    """
    
    def print_footnotes(outputfile):
        """Add footnotes to the output file.
            
        Args:
                outputfile (file): The file to which the footnotes will be written.
        """
        outputfile.write('Footnotes:\n\nSTATISTICS Section Legend:\nN: Number of residues\nMin: Minimum B-factor value\nMax: Maximum B-factor value\nMean: Mean B-factor value\nMode: Mode B-factor value\nMedian: Median B-factor value\nSTDDev: Standard Deviation of the B-factor values\n')
        outputfile.write('\nRMSF Plane Visualization:\nDownload plane.py file from https://pymolwiki.org/index.php/Plane_Wizard and place it in the PyMol working directory if the RMSF plane needs to be visualized.\n')
        return True

    try:
        if(args.chain):
            Backbone = [item for sublist in mol[0][args.chain].get_backbone() for item in sublist]
            alpha_shape= predict_hinge(Backbone, args.outputfile, Alpha=float(args.alpha), filename=args.filename,nclusters=args.e_clusters,MinimumHingeLength=args.minhnglen)
            if(args.generateobj is not None):
                WriteOBJ(Backbone, alpha_shape, args.generateobj)
        else:
            for i in mol[0].get_chains():
                try:
                    Backbone = [item for sublist in mol[0][i.get_id()].get_backbone() for item in sublist]
                    alpha_shape = predict_hinge(Backbone, args.outputfile, Alpha=float(args.alpha), filename=args.filename,nclusters=args.e_clusters,MinimumHingeLength=args.minhnglen)
                except:
                    logging.warning('Chain '+str(i.get_id())+' has none/missing backbone atoms to calculate the hinge.')
    finally:
        print_footnotes(args.outputfile)
        args.outputfile.flush()
        args.logfile.flush()

        if(args.generateobj is not None):
            args.generateobj.flush()

    return True