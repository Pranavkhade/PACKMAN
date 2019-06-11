#!/usr/bin/env python
'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

import numpy
import Molecule

import operator
import argparse
import os
import sys
import functools
import urllib

import networkx as nx

import scipy.optimize
from scipy.spatial import KDTree
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from scipy import stats

from sklearn.cluster import KMeans
from mlxtend.evaluate import permutation_test

from itertools import groupby, count

'''
##################################################################################################
#                                       PACKMAN Functions                                        #
##################################################################################################
'''


def HingePredict(atoms, outputfile, Alpha=float('Inf'),method='AlphaShape',GenerateKirchoff=False,filename='Output.pdb',MinimumHingeLength=5,nclusters=4):
    """
    """
    if(method=='AlphaShape'):
        def GetCircumsphere(Tetrahydron):
            alpha_mat,gamma_mat,Dx_mat,Dy_mat,Dz_mat=[],[],[],[],[]
            for i in Tetrahydron:
                temp_coords=i.get_location()
                alpha_mat.append([temp_coords[0],temp_coords[1],temp_coords[2],1])
                gamma_mat.append([temp_coords[0]**2+temp_coords[1]**2+temp_coords[2]**2,temp_coords[0],temp_coords[1],temp_coords[2]])
                Dx_mat.append([temp_coords[0]**2+temp_coords[1]**2+temp_coords[2]**2,temp_coords[1],temp_coords[2],1])
                Dy_mat.append([temp_coords[0]**2+temp_coords[1]**2+temp_coords[2]**2,temp_coords[0],temp_coords[2],1])
                Dz_mat.append([temp_coords[0]**2+temp_coords[1]**2+temp_coords[2]**2,temp_coords[0],temp_coords[1],1])
            alpha=numpy.linalg.det(alpha_mat)
            gamma=numpy.linalg.det(gamma_mat)
            Dx=(+numpy.linalg.det(Dx_mat))
            Dy=(-numpy.linalg.det(Dy_mat))
            Dz=(+numpy.linalg.det(Dz_mat))
            Centre=numpy.array([Dx/2*alpha,Dy/2*alpha,Dz/2*alpha])
            Radius=numpy.sqrt(Dx**2 + Dy**2 + Dz**2 - 4*alpha*gamma)/(2*numpy.absolute(alpha))
            DistanceToSurface=numpy.linalg.norm(Centre-AllAtoms[a].get_location())
            return Centre,Radius,DistanceToSurface

        def AlphaTest(atoms,Alpha,Centre,CircumRadius,DistanceToSurface,Tree=None):
            """
            NOTE: Confirm is sphere hollowness should be checked or it is covered by tesselations.
            """
            if(CircumRadius<Alpha):
                return True
            else:
                return False
        
        def GetStats(atoms,SelectedHingeResidues,filename='Output'):
            hinge_atoms=[i.get_backbone() for i in SelectedHingeResidues]
            hinge_atoms=[item for sublist in hinge_atoms for item in sublist]
            non_hinge_atoms=list(set([i for i in atoms])-set(hinge_atoms))
            all_atoms_bfactor=[i.get_bfactor() for i in atoms]
            hinge_atoms_bfactor=[i.get_bfactor() for i in hinge_atoms]
            non_hinge_atoms_bfactor=[i.get_bfactor() for i in non_hinge_atoms]

            #Example\tGroup\tNumberofelements\tMin\tMax\tMean\tMode\tMedian\tStandardDeviation\n

            outputfile.write('\nSTATISTICS\n\t\tN\tMin\tMax\tMean\tMode\tMedian\tSTDDev\n')
            outputfile.write('Total   '+'\t'+str(len(all_atoms_bfactor))+'\t'+str(numpy.min(all_atoms_bfactor))+'\t'+str(numpy.max(all_atoms_bfactor))+'\t'+str(numpy.mean(all_atoms_bfactor))+'\t'+str(stats.mode(all_atoms_bfactor)[0][0])+'\t'+str(numpy.median(all_atoms_bfactor))+'\t'+str(numpy.std(all_atoms_bfactor))+'\n')
            outputfile.write('Hinge   '+'\t'+str(len(hinge_atoms_bfactor))+'\t'+str(numpy.min(hinge_atoms_bfactor))+'\t'+str(numpy.max(hinge_atoms_bfactor))+'\t'+str(numpy.mean(hinge_atoms_bfactor))+'\t'+str(stats.mode(hinge_atoms_bfactor)[0][0])+'\t'+str(numpy.median(hinge_atoms_bfactor))+'\t'+str(numpy.std(hinge_atoms_bfactor))+'\n')
            outputfile.write('NonHinge'+'\t'+str(len(non_hinge_atoms_bfactor))+'\t'+str(numpy.min(non_hinge_atoms_bfactor))+'\t'+str(numpy.max(non_hinge_atoms_bfactor))+'\t'+str(numpy.mean(non_hinge_atoms_bfactor))+'\t'+str(stats.mode(non_hinge_atoms_bfactor)[0][0])+'\t'+str(numpy.median(non_hinge_atoms_bfactor))+'\t'+str(numpy.std(non_hinge_atoms_bfactor))+'\n')
            
            p_value = permutation_test(hinge_atoms_bfactor, non_hinge_atoms_bfactor,method='approximate',num_rounds=10000,seed=0)
            outputfile.write('\np-value:\t'+str(p_value)+'\n')
            return p_value
        
        

        def GetLeastSquarePlane(atoms,HingeResidues,HingePlane):
            """
            """
            hinge_points=[item.get_location() for sublist in [i.get_backbone() for i in HingeResidues] for item in sublist]
            HingeAtoms=[item for sublist in [i.get_backbone() for i in HingeResidues] for item in sublist]
            NonHingeAtoms=[i.get_location() for i in atoms if i not in HingeAtoms]
            
            #Zipping the values
            hinge_xs,hinge_ys,hinge_zs = zip(*hinge_points)
            XS,YS,ZS=zip(*NonHingeAtoms)

            def project_points(x, y, z, a, b, c):
                """
                Projects the points with coordinates x, y, z onto the plane
                defined by a*x + b*y + c*z = 1
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
                """
                """
                xs,ys,zs = zip(*points)
                def plane(x, y, params):
                    a = params[0]
                    b = params[1]
                    c = params[2]
                    z = a*x + b*y + c
                    return z

                def error(params, points):
                    result = 0
                    for (x,y,z) in points:
                        plane_z = plane(x, y, params)
                        diff = abs(plane_z - z)
                        result += diff**2
                    return result

                def cross(a, b):
                    return [a[1]*b[2] - a[2]*b[1],
                            a[2]*b[0] - a[0]*b[2],
                            a[0]*b[1] - a[1]*b[0]]

                fun = functools.partial(error, points=points)
                params0 = [0, 0, 0]
                res = scipy.optimize.minimize(fun, params0)

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
        
        AllAtoms=[i for i in atoms]
        AlphaKirchoff=numpy.zeros((len(AllAtoms),len(AllAtoms)))
        DelaunayTesssellation=Delaunay([i.get_location() for i in AllAtoms])
        SelectedTesselations=[]
        ProteinGraph=nx.Graph()
        if(GenerateKirchoff):
            Tree=KDTree([i.get_location() for i in AllAtoms])
            for a,b,c,d in DelaunayTesssellation.vertices:
                Tetrahydron=[AllAtoms[a],AllAtoms[b],AllAtoms[c],AllAtoms[d]]
                Centre,Radius,DistanceToSurface=GetCircumsphere(Tetrahydron)
                #Alpha Test
                if(AlphaTest(AllAtoms,Alpha,Centre,Radius,DistanceToSurface,Tree=Tree)):
                    ProteinGraph.add_nodes_from([a,b,c,d])
                    ProteinGraph.add_edges_from([(a,b),(a,c),(a,d),(b,c),(b,d),(c,d)])
                    SelectedTesselations.append([AllAtoms[a],AllAtoms[b],AllAtoms[c],AllAtoms[d]])
                    AlphaKirchoff[a][b]=1
                    AlphaKirchoff[a][c]=1
                    AlphaKirchoff[a][d]=1
                    AlphaKirchoff[b][a]=1
                    AlphaKirchoff[b][c]=1
                    AlphaKirchoff[b][d]=1
                    AlphaKirchoff[c][a]=1
                    AlphaKirchoff[c][b]=1
                    AlphaKirchoff[c][d]=1
                    AlphaKirchoff[d][a]=1
                    AlphaKirchoff[d][b]=1
                    AlphaKirchoff[d][c]=1

            for numi,i in enumerate(AlphaKirchoff):
                AlphaKirchoff[numi][numi]=numpy.sum(i)
        else:
            for a,b,c,d in DelaunayTesssellation.vertices:
                Tetrahydron=[AllAtoms[a],AllAtoms[b],AllAtoms[c],AllAtoms[d]]
                Centre,Radius,DistanceToSurface=GetCircumsphere(Tetrahydron)
                #Alpha Test
                if(AlphaTest(AllAtoms,Alpha,Centre,Radius,DistanceToSurface)):
                    ProteinGraph.add_nodes_from([a,b,c,d])
                    ProteinGraph.add_edges_from([(a,b),(a,c),(a,d),(b,c),(b,d),(c,d)])
                    SelectedTesselations.append([AllAtoms[a],AllAtoms[b],AllAtoms[c],AllAtoms[d]])

        centrality=nx.eccentricity(ProteinGraph)
        centrality_sorted_with_keys=numpy.array([float(centrality[j]) for j in sorted([i for i in centrality.keys()])]).reshape(-1, 1)
        #Cluster (4 is like a resolution here)
        km=KMeans(n_clusters=nclusters)
        km.fit(centrality_sorted_with_keys)
        central_nodes=numpy.argwhere(km.labels_==numpy.argmin(km.cluster_centers_)).T[0]
        HingeResidues=list(set([atoms[i].get_parent() for i in central_nodes]))
        SortedHingeResidues=sorted(HingeResidues)
        HingeResiduesID=[i.get_id() for i in sorted(HingeResidues)]

        PredictedHinges=[]
        for i,j in groupby(HingeResiduesID,lambda n, c=count(): n-next(c)):
            Hinge=[k for k in j]
            if(len(Hinge)>MinimumHingeLength):
                PredictedHinges.append(SortedHingeResidues[HingeResiduesID.index(Hinge[0]):HingeResiduesID.index(Hinge[-1])+1])

        #Print part
        outputfile.write('Filename= '+str(filename)+'\t| AlphaValue= '+str(Alpha)+'\t| MinimumHingeLength= '+str(MinimumHingeLength)+'\t| EccentricityClusters= '+str(nclusters)+ '\n')
        outputfile.write('Hindge Residues(Predicted):\n')
        for numi,i in enumerate(PredictedHinges):
            #outputfile.write('\nHinge #'+str(numi+1)+'\nResidues: '+','.join([j.get_name()+'-'+str(j.get_id()) for j in i])+'\n')
            outputfile.write('\nHinge #'+str(numi+1)+'\nResidues: '+i[0].get_name()+'-'+str(i[0].get_id()) +' to '+i[-1].get_name()+'-'+str(i[-1].get_id())+'\n')
            GetStats(atoms,i,filename=filename)
            #outputfile.write('\nPymol Terminal Commands for Visualizing:\ncolor blue, resi '+'+'.join(str(j.get_id()) for j in i)+'\n')
            outputfile.write('\nPymol Terminal Commands for Visualizing:\ncolor blue, resi '+str(i[0].get_id())+':'+str(i[-1].get_id())+'\n')
            GetLeastSquarePlane(atoms,i,numi+1)
            outputfile.write("#--------------------------------------------------#\n")
        #RotateHingeBonds(atoms,HingeResidues)

        if(GenerateKirchoff):
            return AlphaKirchoff,SelectedTesselations,PredictedHinges
        else:
            return SelectedTesselations,PredictedHinges


'''
##################################################################################################
#                                    Non Algorithm Functions                                     #
##################################################################################################
'''

def WriteOBJ(atoms,faces, fh):
    """
    """
    NewIDs={i.get_id():numi+1 for numi,i in enumerate(atoms)}
    fh.write('mtllib master.mtl\ng\n')
    fh.write('usemtl atoms\n')
    for i in atoms:
        x,y,z=i.get_location()
        fh.write("v %f %f %f\n"%(x,y,z))
    
    line='usemtl bonds\nl'
    for i in atoms:
        line=line+" "+str(NewIDs[i.get_id()])
    line=line+'\n'
    fh.write(line)
    
    fh.write('usemtl faces\n')
    for i in faces:
        faces=[NewIDs[j.get_id()] for j in i]
        fh.write("f %i %i %i %i\n"%(faces[0],faces[1],faces[2],faces[3]))
        #fh.write("l %i %i %i %i\n"%(faces[0],faces[1],faces[2],faces[3]))
    return True


def GenerateMovie(atoms,ProteinGraph,SelectedHingeResidues,gamma=1.0,dr=15.0,pf=2,mode_number=6):
    """  
    """
    #Remove the hinge from the graph
    SelectedHingeResiduesID=[i.get_id() for i in SelectedHingeResidues]
    hinge_atoms_index=[i.get_backbone() for i in SelectedHingeResidues]
    hinge_atoms_index=[atoms.index(item) for sublist in hinge_atoms_index for item in sublist]
    
    ProteinGraphHingeRemoved=ProteinGraph.copy()
    ProteinGraphHingeRemoved.remove_nodes_from(hinge_atoms_index)

    coords=numpy.array([i.get_location() for i in atoms])
    n_atoms=len(coords)
    #print "Number of Atoms",n_atoms
    hessian=numpy.zeros((n_atoms*3, n_atoms*3), float)
    distance_mat=numpy.ones((n_atoms*3, n_atoms*3), float)
    for i in range(len(coords)):
        diff=coords[i+1:, :]-coords[i]
        squared_diff = diff*diff
        for j, s_ij in enumerate(squared_diff.sum(1)):
            if s_ij <= dr**2:
                diff_coords = diff[j]
                j = j + i + 1
                if(atoms[i].get_parent().get_id() in SelectedHingeResiduesID and atoms[j].get_parent().get_id() in SelectedHingeResiduesID):
                    #Both atoms in the predicted hinge
                    gamma=0.5
                elif(atoms[i].get_parent().get_id() not in SelectedHingeResiduesID and atoms[j].get_parent().get_id() not in SelectedHingeResiduesID):
                    #Both atoms NOT in the predicted hinge
                    try:
                        if(nx.has_path(ProteinGraphHingeRemoved,i,j)):
                            gamma=1.0
                        else:
                            gamma=0.1
                    except:
                        gamma=0.75
                        '''networkx.exception.NodeNotFound'''
                else:
                    gamma=1.0
                derivative = numpy.outer(diff_coords, diff_coords)*(float(-gamma)/s_ij)
                hessian[i*3:i*3+3, j*3:j*3+3] = derivative
                hessian[j*3:j*3+3, i*3:i*3+3] = derivative
                hessian[i*3:i*3+3, i*3:i*3+3] = hessian[i*3:i*3+3, i*3:i*3+3] - derivative
                hessian[j*3:j*3+3, j*3:j*3+3] = hessian[j*3:j*3+3, j*3:j*3+3] - derivative
                
                d = numpy.sqrt(s_ij)
                lobj = [[d,d,d],[d,d,d], [d,d,d]]
                dmat = numpy.array(lobj)
                distance_mat[i*3:i*3+3, j*3:j*3+3] = dmat
    if pf != None:
        hessian = hessian/(distance_mat**pf)
    eigenvalue,eigenvector = numpy.linalg.eigh(hessian)

    x0=numpy.array([i.get_location() for i in atoms])
    new_coords=[]
    scale=0.1
    n=50
    print eigenvalue
    with open('movie.pdb','w') as fh:
        for j in range(n)+range(n)[::-1]:
            #fh.write('MODEL        ')
            for numi,i in enumerate(x0):
                #print(len(self.eigenvector[mode_number]))
                #print(len(x0))
                #new_x=i[0]+scale*i[0]*j*eigenvector[mode_number][(numi*3)+0]
                #new_y=i[1]+scale*i[1]*j*eigenvector[mode_number][(numi*3)+1]
                #new_z=i[2]+scale*i[2]*j*eigenvector[mode_number][(numi*3)+2]
                new_x=i[0]+scale*i[0]*j*eigenvector[(numi*3)+0][mode_number]
                new_y=i[1]+scale*i[1]*j*eigenvector[(numi*3)+1][mode_number]
                new_z=i[2]+scale*i[2]*j*eigenvector[(numi*3)+2][mode_number]
                
                new_coords.append([new_x,new_y,new_z])
                fh.write("ATOM  %5s  %-4s%3s %1s%4s    %8s%8s%8s%5s%6s      %-4s%2s%2s\n"%(atoms[numi].get_id(),atoms[numi].get_name(),atoms[numi].get_parent().get_name(),atoms[numi].get_parent().get_parent().get_id(),atoms[numi].get_parent().get_id(),round(new_x,3),round(new_y,3),round(new_z,3),atoms[numi].get_occupancy(),atoms[numi].get_bfactor(),'',atoms[numi].get_elment(),''))
            fh.write('TER    %4s      %3s %1s%4s                                                      \n'%(atoms[-1].get_id(),atoms[-1].get_parent().get_name(),atoms[-1].get_parent().get_parent().get_id(),atoms[-1].get_parent().get_id()))
            fh.write('ENDMDL\n')
    return True

def RotateHingeBonds(atoms,PredictedHinges):
    """
    """
    print ConvexHull([i.get_location() for i in atoms]).area
    def CalculateDihedral(atom1,atom2,atom3,atom4):
        p0=atom1.get_location()
        p1=atom2.get_location()
        p2=atom3.get_location()
        p3=atom4.get_location()
        b0 = -1.0*(p1 - p0)
        b1 = p2 - p1
        b2 = p3 - p2
        b1 /= numpy.linalg.norm(b1)
        v = b0 - numpy.dot(b0, b1)*b1
        w = b2 - numpy.dot(b2, b1)*b1
        x = numpy.dot(v, w)
        y = numpy.dot(numpy.cross(b1, v), w)
        radang=numpy.arctan2(y, x)
        return radang

    for i in PredictedHinges:
        Backbone=[item for sublist in [j.get_backbone() for j in i] for item in sublist if item.get_name()!='O']
        psi=[]
        phi=[]
        for numj,j in enumerate(Backbone):
            try:    
                #print Backbone[numj].get_name()+'-'+Backbone[numj+1].get_name()+'-'+Backbone[numj+2].get_name()+'-'+Backbone[numj+3].get_name()
                if(Backbone[numj+1].get_name()+'-'+Backbone[numj+2].get_name()=='N-CA'):
                    phi.append(CalculateDihedral(Backbone[numj],Backbone[numj+1],Backbone[numj+2],Backbone[numj+3]))
                if(Backbone[numj+1].get_name()+'-'+Backbone[numj+2].get_name()=='CA-C'):
                    psi.append(CalculateDihedral(Backbone[numj],Backbone[numj+1],Backbone[numj+2],Backbone[numj+3]))
            except IndexError:
                break
        print 'PSI:',numpy.mean(psi)
        print 'PHI:',numpy.mean(phi),'\n<<<<'
    return True

'''
##################################################################################################
#                                          Interface                                             #
##################################################################################################
'''

def IO():
    """
    INFO: Argument parser to the program for now.
    """
    parser=argparse.ArgumentParser(description='PACKMAN: PACKing and Motion ANalysis. (https://github.com/Pranavkhade/PACKMAN)')

    parser.add_argument('-pdbid','--pdbid', metavar='PDB_ID', type=str, help='If provided, the PBD with this ID will be downloaded and saved to FILENAME.')

    parser.add_argument('alpha', metavar='AlphaValue', help='Recommended: 2.8 for closed; 4.5 for open form, Please refer to the paper for more details')
    parser.add_argument('filename', metavar='FILENAME', help='Path and filename of the PDB file.')

    parser.add_argument('--e_clusters',metavar='NumberOfEccentricityClusters',type=int,default=4,help='Recommended: 4, Please refer to the paper for more details')
    parser.add_argument('--minhnglen',metavar='MinimumHingeLength',type=int,default=5,help='Recommended: 5, Please refer to the paper for more details')
    parser.add_argument("--chain", help='Enter The Chain ID')
    parser.add_argument('--generateobj', type=argparse.FileType('wb', 0), help='Path and filename to save the .obj file at. Ignored unless --chain is provided.')

    # web server parameters
    web_server_group = parser.add_argument_group('Web server parameters', 'Used by the web form')
    web_server_group.add_argument('--outputfile', type=argparse.FileType('wb', 0), default=sys.stdout, help='Path and filename write output to')
    web_server_group.add_argument('--logfile', type=argparse.FileType('wb', 0), default=sys.stderr, help='Path and filename write log messages to')
    web_server_group.add_argument('--callbackurl', type=str, help='Optional callback url if this script was called from Drupal.')
    web_server_group.add_argument('--nodeid', type=int, help='Optional node id if this script was called from Drupal.')
    
    args=parser.parse_args()
    return args

'''
##################################################################################################
#                                             Main                                               #
##################################################################################################
'''

def main():
    """
    """
    args=IO()

    if(args.pdbid is not None):
        Molecule.DownloadPDB(args.pdbid, args.filename)

    mol = Molecule.LoadPDB(args.filename)

    try:

        if(args.chain):
            Backbone = [item for sublist in mol[0][args.chain].get_backbone() for item in sublist]
            SelectedTesselations, PredictedHinges = HingePredict(Backbone, args.outputfile, Alpha=float(args.alpha), filename=args.filename,nclusters=args.e_clusters,MinimumHingeLength=args.minhnglen)
            if(args.generateobj is not None):
                WriteOBJ(Backbone, SelectedTesselations, args.generateobj)
        else:
            for i in mol[0].get_chains():
                Backbone = [item for sublist in mol[0].get_backbone() for item in sublist]
                SelectedTesselations, PredictedHinges = HingePredict(Backbone, args.outputfile, Alpha=float(args.alpha), filename=args.filename,nclusters=args.e_clusters,MinimumHingeLength=args.minhnglen)
        
        if args.nodeid is not None:
            urllib.urlopen(args.callbackurl + '/' + str(args.nodeid) + "/0")

    except Exception:

        if args.nodeid is not None:
            urllib.urlopen(args.callbackurl + '/' + str(args.nodeid) + "/1")
    
    finally:

        print_footnotes(args.outputfile)
        args.outputfile.flush()
        args.logfile.flush()

        if(args.generateobj is not None):
            args.generateobj.flush()
            
    return True

'''
##################################################################################################
#                                             End                                                #
##################################################################################################
'''
def print_footnotes(outputfile):
    outputfile.write('Footnotes:\n\nSTATISTICS Section Legend:\nN: Number of residues\nMin: Minimum B-factor value\nMax: Maximum B-factor value\nMean: Mean B-factor value\nMode: Mode B-factor value\nMedian: Median B-factor value\nSTDDev: Standard Deviation of the B-factor values\n')
    outputfile.write('\nRMSF Plane Visualization:\nDownload plane.py file from https://pymolwiki.org/index.php/Plane_Wizard and place it in the PyMol working directory if the RMSF plane needs to be visualized.\n')
    return True



if(__name__=='__main__'):
    main()
