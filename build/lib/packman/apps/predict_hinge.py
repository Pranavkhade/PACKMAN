'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

import numpy
import functools

import networkx as nx

from scipy.spatial import Delaunay
from scipy.spatial import KDTree
from scipy.stats import mode
from scipy.optimize import minimize

from sklearn.cluster import KMeans
from mlxtend.evaluate import permutation_test
from itertools import groupby, count

from ..molecule import Hinge



def predict_hinge(atoms, outputfile, Alpha=float('Inf'),method='AlphaShape',GenerateKirchoff=False,filename='Output.pdb',MinimumHingeLength=5,nclusters=4):
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
            '''
            '''
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
            p_value,hstats=GetStats(atoms,i,filename=filename)
            outputfile.write('\nPymol Terminal Commands for Visualizing:\ncolor blue, resi '+str(i[0].get_id())+':'+str(i[-1].get_id())+'\n')
            GetLeastSquarePlane(atoms,i,numi+1)
            outputfile.write("#--------------------------------------------------#\n")
            #HingeObject
            Hinges.append( Hinge(numi,i,hstats,p_value))
        
        #Assigning domain IDs
        AllChainResidues=[i for i in HingeResidues[0].get_parent().get_residues()]
        #If sorting is needed
        #AllChainResiduesID=[i.get_id() for i in AllChainResidues]
        #print([i.get_id() for i in AllChainResidues])
        #AllChainResidues=[x for _,x in sorted(zip(AllChainResiduesID,AllChainResidues))]
        DomainNumber,flag=0,True
        for i in AllChainResidues:
            if(i.get_domain_id()==None):
                i.set_domain_id('DM'+str(DomainNumber))
                flag=True
            elif(i.get_domain_id()[:2]=='FL' and flag):
                DomainNumber+=1
                flag=False
        
        AllChainResidues[0].get_parent().set_hinges(Hinges)

        if(GenerateKirchoff):
            return AlphaKirchoff,SelectedTesselations
        else:
            return SelectedTesselations
