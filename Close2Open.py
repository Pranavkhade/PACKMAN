import PACKMAN
import Molecule
import numpy
import copy
import scipy
import functools
import random

from scipy import optimize

def get_leastsquareplane(HingeResidues):
    """
    """
    hinge_points=[item.get_location() for sublist in [i.get_backbone() for i in HingeResidues] for item in sublist]
    HingeAtoms=[item for sublist in [i.get_backbone() for i in HingeResidues] for item in sublist]
    
    #Zipping the values
    hinge_xs,hinge_ys,hinge_zs = zip(*hinge_points)

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
        return a,b,c,d
    
    a,b,c,d=FitPlane(hinge_points)
    projected=project_points(hinge_xs,hinge_ys,hinge_zs,a,b,c)
    Difference=0
    for numi,i in enumerate(hinge_points):Difference+=numpy.linalg.norm(hinge_points[numi]-projected[numi])
    return Difference


def rotate_bond(Atom1,Atom2,Atom3,Atom4,model,theta):
    p0=Atom1.get_location()
    p1=Atom2.get_location()
    p2=Atom3.get_location()
    p3=Atom4.get_location()
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= numpy.linalg.norm(b1)
    v = b0 - numpy.dot(b0, b1)*b1
    w = b2 - numpy.dot(b2, b1)*b1
    x = numpy.dot(v, w)
    y = numpy.dot(numpy.cross(b1, v), w)
    radang=numpy.arctan2(y, x)
    rotang=theta-radang
    sn=numpy.sin(rotang)
    cs=numpy.cos(rotang)
    t=1-cs

    v2x=Atom1.get_location()[0]-Atom2.get_location()[0]
    v2y=Atom1.get_location()[1]-Atom2.get_location()[1]
    v2z=Atom1.get_location()[2]-Atom2.get_location()[2]
    #Normalize
    mag=numpy.sqrt(numpy.square(v2x)+numpy.square(v2y)+numpy.square(v2z))
    x = float(v2x)/mag
    y = float(v2y)/mag
    z = float(v2z)/mag
    #set up the rotation matrix
    m=numpy.zeros(9)
    m[0]= t*x*x + cs
    m[1]= t*x*y + sn*z
    m[2]= t*x*z - sn*y
    m[3]= t*x*y - sn*z
    m[4]= t*y*y + cs
    m[5]= t*y*z + sn*x
    m[6]= t*x*z + sn*y
    m[7]= t*y*z - sn*x
    m[8]= t*z*z + cs

    #Rotate The Atoms
    tx = Atom1.get_location()[0]
    ty = Atom1.get_location()[1]
    tz = Atom1.get_location()[2]
    for i in model.get_atoms():
        #Imporve this later on, this skips some atoms
        if(i.get_parent().get_id()>Atom2.get_parent().get_id()):
            Coordinates=i.get_location()
            Coordinates[0]-=tx
            Coordinates[1]-=ty
            Coordinates[2]-=tz
            x=Coordinates[0]*m[0]+Coordinates[1]*m[1]+Coordinates[2]*m[2]
            y=Coordinates[0]*m[3]+Coordinates[1]*m[4]+Coordinates[2]*m[5]
            z=Coordinates[0]*m[6]+Coordinates[1]*m[7]+Coordinates[2]*m[8]
            Coordinates[0]=x
            Coordinates[1]=y
            Coordinates[2]=z
            Coordinates[0]+=tx
            Coordinates[1]+=ty
            Coordinates[2]+=tz
            i.set_location(numpy.around(numpy.array(Coordinates),decimals=3))
    return model

def get_torsion_angles(Backbone):
    omega,phi,psi=[],[],[]
    for i in range(0,len(Backbone)-3):
        combination=Backbone[i].get_name()+'-'+Backbone[i+1].get_name()+'-'+Backbone[i+2].get_name()+'-'+Backbone[i+3].get_name()
        if(combination=='CA-C-N-CA'):
            omega.append([Backbone[i],Backbone[i+1],Backbone[i+2],Backbone[i+3]])
        elif(combination=='C-N-CA-C'):
            phi.append([Backbone[i],Backbone[i+1],Backbone[i+2],Backbone[i+3]])
        elif(combination=='N-CA-C-N'):
            psi.append([Backbone[i],Backbone[i+1],Backbone[i+2],Backbone[i+3]])
    return omega,phi,psi


#GLOBAL VARIABLES
lowest_diff=float('Inf')
models=[]
native_clashes=None

def SetAngles(features,model,angles):
    omega,phi,psi=angles
    phipsi=phi+psi
    for numi,i in enumerate(phipsi):
        rotate_bond(i[0],i[1],i[2],i[3],model,features[numi])
        #print('angle:',features[numi])
    Difference=get_leastsquareplane(model['A'].get_hinges()[0].get_elements())
    global lowest_diff,native_clashes
    if(model.check_clashes()>native_clashes):
        Difference=Difference*2
    if(Difference<lowest_diff):
        models.append(copy.deepcopy(model))
        lowest_diff=Difference
    print(Difference)
    return Difference


def Close2Open(model,chain='A',hinge_number=0):
    HingeBackbone=[item for sublist in [i.get_backbone() for i in model[chain].get_hinges()[0].get_elements()] for item in sublist if item.get_name()!='O']
    omega,phi,psi=get_torsion_angles(HingeBackbone)

    numberofparameters=len(phi)+len(psi)
    RandomFeatureValues=[]
    for _ in range(0,numberofparameters):RandomFeatureValues.append(random.randrange(-180,180,10))
    ranges=list(zip([-180]*numberofparameters,[180]*numberofparameters))
    result=optimize.minimize(SetAngles,x0=RandomFeatureValues,args=(model,[omega,phi,psi]),bounds=ranges,method='L-BFGS-B')
    #,options={'maxiter':1}

    new_mol=Molecule.Protein('test','testname',models)
    Molecule.WritePDB(new_mol,'test.pdb')
    return True


def main():
    mol=Molecule.LoadPDB('1prw.pdb')
    Backbone = [item for sublist in mol[0].get_backbone() for item in sublist]
    SelectedTesselations=PACKMAN.HingePredict(Backbone,open('test.txt','w'),Alpha=2.8)
    global native_clashes
    native_clashes=mol[0].check_clashes()
    print(native_clashes)
    Close2Open(mol[0])
    return True


if(__name__=='__main__'):
    main()