'''
Author: Pranav Khade(pranavk@iastate.edu)
SourceofInformation: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
'''
import numpy
import argparse
import os
import warnings

'''
##################################################################################################
#                                     Molecule Objects                                           #
##################################################################################################
'''

class Atom():
    """
    This class contains the information about the Atom.
    """
    def __init__(self,id,AtomName,Coordinates,Occupancy,bfactor,Element,Charge,parent):
        self.__id=id
        self.__AtomName=AtomName
        self.__AlternateLocationIndicator=None #remove it later
        self.__parent=parent
        self.__Coordinates=Coordinates
        self.__Occupancy=Occupancy
        self.__bfactor=bfactor
        self.__SegmentIdentifier=None
        self.__Element=Element

    #Get Functions
    def get_id(self):
        """
        :returns: ID of an 'Atom' object
        """
        return self.__id
    
    def get_name(self):
        """
        :returns: Name of an 'Atom' object
        """
        return self.__AtomName
    
    def get_alternatelocationindicator(self):
        """
        :returns: Alternate Location Indicator
        """
        return self.__AlternateLocationIndicator

    def get_parent(self):
        """
        :returns: Parent residue of an 'Atom' Object. (If available)
        """
        return self.__parent
    
    def get_location(self):
        """
        :returns: Cartesian Coordinates of an 'Atom'
        """
        return self.__Coordinates
    
    def get_occupancy(self):
        """
        :returns: Occupancy of an 'Atom'
        """
        return self.__Occupancy
    
    def get_bfactor(self):
        """
        :returns: B-Factor value of an 'Atom'
        """
        return self.__bfactor
    
    def get_segmentidentifier(self):
        """
        :returns: Segment Identifier of an 'Atom'
        """
        return self.__SegmentIdentifier
    
    def get_element(self):
        """
        :returns: Element of an 'Atom'
        """
        return self.__Element
    
    def get_domain_id(self):
        """
        :returns: Domain ID of the parent Residue
        """
        return self.get_parent().get_domain_id()

    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    def set_name(self,new_name):
        """
        :param new_name: New name to be assigned
        """
        self.__AtomName=new_name
    
    def set_alternatelocationindicator(self,new_alternatelocationindicator):
        """
        :param new_alternatelocationindicator: New Alternate Location Indicator to be assigned
        """
        self.__AlternateLocationIndicator=new_alternatelocationindicator

    def set_parent(self,new_parent):
        """
        :param new_parent: New parent to be assigned
        """
        self.__parent=new_parent
    
    def set_location(self,new_location):
        """
        :param new_location: New Cartesian Coordinates to be assigned
        """
        self.__Coordinates=new_location
    
    def set_occupancy(self,new_occupancy):
        """
        :param new_occupancy: New occupancy value to be assigned
        """
        self.__Occupancy=new_occupancy
    
    def set_bfactor(self,new_bfactor):
        """
        :param new_bfactor: New B-Factor to be assigned
        """
        self.__bfactor=new_bfactor
    
    def set_segmentidentifier(self,new_segmentidentifier):
        """
        :param new_segmentidentifier: New Segment Identifier to be assigned
        """
        self.__SegmentIdentifier=new_segmentidentifier
    
    def set_elment(self,new_element):
        """
        :param new_element: New element to be assigned
        """
        self.__Element=new_element

    #Calculation Functions
    def calc_dist(self,another_atom):
        """
        :param another_atom: 'Atom' object to which the distance should be calculated.
        """
        return numpy.linalg.norm(self.__Coordinates-another_atom.get_location())

class HetAtom():
    """
    This class contains the information about the HetAtom.
    """
    def __init__(self,id,AtomName,Coordinates,Occupancy,bfactor,Element,Charge,parent):
        self.__id=id
        self.__AtomName=AtomName
        self.__AlternateLocationIndicator=None #remove it later
        self.__parent=parent
        self.__Coordinates=Coordinates
        self.__Occupancy=Occupancy
        self.__bfactor=bfactor
        self.__SegmentIdentifier=None
        self.__Element=Element

    #Get Functions
    def get_id(self):
        """
        :returns: ID of an 'Atom' object
        """
        return self.__id
    
    def get_name(self):
        """
        :returns: Name of an 'Atom' object
        """
        return self.__AtomName
    
    def get_alternatelocationindicator(self):
        """
        :returns: Alternate Location Indicator
        """
        return self.__AlternateLocationIndicator

    def get_parent(self):
        """
        :returns: Parent residue of an 'Atom' Object. (If available)
        """
        return self.__parent
    
    def get_location(self):
        """
        :returns: Cartesian Coordinates of an 'Atom'
        """
        return self.__Coordinates
    
    def get_occupancy(self):
        """
        :returns: Occupancy of an 'Atom'
        """
        return self.__Occupancy
    
    def get_bfactor(self):
        """
        :returns: B-Factor value of an 'Atom'
        """
        return self.__bfactor
    
    def get_segmentidentifier(self):
        """
        :returns: Segment Identifier of an 'Atom'
        """
        return self.__SegmentIdentifier
    
    def get_elment(self):
        """
        :returns: Element of an 'Atom'
        """
        return self.__Element

    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    def set_name(self,new_name):
        """
        :param new_name: New name to be assigned
        """
        self.__AtomName=new_name
    
    def set_alternatelocationindicator(self,new_alternatelocationindicator):
        """
        :param new_alternatelocationindicator: New Alternate Location Indicator to be assigned
        """
        self.__AlternateLocationIndicator=new_alternatelocationindicator

    def set_parent(self,new_parent):
        """
        :param new_parent: New parent to be assigned
        """
        self.__parent=new_parent
    
    def set_location(self,new_location):
        """
        :param new_location: New Cartesian Coordinates to be assigned
        """
        self.__Coordinates=new_location
    
    def set_occupancy(self,new_occupancy):
        """
        :param new_occupancy: New occupancy value to be assigned
        """
        self.__Occupancy=new_occupancy
    
    def set_bfactor(self,new_bfactor):
        """
        :param new_bfactor: New B-Factor to be assigned
        """
        self.__bfactor=new_bfactor
    
    def set_segmentidentifier(self,new_segmentidentifier):
        """
        :param new_segmentidentifier: New Segment Identifier to be assigned
        """
        self.__SegmentIdentifier=new_segmentidentifier
    
    def set_elment(self,new_element):
        """
        :param new_element: New element to be assigned
        """
        self.__Element=new_element

    #Calculation Functions
    def calc_dist(self,another_atom):
        """
        :param another_atom: 'Atom' object to which the distance should be calculated.
        """
        return numpy.linalg.norm(self.__Coordinates-another_atom.get_location())

    

class Residue():
    '''
    Domain Identifiers: (Obtained after running PACKMAN)
    FL: Flexible Linkers (Hinges); Followed by a number
    DM: Domain; Followed by a number
    '''
    def __init__(self,id,name,parent):
        self.__id=id
        self.__name=name
        self.__parent=parent
        self.__Atoms=None
        self.__domain_id=None
    
    def __setitem__(self,id,Atom):
        try:
            self.__Atoms[id]=Atom
        except:
            self.__Atoms={}
            self.__Atoms[id]=Atom
    
    #Get Functions
    def get_id(self):
        """
        :returns: ID of an 'Residue' object
        """
        return self.__id
    
    def get_name(self):
        """
        :returns: Name of an 'Residue' object
        """
        return self.__name

    def get_parent(self):
        """
        :returns: Parent residue of an 'Residue' Object. (If available)
        """
        return self.__parent
    
    def get_atoms(self):
        """
        :returns: Generator yielding all the 'Atom' objects in a residue.
        """
        for i in sorted(self.__Atoms.keys()):yield self.__Atoms[i]
    
    def get_domain_id(self):
        '''
        :returns: The assigned domain ID of the residue
        '''
        return self.__domain_id
    
    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    def set_name(self,new_name):
        """
        :param new_name: New name to be assigned
        """
        self.__name=new_name
    
    def set_parent(self,new_parent):
        """
        :param new_parent: New parent to be assigned
        """
        self.__parent=new_parent

    def set_domain_id(self,new_domain_id):
        self.__domain_id=new_domain_id

    #Calculation Functions
    def get_calpha(self):
        """
        :returns: A 'C-Alpha' atom as an 'Atom' object.
        """
        return [i for i in self.get_atoms() if i.get_name()=='CA'][0]
    
    def get_centerofgravity(self):
        """
        NOTE: Yet to add the atomic masses.
        :returns: Cartesian Coordinates of the centre of the gravity.
        """
        atoms=self.get_atoms()
        AtomicMass=1
        XYZ_M=[0,0,0]
        MassofAA=0
        for i in atoms:
            XYZ_M[0]+=i.Coordinates[0]*AtomicMass
            XYZ_M[1]+=i.Coordinates[1]*AtomicMass
            XYZ_M[2]+=i.Coordinates[2]*AtomicMass
            MassofAA=MassofAA+AtomicMass
        return numpy.array([i/MassofAA for i in XYZ_M])
    
    def get_backbone(self):
        """
        NOTE: Some files like 1k20 are showing multiple backbone atoms because seperation of models is based on [Model] in this tool
        :returns: 'Backbone' atoms as an 'Atom' object.
        """
        return [i for i in self.get_atoms() if(i.get_name()=='N' or i.get_name()=='CA' or i.get_name()=='C' or i.get_name()=='O')]
    
    def get_tip(self):
        """
        :returns: Tip atom of an amino acid as an 'Atom' object.
        """
        CAlpha=self.get_calpha()
        resname=self.get_name()
        TipofAA=None
        if(resname=='ALA' or resname=='GLY'):
            TipofAA=CAlpha
        else:
            MaxDistance=0
            for i in self.get_atoms():
                tempdistance=CAlpha.calc_dist(i)
                if(tempdistance>MaxDistance):
                    MaxDistance=tempdistance
                    TipofAA=i
        return TipofAA



class HetMol():
    def __init__(self,id,name,parent):
        self.__id=id
        self.__name=name
        self.__parent=parent
        self.__Atoms=None
    
    def __setitem__(self,id,Atom):
        try:
            self.__Atoms[id]=Atom
        except:
            self.__Atoms={}
            self.__Atoms[id]=Atom
    
    #Get Functions
    def get_id(self):
        """
        :returns: ID of an 'Residue' object
        """
        return self.__id
    
    def get_name(self):
        """
        :returns: Name of an 'Residue' object
        """
        return self.__name

    def get_parent(self):
        """
        :returns: Parent residue of an 'Residue' Object. (If available)
        """
        return self.__parent
    
    def get_hetatoms(self):
        """
        :returns: Generator yielding all the 'Atom' objects in a residue.
        """
        for i in sorted(self.__Atoms.keys()):yield self.__Atoms[i]
    
    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    def set_name(self,new_name):
        """
        :param new_name: New name to be assigned
        """
        self.__name=new_name
    
    def set_parent(self,new_parent):
        """
        :param new_parent: New parent to be assigned
        """
        self.__parent=new_parent

    def get_centerofgravity(self):
        """
        NOTE: Yet to add the atomic masses.
        :returns: Cartesian Coordinates of the centre of the gravity.
        """
        atoms=self.get_hetatoms()
        AtomicMass=1
        XYZ_M=[0,0,0]
        MassofAA=0
        for i in atoms:
            XYZ_M[0]+=i.Coordinates[0]*AtomicMass
            XYZ_M[1]+=i.Coordinates[1]*AtomicMass
            XYZ_M[2]+=i.Coordinates[2]*AtomicMass
            MassofAA=MassofAA+AtomicMass
        return numpy.array([i/MassofAA for i in XYZ_M])



class Chain():
    '''
    NOTE1: As of now, iterating over chain only fetches the residues not the hetero atoms
    NOTE2: The object fetches the atoms and residues which are not in order as they appear in the PDB. Find a way to fix this
    NOTE3: Chain has no get parent as of now | Resolved
    '''
    def __init__(self,id):
        self.__id=id
        self.__Residues=None
        self.__HetMols=None
        self.__parent=None
        #More Features
        self.__Hinges=[]
    
    def __setitem__(self,Number,Entity,Type):
        if(Type=='Residue'):
            try:
                self.__Residues[Number]=Entity
            except:
                self.__Residues={}
                self.__Residues[Number]=Entity
        elif(Type=='HetMol'):
            try:
                self.__HetMols[Number]=Entity
            except:
                self.__HetMols={}
                self.__HetMols[Number]=Entity

    def __getitem__(self,Number,Type=None):
        try:
            return self.__Residues[Number]
        except:
            return self.__HetMols[Number]
    
    #Get Functions
    def get_id(self):
        """
        :returns: ID of an 'Chain' object
        """
        return self.__id
    
    def get_residues(self):
        """
        :returns: Generator yielding all the 'Residue' objects in a Chain.
        """
        for i in sorted(self.__Residues.keys()):yield self.__Residues[i]
    
    def get_calpha(self):
        """
        :returns: Generator yielding all the 'C-Alpha' Atom objects in a Chain.
        """
        for i in self.get_residues():yield i.get_calpha()
    
    def get_hetmols(self):
        """
        :returns: Generator yielding all the 'HetMol' objects in a Chain.
        """
        for i in sorted(self.__HetMols.keys()):yield self.__HetMols[i]

    
    def get_hinges(self):
        """
        :returns: The hinge residues
        """
        return self.__Hinges
    
    def get_parent(self):
        """
        """
        return self.__parent

    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    def set_hinges(self,new_hinges):
        """
        """
        self.__Hinges=self.__Hinges+new_hinges
    
    def set_parent(self,parent):
        """
        """
        self.__parent=parent
    
    #Compute Functions
    def get_backbone(self):
        """
        :returns: All 'Backbone' atoms as an 'Atom object.'
        """
        return [i.get_backbone() for i in self.get_residues()]



class Model():
    def __init__(self,id,AllAtoms,AllResidues,AllChains,AllHetAtoms,AllHetMols):
        self.__id=id
        self.__AllAtoms=AllAtoms
        self.__AllResidues=AllResidues
        self.__AllChains=AllChains
        self.__AllHetAtoms=AllHetAtoms
        self.__AllHetMols=AllHetMols

    def __getitem__(self,ChainID):
        return self.__AllChains[ChainID]
    
    #Get Functions
    def get_id(self):
        """
        :returns: ID of an 'Model' object
        """
        return self.__id

    def get_chains(self):
        """
        :returns: Constructer for all the chains as a 'Chain' object
        """
        for i in sorted(self.__AllChains.keys()):yield self.__AllChains[i]

    def get_residues(self):
        """
        :returns: Generator yielding all the 'Residue' objects in a Model
        """
        return [j for i in self.__AllChains.keys() for j in self.__AllChains[i].get_residues()]
    
    def get_atoms(self):
        """
        :returns: Generator yielding all the 'Atom' objects in a Model
        """
        for i in sorted(self.__AllAtoms.keys()):yield self.__AllAtoms[i]
    
    def get_chain(self,ChainID):
        """
        :returns: All chains as a 'Chain' objects in a Model
        """
        return self.__AllChains[ChainID]
    
    def get_atom_byid(self,query_atom_id):
        """
        :param query_atom_id: ID of the atom to be returned
        :returns: Atom in from of 'Atom'.
        """
        return self.__AllAtoms[query_atom_id]
    
    def get_hetmols(self):
        """
        :returns: Generator yielding all the 'HetMol' objects in a Model
        """
        for i in sorted(self.__AllHetMols.keys()):yield self.__AllHetMols[i]
     
    def get_hetatoms(self):
        """
        :returns: Generator yielding all the 'HetAtom' objects in a Model
        """
        for i in sorted(self.__AllHetAtoms.keys()):yield self.__AllHetAtoms[i]
    
    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    def set_torsion(self,Atom1,Atom2,Atom3,Atom4,theta):
        """
        NOTE: Currently we can only set torsion for the backbone atoms (Psi and Phi)
        """
        radang=self.get_torsion(Atom1,Atom2,Atom3,Atom4)
        Atom1=self.get_atom_byid(Atom1)
        Atom2=self.get_atom_byid(Atom2)
        Atom3=self.get_atom_byid(Atom3)
        Atom4=self.get_atom_byid(Atom4)
        torsion_type=Atom1.get_name()+Atom2.get_name()+Atom3.get_name()+Atom4.get_name()

        if(torsion_type=='NCACN'):
            rotang=numpy.deg2rad(theta)-radang
        elif(torsion_type=='CNCAC'):
            rotang=radang-numpy.deg2rad(theta)
        else:
            raise Exception('Not a Psi/Phi angle. Current version of PACKMAN only supports Psi/Phi rotation.')
        
        sn=numpy.sin(rotang)
        cs=numpy.cos(rotang)
        t=1-cs

        v2x=Atom2.get_location()[0]-Atom3.get_location()[0]
        v2y=Atom2.get_location()[1]-Atom3.get_location()[1]
        v2z=Atom2.get_location()[2]-Atom3.get_location()[2]
        #Normalize the rotation vector
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
        tx = Atom2.get_location()[0]
        ty = Atom2.get_location()[1]
        tz = Atom2.get_location()[2]

        import copy
        for i in self.get_atoms():
            Coordinates=i.get_location()
            #Imporve this later on, this skips some atoms
            if(torsion_type=='CNCAC' and i.get_parent().get_id()<=Atom1.get_parent().get_id()):
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
                i.set_location(numpy.array(Coordinates))
            if(torsion_type=='NCACN' and i.get_parent().get_id()>=Atom4.get_parent().get_id()):
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
                i.set_location(numpy.array(Coordinates))
        return True
    
    #Compute Functions
    def get_calpha(self):
        """
        :returns: All 'C-Alpha' atoms as an 'Atom' object.
        """
        return [i.get_calpha() for i in self.get_residues()]
    
    def get_backbone(self):
        """
        :returns: All 'Backbone' atoms as an 'Atom object.'
        """
        return [i.get_backbone() for i in self.get_residues()]
    
    def get_torsion(self,Atom1,Atom2,Atom3,Atom4):
        """
        """
        Atom1=self.get_atom_byid(Atom1)
        Atom2=self.get_atom_byid(Atom2)
        Atom3=self.get_atom_byid(Atom3)
        Atom4=self.get_atom_byid(Atom4)
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
        return radang
    
    #Check Function
    def check_clashes(self,distance=0.77):
        from scipy.spatial import KDTree
        all_atoms=[i for i in self.get_atoms()]
        T=KDTree([i.get_location() for i in all_atoms])
        return len(T.query_pairs(distance))


class Protein():
    def __init__(self,id,name,Models):
        self.id=id
        self.name=name
        self.Models=Models
    
    def __getitem__(self,ModelNumber):
        return self.Models[ModelNumber]
    
    def write_pdb(self,filename):
        #sys.stdout.write("%-6s %-50s %-25s\n" % (code, name, industry))
        open(filename,'w').write('')
        fh=open(filename,'a')
        for num_,_ in enumerate(self):
            fh.write("Model\t"+str(num_)+'\n')
            for i in _.get_atoms():
                fh.write("ATOM  %5s %-4s %3s %1s%4s    %8s%8s%8s%6s%6s         %-4s%2s%2s\n"%(i.get_id(),i.get_name(),i.get_parent().get_name(),i.get_parent().get_parent().get_id(),i.get_parent().get_id(),numpy.around(i.get_location()[0],decimals=3),numpy.around(i.get_location()[1],decimals=3),numpy.around(i.get_location()[2],decimals=3),i.get_occupancy(),i.get_bfactor(),'',i.get_element(),''))
            fh.write("ENDMDL\n")
        return True

'''
##################################################################################################
#                                 Structural Super Elements                                      #
##################################################################################################
'''

class Hinge():
    '''
    Information about the hinge is stored here
    '''
    def __init__(self,hid,elements,stats,p):
        self.__id=hid
        self.__elements=elements
        self.__stats=stats
        self.__p=p

    #Get functions
    def get_id(self):
        return self.__id

    def get_elements(self):
        return self.__elements

    def get_stats(self):
        return self.__stats

    def get_pvalue(self):
        return self.__p

'''
##################################################################################################
#                                          File Load                                             #
##################################################################################################
'''

def LoadPDB(filename):
    Models=[]
    start=0

    fh=open(filename).read()

    frames=fh.split('\nMODEL')

    if(len(frames)>1):
        start=1
    else:
        start=0
    for FrameNumber,frame in enumerate(frames[start:]):
        #Map
        AllAtoms={}
        AllResidues={}
        AllChains={}

        AllHetAtoms={}
        AllHetMols={}

        lines=frame.split('\n')   
        for _ in lines:
            if(_[0:4]=='ATOM'):
                #NOTE: MAPS CAN BE REMOVED SAFELY
                #Chain Defined
                ChainID=_[21]
                if(ChainID not in AllChains.keys()):AllChains[ChainID]=Chain(ChainID)

                #Residue Defined
                ResidueNumber=int(_[22:26].strip())
                ResidueName=_[17:20]
                if(str(ResidueNumber)+ChainID not in AllResidues.keys()):AllResidues[str(ResidueNumber)+ChainID]=Residue(ResidueNumber,ResidueName,AllChains[ChainID])

                #Residue Added to the chain
                AllChains[ChainID].__setitem__(ResidueNumber,AllResidues[str(ResidueNumber)+ChainID],Type='Residue')

                #Atom Defined
                id=int(_[6:11])
                AtomName=_[12:16].strip()
                Coordinates=numpy.array([float(_[30:38]),float(_[38:46]),float(_[46:54])])
                Occupancy=float(_[54:60])
                bfactor=float(_[60:66])
                Element=_[76:78].strip()
                Charge=_[78:80]
                AllAtoms[id]=Atom(id,AtomName,Coordinates,Occupancy,bfactor,Element,Charge,AllResidues[str(ResidueNumber)+ChainID])
                
                #Atom added to the residue
                AllResidues[str(ResidueNumber)+ChainID].__setitem__(id,AllAtoms[id])

                #What to do with these?
                AlternateLocationIndicator=_[16]
                CodeForInsertions=_[26]
                SegmentIdentifier=_[72:76]
            elif(_[0:6]=='HETATM'):
                ChainID=_[21]
                if(ChainID not in AllChains.keys()):AllChains[ChainID]=Chain(ChainID)

                #HetMol Defined
                HetMolNumber=int(_[22:26].strip())
                HetMolName=_[17:20]
                if(str(HetMolNumber)+ChainID not in AllHetMols.keys()):AllHetMols[str(HetMolNumber)+ChainID]=HetMol(HetMolNumber,HetMolName,AllChains[ChainID])

                #HetMol Added to the chain
                AllChains[ChainID].__setitem__(HetMolNumber,AllHetMols[str(HetMolNumber)+ChainID],Type='HetMol')

                #HetAtom Defined
                id=int(_[6:11])
                AtomName=_[12:16].strip()
                Coordinates=numpy.array([float(_[30:38]),float(_[38:46]),float(_[46:54])])
                Occupancy=float(_[54:60])
                bfactor=float(_[60:66])
                Element=_[76:78].strip()
                Charge=_[78:80]
                AllHetAtoms[id]=HetAtom(id,AtomName,Coordinates,Occupancy,bfactor,Element,Charge,AllHetMols[str(HetMolNumber)+ChainID])
                
                #Atom added to the residue
                AllHetMols[str(HetMolNumber)+ChainID].__setitem__(id,AllHetAtoms[id])

                #What to do with these?
                AlternateLocationIndicator=_[16]
                CodeForInsertions=_[26]
                SegmentIdentifier=_[72:76]

        #print
        Models.append(Model(FrameNumber,AllAtoms,AllResidues,AllChains,AllHetAtoms,AllHetMols))
        for i in AllChains:AllChains[i].set_parent(Models[FrameNumber])

    if(len(Models)>1):
        #NMR DETECTED
        print("Caution: NMR Structure")
        warnings.warn('Multiple models are detected. B-Factor field is turned to a calculated parameter.',UserWarning)
        All_Coords=[]
        for i in Models:
            All_Coords.append(numpy.array([j.get_location() for j in i.get_atoms()]))
        All_Coords=numpy.array(All_Coords)
        
        flattened_std=[]
        for i in range(0,All_Coords.shape[1]):
            xyz_var=0
            for j in All_Coords[:,i].T:
                xyz_var=xyz_var+numpy.var(j)
            flattened_std.append(numpy.sqrt(xyz_var))
        
        for i in Models:
            for numj,j in enumerate(i.get_atoms()):
                j.set_bfactor(flattened_std[numj])
        
    return Protein(filename,None,Models)


'''
##################################################################################################
#                                   Download and Write                                           #
##################################################################################################
'''

def DownloadPDB(pdbid,save_name='Download.pdb'):
    '''
    INFO: This class is used to download a PDB stucture from RCSB PDB
    '''
    import urllib.request as ur
    response=ur.urlopen('https://files.rcsb.org/view/'+pdbid+'.pdb')
    try:
        open(save_name,'wb').write(response.read())
    except(IOError):
        None
    return True


'''
##################################################################################################
#                                         Interface                                              #
##################################################################################################
'''

def IO():
    '''
    INFO: Argument parser to the program for now.
    '''
    parser=argparse.ArgumentParser()
    parser.add_argument('filename',metavar='PDBID')
    args=parser.parse_args()
    return args


def main():
    print("NOTE: Use this mode as a")
    return True

if(__name__=='__main__'):
  main()