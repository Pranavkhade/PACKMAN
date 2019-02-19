'''
Author: Pranav Khade(pranavk@iastate.edu)
SourceofInformation: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
'''
import numpy
import argparse
import os
import warnings

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
    
    def get_atoms(self):
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


class Chain():
    def __init__(self,id):
        self.__id=id
        self.__Residues=None
    
    def __setitem__(self,ResidueNumber,Residue):
        try:
            self.__Residues[ResidueNumber]=Residue
        except:
            self.__Residues={}
            self.__Residues[ResidueNumber]=Residue

    def __getitem__(self,ResidueNumber):
        return self.__Residues[ResidueNumber]
    
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
    
    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
    #Compute Functions
    def get_backbone(self):
        """
        :returns: All 'Backbone' atoms as an 'Atom object.'
        """
        return [i.get_backbone() for i in self.get_residues()]

    

class Model():
    def __init__(self,id,AllAtoms,AllResidues,AllChains):
        self.__id=id
        self.__AllAtoms=AllAtoms
        self.__AllResidues=AllResidues
        self.__AllChains=AllChains

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
        for i in sorted(self.__AllResidues.keys()):yield self.__AllResidues[i]
    
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
    
    #Set Functions
    def set_id(self,new_id):
        """
        :param new_id: New ID to be assigned
        """
        self.__id=new_id
    
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

class Protein():
    def __init__(self,id,name,Models):
        self.id=id
        self.name=name
        self.Models=Models
    
    def __getitem__(self,ModelNumber):
        return self.Models[ModelNumber]

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
                AllChains[ChainID].__setitem__(ResidueNumber,AllResidues[str(ResidueNumber)+ChainID])

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
        #print
        Models.append(Model(FrameNumber,AllAtoms,AllResidues,AllChains))

    if(len(Models)>1):
        #NMR DETECTED
        print "Caution: NMR Structure"
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


#####-----#####-----#####-----#####-----#####-----

def DownloadPDB(pdbid,save_name='Download.pdb'):
    '''
    INFO: This class is used to download a PDB stucture from RCSB PDB
    '''
    import urllib
    response=urllib.urlopen('https://files.rcsb.org/view/'+pdbid+'.pdb')
    try:
        open(save_name,'w').write(response.read())
    except(IOError):
        None
    return True

def IO():
    '''
    INFO: Argument parser to the program for now.
    '''
    parser=argparse.ArgumentParser()
    parser.add_argument('filename',metavar='PDBID')
    args=parser.parse_args()
    return args


def main():
    print "NOTE: Use this mode as a"
    return True

if(__name__=='__main__'):
  main()