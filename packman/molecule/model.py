# -*- coding: utf-8 -*-
"""The 'Model' object host file.

This is file information, not the class information. This information is only for the API developers.
Please read the 'Model' object documentation for details.

Example::

    from packman.molecule import Model
    help( Model )

Note:
    * The models are nothing but frames of the PDB file.
    
Todo:
    * Finish writing up the documentation.
    * Finish error handling.
    * Finish optimizing the performance.

Authors:
    * Pranav Khade(https://github.com/Pranavkhade)
"""

from ..entropy import PackingEntropy

import numpy
import logging


class Model():
    """This class contains the information about the 'Chain' object (packman.molecule.Chain).

        This class contains all the information available about the Chain and stores the corresponding 'Residue' and 'Hetmol' objects in itself. The Chain class is the third lowest in the hierarchy of the 'molecule' API classes.
        the order of hierarchy being: Protein> Model> Chain> Residue> Atom. This class is also the component of the 'molecule' module API.
        Please read the Tutorials and Documentation for more details.

        Notes:
            * Please refer to the [https://web.archive.org/web/20080905024351/http://www.wwpdb.org/docs.html] for the description of the arguments.
        
        Args:
            id (int)                                : Model ID from the PDB file ordered from first to the last. Each Model in a PDB file has unique ID. (essential)
            AllAtoms ({packman.molecule.Atom})      : Dictionary of all the 'Atom' in the given model.
            AllResidues ({packman.molecule.Residue}): Dictionary of all the 'Residue' in the given model.
            AllChains ({packman.molecule.Chain})    : Dictionary of all the 'Chain' in the given model.
            AllHetAtoms ({packman.molecule.HetAtom}): Dictionary of all the 'HetAtom' in the given model.
            AllHetMols ({packman.molecule.HetMol})  : Dictionary of all the 'HetMol' in the given model.

        """
        
    def __init__(self,id,AllAtoms,AllResidues,AllChains,AllHetAtoms,AllHetMols):                
        self.__id=id
        self.__AllAtoms=AllAtoms
        self.__AllResidues=AllResidues
        self.__AllChains=AllChains
        self.__AllHetAtoms=AllHetAtoms
        self.__AllHetMols=AllHetMols
        self.__parent = None
        
        #Properties are the entities that are not included in the PDB files and are obtained by calculations
        self.__properties = {}

    def __getitem__(self,ChainID):
        return self.__AllChains[ChainID]
    
    #Get Functions
    def get_id(self):
        """Get the ID of the 'Model'

        Returns:
            int if successful, None otherwise.
        """
        return self.__id

    def get_chains(self):
        """Get the list of corresponding 'Chain' objects of the 'Model'

        Returns:
            [packman.molecule.Chain] if successful, None otherwise.
        """
        for i in sorted(self.__AllChains.keys()):yield self.__AllChains[i]

    def get_residues(self):
        """Get the generator of corresponding 'Residue' objects of the 'Model'

        Returns:
            array of 'Residue' objects if successful, None otherwise.
        """
        #return [j for i in self.__AllChains.keys() for j in self.__AllChains[i].get_residues()]
        residues = []
        for i in self.__AllChains.keys():
            try:
                residues.extend( self.__AllChains[i].get_residues() )
            except:
                logging.warning("Chain "+str(i)+" either doesn't have residues or an error occurred; Model.get_residues() may have loaded other chains.")
        return residues
    
    def get_atoms(self):
        """Get the generator of corresponding 'Atom' objects of the 'Model'

        Returns:
            generator of 'Atom' objects if successful, None otherwise.
        """
        for i in sorted(self.__AllAtoms.keys()):yield self.__AllAtoms[i]
    
    def get_chain(self,ChainID):
        """Get the corresponding 'Chain' object

        Returns:
            'Chain' object if successful, None otherwise.
        """
        return self.__AllChains[ChainID]
    
    def get_atom_byid(self,query_atom_id):
        """Get the 'Atom' with corresponding 'Atom' ID

        Returns:
            'packman.molecule.Atom' object if successful, None otherwise.
        """
        return self.__AllAtoms[query_atom_id]
    
    def get_hetmols(self):
        """Get the generator of corresponding 'HetMol' objects of the 'Model'

        Returns:
            generator of 'packman.molecule.HetMol' objects if successful, None otherwise.
        """
        for i in sorted(self.__AllHetMols.keys()):yield self.__AllHetMols[i]
     
    def get_hetatoms(self):
        """Get the generator of corresponding 'HetAtom' objects of the 'Model'

        Returns:
            generator of 'packman.molecule.HetAtom' objects if successful, None otherwise.
        """
        for i in sorted(self.__AllHetAtoms.keys()):yield self.__AllHetAtoms[i]
    
    def get_parent(self):
        """Get the 'Protein' parent of the 'Model' object.

        Returns:
            'packman.molecule.Protein' object if successful, None otherwise.
        """
        return self.__parent
    
    def get_entropy(self,entropy_type):
        """Get the Packing Entropy of the given 'Chain'.

        Args:
            type (str): Type of entropy (Allowed Values: 1. PackingEntropy)

        Note:
            - More type of Entropies might be added in the future.
        """
        EntropyTypes = ['PackingEntropy']
        try:
            return numpy.sum( [self.__AllResidues[i].get_entropy(entropy_type) for i in self.__AllResidues] )
        except:
            if(entropy_type in EntropyTypes):
                logging.warning('This Entropy type might not be calculated for all the chains/residues. Please check the "calculate_entropy" function in the documentation for the details.')
            else:
                logging.warning('The Entropy type provided is invalid. Please check the documentation for the details.')

    def get_property(self,property_name):
        """Get the Property of the given 'Chain'.

        Property is any key and value combination that can be assigned to this object. This (along with the set_property) feature is mainly useful for the user customization.
        Properties are like pinboards. You can pin anything to the object with a key as a pin.

        Args:
            property_name (object): The 'Key' or a name the user wants to assign to to the property
        
        Note:
            - Users can add custom annotations; for example: If particular chain becomes disordered, it can be annotated with this feature.
        """
        try:
            return self.__properties[property_name]
        except:
            logging.warning('The Property Name provided is not assigned.')


    #Compute Functions
    def get_calpha(self):
        """Get the C-Alpha atom of the 'Model' as an 'Atom' object.

        Returns:
            list of packman.molecule.Atom if successful, None otherwise.
        """
        return [i.get_calpha() for i in self.get_residues()]
    
    def get_backbone(self):
        """Get the Backbone atoms of the given 'Model' as a list of 'Atom' object

        Note:
            * Backbone Atoms: CA, O, N, C

        Returns:
            list of packman.molecule.Atom if successful, None otherwise.
        """
        return [i.get_backbone() for i in self.get_residues()]
    
    def get_torsion(self,Atom1,Atom2,Atom3,Atom4):
        """Get the torsion angle between the two planes defined by four atoms.

        Notes:
            * This function will be moved to the molecule manipulation package later

        
        Args:
            Atom1 (packman.molecule.Atom): First reference 'Atom'
            Atom2 (packman.molecule.Atom): Second reference 'Atom'
            Atom3 (packman.molecule.Atom): Third reference 'Atom'
            Atom4 (packman.molecule.Atom): Fourth reference 'Atom'
        
        Returns:
            theta (float)                : Angle user wishes to set (radian)
        
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
    
        #Set Functions
    def set_id(self,new_id):
        """Set the ID of the given 'Model'

        Args:
            new_id (int): The ID User wishes to assign to the given 'Model'
        """
        self.__id=new_id
    
    def set_parent(self, new_parent):
        """Set the 'Protein' object as a parent to the 'Model' object.

        Args:
            new_parent (packman.molecule.Protein): The 'Protein' object as a parent to the given 'Model'
        """
        self.__parent = new_parent
    
    def set_torsion(self,Atom1,Atom2,Atom3,Atom4,theta):
        """Set the torsion angle between the two planes defined by four atoms.

        Notes:
            * This function will be moved to the molecule manipulation package later

        
        Args:
            Atom1 (packman.molecule.Atom): First reference 'Atom'
            Atom2 (packman.molecule.Atom): Second reference 'Atom'
            Atom3 (packman.molecule.Atom): Third reference 'Atom'
            Atom4 (packman.molecule.Atom): Fourth reference 'Atom'
            theta (float)                : Angle user wishes to set (Degree)
        
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

    def set_property(self,property_name,value):
        """Set the Property of the given 'Model'.

        Property is any key and value combination that can be assigned to this object. This (along with the get_property) feature is mainly useful for the user customization.
        Properties are like pinboards. You can pin anything to the object with a key as a pin.
        
        Args:
            property_name (object): The 'Key' or a name the user wants to assign to to the property
            value (object):         The value the user wants to assign to the property
        
        Note:
            - Users can add custom annotations; for example: If particular amino acid becomes disordered, it can be annotated with this feature.
        """
        try:
            self.__properties[property_name] = value
        except:
            logging.warning('Please check the property name. Check the allowed Python dictionary key types for more details.')

    #Calculate Functions
    def calculate_entropy(self,entropy_type,chains=None, probe_size=1.4, onspherepoints=30):
        """Calculate the entropy for the each amino acid will be returned.
    
        The 'chains' argument should be used when the user wants to restrict the analysis to a chain or group of chains rather than the whole structure.

        Args:
            entropy_type (str)              : Type of entropy to be calculated (Options: PackingEntropy)
            chains ([str]/str)              : Chain IDs for the Entropy calculation (None means all the chains are included; single string means only one chain ID; multiple chains should be an array of strings).
            probe_size (float)              : Radius of the probe to generate the surface points (This value should not be less than 1;Read the Publication for more details)
            onspherepoints (int)            : Number of points to be generated around each point for the surface (Read the Publication for more details)
        """
        if(entropy_type=='PackingEntropy'):
            PackingEntropy(self.get_atoms(), chains=chains, probe_size=probe_size, onspherepoints=onspherepoints)
        else:
            logging.warning("Please provide the valid type for the Entropy Calculation.")

    #Check Function
    def check_clashes(self,distance=0.77):
        """Check if any atoms are too close to each other. This is important since too close atoms in the elastic network models can be very bad for the results.

        Notes:
            * This function will be moved to the molecule manipulation package later

        
        Args:
            distance (float): The distance cutoff user wishes to defined as a clash radius (default:0.77; max bond length)
        
        Returns:
            clashes (int)   : Number of clashes present according to the set cutoff.
        """
        from scipy.spatial import KDTree
        all_atoms=[i for i in self.get_atoms()]
        T=KDTree([i.get_location() for i in all_atoms])
        return len(T.query_pairs(distance))
