'''
Author: Pranav Khade(pranavk@iastate.edu)
'''

from numpy import around

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
                fh.write("ATOM  %5s %-4s %3s %1s%4s    %8s%8s%8s%6s%6s         %-4s%2s%2s\n"%(i.get_id(),i.get_name(),i.get_parent().get_name(),i.get_parent().get_parent().get_id(),i.get_parent().get_id(),around(i.get_location()[0],decimals=3),around(i.get_location()[1],decimals=3),around(i.get_location()[2],decimals=3),i.get_occupancy(),i.get_bfactor(),'',i.get_element(),''))
            fh.write("ENDMDL\n")
        return True