'''
Author: Pranav Khade(pranavk@iastate.edu)
'''


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
