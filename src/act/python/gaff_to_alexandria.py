#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os
from actutils     import *
from get_csv_rows import *

class GaffToAlexandria:
    '''Class to facilitate translation from GAFF names
    to Alexandria force fields.
    '''
    def __init__(self):
        self.alex_rename_table = {}
        self.elem_table = {}
        self.read()
        
    def read(self):
        map_file = act_library_filename("atomtype_mapping.csv")
        if os.path.exists(map_file):
            for words in get_csv_rows(map_file, 3):
                self.alex_rename_table[words[0]] = words[1]
                self.elem_table[words[0]] = words[2]

    def rename(self, atype:str):
        if atype in self.alex_rename_table:
            return self.alex_rename_table[atype]
        else:
            return atype
        
    def get_elem(self, atype: str):
        if atype in self.elem_table:
            return self.elem_table[atype]
        else:
            return atype
        
