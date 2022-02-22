#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os

class GaffToAlexandria:
    '''Class to facilitate translation from GAFF names
    to Alexandria force fields.
    '''
    def __init__(self):
        self.alex_rename_table = {}
        self.read()
        
    def read(self):
        actdata      = "ACTDATA"
        if actdata in os.environ:
            map_file = os.environ[actdata] + "/top/atomtype_mapping.dat"
            if os.path.exists(map_file):
                with open(map_file, "r") as inf:
                    for line in inf:
                        words = line.strip().split()
                        if len(words) == 2:
                            self.alex_rename_table[words[0]] = words[1]

    def rename(self, atype:str):
        if atype in self.alex_rename_table:
            return self.alex_rename_table[atype]
        else:
            return atype
        
