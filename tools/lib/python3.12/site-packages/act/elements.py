#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os, sys
from actutils     import *
from get_csv_rows import *

atomprops = {}
def get_atomprops():
    adata = act_library_filename("atomprops.csv")
    for row in get_csv_rows(adata, 7, delim=","):
        try:
            atomprops[row[0]] = { "symbol": row[1], "name": row[2], "atomicnumber": int(row[3]), "mass": float(row[4]),
                                  "charge": int(row[5]), "mult": int(row[6]) }
        except ValueError:
            print("Could not understand line {}".format(row))

def AtomNumberToAtomName(atomnumber:int) -> str:
    if len(atomprops) == 0:
        get_atomprops()
    for k in atomprops.keys():
        if atomprops[k]["atomicnumber"] == atomnumber and atomprops[k]["charge"] == 0:
            return k
    sys.exit("Invalid atomnumber %d" % atomnumber)

def AtomNameToAtomNumber(atomname:str) -> int:
    if len(atomprops) == 0:
        get_atomprops()
    for k in atomprops.keys():
        if atomname.upper() == k.upper():
            return atomprops[k]["atomicnumber"]
    sys.exit("Invalid atomname %s" % atomname)

def ElementName(atomname:str) -> str:
    if len(atomprops) == 0:
        get_atomprops()
    for k in atomprops.keys():
        if atomname.upper() == k.upper():
            return AtomNumberToAtomName(AtomNameToAtomNumber(atomname))
    sys.exit("Invalid element %s" % atomname)

def StringIsElement(atomname:str) -> bool:
    return atomname in atomprops

