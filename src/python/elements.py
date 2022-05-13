#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os, sys
from actutils     import *
from get_csv_rows import *

numberToName = {}
nameToNumber = {}
def readElementsData():
    edata = act_library_filename("elements.csv")
    for words in get_csv_rows(edata, 2):
        if len(words) >= 2:
            try:
                atomicnumber = int(words[1])
            except ValueError:
                sys.exit("Incorrect words '%s|%s' in %s" % ( words[0], words[1], edata ))
            nameToNumber[words[0]] = atomicnumber
            nameToNumber[words[0].upper] = atomicnumber
            # There are the short and the long name in the file, only store the short one
            if not atomicnumber in numberToName:
                numberToName[atomicnumber] = words[0]

def AtomNumberToAtomName(atomnumber:int) -> str:
    if len(numberToName) == 0:
        readElementsData()
    if atomnumber in numberToName:
        return numberToName[atomnumber]
    sys.exit("Invalid atomnumber %d" % atomnumber)

def AtomNameToAtomNumber(atomname:str) -> int:
    if len(nameToNumber) == 0:
        readElementsData()
    if atomname in nameToNumber:
        return nameToNumber[atomname]
    sys.exit("Invalid atomname %s" % atomname)
    
def StringIsElement(atomname:str) -> bool:
    return atomname in nameToNumber
    
