#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os

from get_csv_rows import *

numberToName = {}
nameToNumber = {}
def readElementsData():
    actdata = "ACTDATA"
    if actdata in os.environ:
        edata = ("%s/top/elements.csv" % os.environ[actdata])
        for words in get_csv_rows(edata, 2):
            if len(words) >= 2:
                atomicnumber = int(words[1])
                nameToNumber[words[0]] = atomicnumber
                # There are the short and the long name in the file, only store the short one
                if not atomicnumber in numberToName:
                    numberToName[atomicnumber] = words[0]
    else:
        sys.exit("No environment variable %s" % actdata)

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
    
