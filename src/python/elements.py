#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os

numberToName = {}
nameToNumber = {}
def readElementsData():
    actdata = "ACTDATA"
    if actdata in os.environ:
        edata = ("%s/top/elements.dat" % os.environ[actdata])
        try:
            with open(edata, "r") as inf:
                for line in inf:
                    words = line.strip().split()
                    if len(words) == 3:
                        atomicnumber = int(words[2])
                        nameToNumber[words[1]] = atomicnumber
                        # There are the short and the long name in the file, only store the short one
                        if not atomicnumber in numberToName:
                            numberToName[atomicnumber] = words[1]
        except FileNotFound:
            sys.exit("Cannot find %s" + edata)
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
    
