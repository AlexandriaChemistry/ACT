1#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import math, os, sys
from actutils     import *
from get_csv_rows import *

debug = False

def compute_dhform(energyHF:float, atomtypes, g2a, ahof,
                   leveloftheory:str, temperature:float) -> float:
    eatom = 0
    for a in atomtypes:
        myelem = g2a.get_elem(a)
        ae     = ahof.get_atomization(myelem, leveloftheory, temperature)
        if None == ae:
            sys.exit("Cannot find atomization energy for %s with %s at %f K" % ( myelem, leveloftheory, temperature))
        eatom += ae
    return energyHF - eatom
        
def UnitToConversionFactor(unit:str):
    KCAL = 4.184
    if unit == "eV" or unit == "electronvolt":
        return 23.060538*KCAL
    elif unit == "kJ/mol" or unit == "J/mol K":
        return 1.0
    elif unit == "kcal/mol":
        return KCAL
    elif unit == "Hartree":
        return 627.509469*KCAL
    else:
        sys.exit("Unknown unit %s" % unit)

class AtomicHOF:
    def __init__(self, method, temperature, verbose=False):
        self.method      = method
        self.temperature = temperature
        self.verbose     = verbose
        self.read()
        
    def read(self):
        for datafile in [ act_library_filename("atomization-energies.csv"),
                          act_library_filename("atomization-energies-dft.csv") ]:
            if not os.path.exists(datafile):
                sys.exit("Cannot find %s, please reinstall ACT" % datafile)
            rows = get_csv_rows(datafile, 8)
            self.ahof = {}
            for row in rows:
                try:
                    mytemp = float(row[4])
                    if ((row[2] == "exp" or row[2] == self.method) and 
                        (mytemp == 0 or abs(mytemp-self.temperature) <= 1e-2)):
                        akey = row[0]+"|"+row[1]
                        if not akey in self.ahof:
                            self.ahof[akey] = []
                        self.ahof[akey].append({ "Method": row[2],
                                                 "Desc": row[3],
                                                 "Temp": float(row[4]),
                                                 "Value": float(row[5]),
                                                 "Multiplicity": row[6] , 
                                                 "Unit": row[7]})
                except ValueError:
                    print("It seems that '%s' is not a number" % row[4])

    def get_atomization(self, elem, method, temp):
        akey = elem + "|0"
        if not akey in self.ahof:
            sys.exit("Cannot find key %s in the Atomic Heat of Formation table. Method is %s" % ( akey, self.method ))
        for p in range(len(self.ahof[akey])):
            thisprop = self.ahof[akey][p]
            if thisprop["Method"] == method and thisprop["Temp"] == temp:
                return thisprop["Value"]
        return None
        
    def get(self, elem, temp, charge):
        akey = elem + "|" + str(charge)
        if not akey in self.ahof:
            sys.exit("Cannot find key %s in the Atomic Heat of Formation table. Method is %s" % ( akey, self.method ))
        HexpT  = None
        S0     = None
        Vdhf   = None
        Vmodel = None
        for p in range(len(self.ahof[akey])):
            thisprop = self.ahof[akey][p]
            if self.verbose and debug:
                print("prop %d method %s desc %s temp %g" % 
                      ( p, thisprop["Method"], thisprop["Desc"], thisprop["Temp"]))
            eFac = UnitToConversionFactor(thisprop["Unit"])
            if abs(temp - thisprop["Temp"]) <= 1e-2:
                if thisprop["Method"] == "exp":
                    if thisprop["Desc"] == "H(0)-H(T)":
                        HexpT = thisprop["Value"]
                        if self.verbose:
                            print("found HexpT = %g" % HexpT)
                    elif thisprop["Desc"] == "S0(T)":
                        S0 = thisprop["Value"]
                        if self.verbose:
                            print("found S0 = %g" % S0)
            elif thisprop["Temp"] == 0:
                if thisprop["Method"] == "exp":
                    if thisprop["Desc"] == "DHf(T)":
                        Vdhf = thisprop["Value"]
                        if self.verbose:
                            print("found Vdhf = %g" % Vdhf)
                else:
                    Vmodel = thisprop["Value"]*eFac
                    if self.verbose:
                        print("found Vmodel = %g" % Vmodel)
        # Check whether we found everything
        if HexpT and S0 and Vdhf and Vmodel:
            dhof0 = Vdhf-Vmodel
            return dhof0, (dhof0-HexpT), -S0
        else:
            return None, None, None
