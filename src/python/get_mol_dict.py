#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
try:
    import openbabel as ob
except:
    print("OpenBabel not found. Check your PYTHONPATH environment variable.")
    print("Proceed at your own risk")
import os, sys
from gaff_to_alexandria import *

debug = False

class MoleculeDict:
    '''Class to hold and extract molecule information'''
    def __init__(self):
        self.molecule = {}
        self.atoms    = {}
        self.bonds    = {}
        
    def check_bondorder(self):
        my_pairs = { "no": "on", "cm": "om", "p5": "om", "s6": "om",
                     "py": "om", "s3": "om", "cz": "n", "s4": "om"  }
        for mp in my_pairs.keys():
            for i in range(1, 1+len(self.atoms)):
                if mp == self.atoms[i]["obtype"]:
                    bondIndex = []
                    for (ai,aj) in self.bonds.keys():
                        if ((i == ai and my_pairs[mp] == self.atoms[aj]["obtype"]) or
                            (i == aj and my_pairs[mp] == self.atoms[ai]["obtype"])):
                            bondIndex.append((ai, aj))
                    if len(bondIndex) >= 2:
                        for b in bondIndex:
                            self.bonds[b] = 1.5

    def analyse_obmol(self, obconversion, obmol, forcefield=None):
        self.title      = obmol.GetTitle()
        self.mol_weight = obmol.GetMolWt()
        self.numb_atoms = obmol.NumAtoms()
        self.formula    = obmol.GetFormula()
        self.charge     = obmol.GetTotalCharge()
        obmol.AssignTotalChargeToAtoms(self.charge)
        obmol.SetAromaticPerceived(False) 
        obconversion.SetOutFormat("inchi")
        self.inchi      = obconversion.WriteString(obmol).strip()
        obconversion.SetOutFormat("inchikey")
        self.inchikey   = obconversion.WriteString(obmol).strip()

        if None != forcefield:
            ff = ob.OBForceField.FindForceField(forcefield)
            if not ff:
                print("No OpenBabel support for force field %s" % forcefield)
                return False
            if not ff.Setup(obmol):
                print("Could not setup the force field %s for %s" % ( forcefield, filename))
                return False
            else:
                if not ff.GetAtomTypes(obmol):
                    print("Could not get atomtypes from force field %s for %s" % ( forcefield, filename))
                    return False

        g2a      = GaffToAlexandria()
        # Add the atoms
        for atom in ob.OBMolAtomIter(obmol):
            index      = atom.GetIdx()
            atomtype   = "Z"
            ffatomtype = "FFAtomType"
            obtype     = ""
            if forcefield and atom.HasData(ffatomtype):
                atp    = atom.GetData(ffatomtype)
                obtype = atp
                if atp:
                    atomtype = g2a.rename(atp.GetValue())
            X = atom.GetX()
            Y = atom.GetY()
            Z = atom.GetZ()
            atomic_number = atom.GetAtomicNum()
            mass = atom.GetExactMass()
            self.atoms.update({index: {}})
            self.atoms[index] = {"atomic_number": atomic_number, "obtype": obtype.GetValue(),
                                 "atomtype": atomtype, "mass": mass,
                                 "X": X, "Y": Y, "Z": Z}
    
        # Add the bonds
        for bond in ob.OBMolBondIter(obmol):
            bbb = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            self.bonds[bbb] = bond.GetBondOrder()
            if bond.IsAromatic() and forcefield == "alexandria":
                self.bonds[bbb] = 1.5
        if forcefield == "alexandria":
            self.check_bondorder()
        if debug:
            print(self)
            
        return True

    def read(self, filename, fileformat, forcefield="alexandria"):
        obconversion = ob.OBConversion()
        obconversion.SetInFormat(fileformat)
        obmol    = ob.OBMol()
        notatend = obconversion.ReadFile(obmol,filename)
        return self.analyse_obmol(obconversion, obmol, forcefield)
        
    def from_coords_elements(self, elements, coords, forcefield="alexandria"):
        if len(coords) != len(elements):
            print("Inconsistent input: There are %d coordinates but %d elements." % ( len(coords), len(elements) ))
            return False
        elif len(coords) == 0:
            print("No coordinates")
            return False
        atomnumber = len(elements)
        xyzstring = ("%d\nCoordinates\n" % atomnumber)
        for i in range(len(elements)):
            xyzstring += ("%2s%22.3f%22.3f%22.3f\n" % ( elements[i], coords[i][0], coords[i][1], coords[i][2] ))

        obConversion = ob.OBConversion()
        obConversion.SetInFormat("xyz")
        obmol = ob.OBMol()

        obConversion.ReadString(obmol, xyzstring)
        return self.analyse_obmol(obConversion, obmol, forcefield)
