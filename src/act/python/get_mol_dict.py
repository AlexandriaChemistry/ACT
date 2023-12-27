#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os, sys
try:
    import openbabel as ob
#    from pybel import *
except:
    pp = "PYTHONPATH"
    if pp in os.environ:
        print("OpenBabel not found. Check your %s environment variable, right now it is '%s'." % ( pp, os.environ[pp] ))
    else:
        print("OpenBabel not found and your %s environment variable is empty." % ( pp ))
    sys.exit("Cannot continue")
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

        if None != forcefield:
            ff = ob.OBForceField.FindForceField(forcefield)
            if not ff:
                print("No OpenBabel support for force field %s" % forcefield)
                return False
            if not ff.Setup(obmol):
                print("Could not setup the force field %s for %s" % ( forcefield, self.formula))
                return False
            else:
                if not ff.GetAtomTypes(obmol):
                    print("Could not get atomtypes from force field %s for %s" % ( forcefield, filename))
                    return False
            del ff

        g2a = GaffToAlexandria()
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
            if debug:
                print("atomtype %s index %d" % ( atomtype, index))
            X = atom.GetX()
            Y = atom.GetY()
            Z = atom.GetZ()
            atomic_number = atom.GetAtomicNum()
            mass = atom.GetExactMass()
            self.atoms.update({index: {}})
            self.atoms[index] = {"atomic_number": atomic_number, "obtype": obtype.GetValue(),
                                 "atomtype": atomtype, "mass": mass,
                                 "X": X, "Y": Y, "Z": Z}
        g2a = None
        # Add the bonds
        if debug:
            print("There are %d bonds" % obmol.NumBonds())
        for bond in ob.OBMolBondIter(obmol):
            bbb = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            self.bonds[bbb] = bond.GetBondOrder()
            if bond.IsAromatic() and forcefield == "alexandria":
                self.bonds[bbb] = 1.5
            if debug:
                print("bond %d-%d order %g" % ( bbb[0], bbb[1], self.bonds[bbb]))
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
        success  = self.analyse_obmol(obconversion, obmol, forcefield)
        # Help garbage collecting
        del obconversion
        obmol.Clear()
        del obmol
        return success
        
    def from_coords_elements_obc(self, elements, coords, obConversion, forcefield="alexandria"):
        if len(coords) != len(elements):
            print("Inconsistent input: There are %d coordinates but %d elements." % ( len(coords), len(elements) ))
            return False
        elif len(coords) == 0:
            print("No coordinates")
            return False
        atomnumber = len(elements)
        xyzstring = ("%d\nCoordinates\n" % atomnumber)
        for i in range(len(elements)):
            if len(coords[i]) != 3:
                print("coords[%d] has $d elements" % ( i, len(coords[i])))
                return False
            for m in range(3):
                if None == coords[i][m]:
                    print("Invalid coords {}".format(coords[i]))
                    return False
            try:
                xyzstring += ("%2s%22.3f%22.3f%22.3f\n" % ( elements[i], coords[i][0], coords[i][1], coords[i][2] ))
            except ValueError:
                print("Coordinates messed up {}".format(coords[i]))
                return False
        obmol = ob.OBMol()
        obConversion.ReadString(obmol, xyzstring)
        success      = self.analyse_obmol(obConversion, obmol, forcefield)
        # Help garbage collecting
        del obmol
        del xyzstring
        return success
        
    def from_coords_elements(self, elements, coords, forcefield="alexandria"):
        obConversion = ob.OBConversion()
        obConversion.SetInFormat("xyz")
        success = self.from_coords_elements_obc(elements, coords, obConversion, forcefield)
        # Help garbage collecting
        del obConversion
        return success
        
    def from_smiles(self, smiles:str, addH:bool, forcefield="alexandria"):
        obConversion = ob.OBConversion()
        obConversion.SetInFormat("smi")
        obmol = ob.OBMol()
        obConversion.ReadString(obmol, smiles)
        if addH:
            obmol.AddHydrogens()
        gen3d = ob.OBOp.FindType("Gen3D")
        gen3d.Do(obmol, "--best")
        success = self.analyse_obmol(obConversion, obmol, forcefield)
        # Help garbage collecting
        del obmol
        del obConversion
        return success
