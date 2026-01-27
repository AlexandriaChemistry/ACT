#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os, sys, tempfile
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from gaff_to_alexandria import *

debug = False

def get_order(order:str)->float:
    orders = { "SINGLE": 1, "DOUBLE": 2, "TRIPLE": 3, "AROMATIC": 1.5 }
    if order in orders:
        return orders[order]
    sys.exit("Unknown bond order '%s'" % order)

def get_atype(hybridization:str)->str:
    hybrid = { "SP": 1, "SP2": 2, "SP2D": 2, "SP3": 3, "SP3D": 3, "SP3D2": 3 }
    if hybridization in hybrid:
        return str(hybrid[hybridization])
    sys.exit("Unknown hybridization '%s'" % hybridization)

class MoleculeDict:
    '''Class to hold and extract molecule information'''
    def __init__(self, verbose=False):
        self.molecule = {}
        self.atoms    = []
        self.bonds    = {}
        self.verbose  = verbose
        self.charge   = 0

    def check_bondorder(self):
        # TODO: do we need to implement more stuff to use the atom_bond.xml info?
        my_pairs = { "no": "on", "cm": "om", "p5": "om", "s6": "om",
                     "py": "om", "s3": "om", "cz": "n", "s4": "om"  }
        for mp in my_pairs.keys():
            for i in range(len(self.atoms)):
                if mp == self.atoms[i]["obtype"]:
                    bondIndex = []
                    for (ai,aj) in self.bonds.keys():
                        if ((i == ai and my_pairs[mp] == self.atoms[aj]["obtype"]) or
                            (i == aj and my_pairs[mp] == self.atoms[ai]["obtype"])):
                            bondIndex.append((ai, aj))
                    if len(bondIndex) >= 2:
                        for b in bondIndex:
                            self.bonds[b] = 1.5

    def analyse(self, mol, molname:str, charge=None)->bool:
        self.title      = molname
        self.mol_weight = 0
        self.numb_atoms = len(mol.GetAtoms())
        self.formula    = Chem.rdMolDescriptors.CalcMolFormula(mol)
        if None != charge:
            self.charge     = charge
        if self.charge != Chem.GetFormalCharge(mol):
            print(f"Warning: user passed charge {charge} for {molname} but RDKit says charge is {Chem.GetFormalCharge(mol)}")
        self.inchi      = Chem.rdinchi.MolToInchi(mol)
        self.mult       = 1
        coords          = mol.GetConformer().GetPositions()
        # First do the atoms
        ii = 0
        mw = 0
        for atom in mol.GetAtoms():
            if self.verbose:
                print("Atom %d atomic number %d valence %s hybridization %s" %
                      ( ii, atom.GetAtomicNum(),
                        atom.GetValence(Chem.ValenceType.EXPLICIT),
                        atom.GetHybridization() ) )
            fftype = atom.GetSymbol()
            if atom.GetFormalCharge() == 0:
                fftype = fftype.lower()
                if fftype in [ "c", "n", "o", "p", "s" ]:
                    fftype += get_atype(str(atom.GetHybridization()))
            else:
                fftype += str(atom.GetFormalCharge())
            self.atoms.append( { "atomic_number": atom.GetAtomicNum(), "obtype": fftype,
                                 "atomtype": fftype, "mass": atom.GetMass(), "element": atom.GetSymbol(),
                                 "X": coords[ii][0], "Y": coords[ii][1], "Z": coords[ii][2] } )
            mw += atom.GetMass()
            ii += 1
        if not self.numb_atoms == ii:
            print(f"self.numb_atoms {self.numb_atoms} != len(self.atoms) {ii}")
            return False

        self.mol_weight = mw

        # Now do the bonds
        for bb in mol.GetBonds():
            ai    = bb.GetBeginAtomIdx()
            aj    = bb.GetEndAtomIdx()
            order = get_order(str(bb.GetBondType()))
            if self.verbose:
                print("Bond %d from %d to %d order %s" % ( ii, ai, aj, order ) )
            self.bonds[(ai, aj)] = order
        self.check_bondorder()

        return True

    def read(self, filename, fileformat=None, charge=None)->bool:
        success = False
        if self.verbose:
            print("Analyzing %s" % filename)
        m = None
        if None == fileformat:
            fileformat = filename[-3:]
        if fileformat == "sdf":
            m = Chem.MolFromMolFile(filename, sanitize=False, removeHs=False)
        elif fileformat == "xyz":
            raw_mol = Chem.MolFromXYZFile(filename)
            m = Chem.Mol(raw_mol)
            Chem.DetermineBonds(m)
        elif fileformat == "pdb":
            m = Chem.MolFromPDBFile(filename, removeHs=False)
        else:
            print("Don't know how to handle %s" % filename)
        if m:
            # TODO remove extension in the correct way
            molname = os.path.basename(os.path.normpath(os.path.splitext(filename)[0]))
            success = self.analyse(m, molname, charge)

        return success
        
    def from_coords_elements(self, elements, coords, charge=None):
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
        success = False
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'something.xyz')
            # use path
            with open(path, "w") as outf:
                outf.write(xyzstring)
            success = self.read(path, "xyz", charge)
            os.unlink(path)
            os.rmdir(tmp)

        return success
        
    def from_smiles(self, smiles:str, addH:bool, charge=None):
        if self.verbose:
            print("Analyzing %s" % smiles)
        m  = Chem.MolFromSmiles(smiles)
        if addH:
            m2 = Chem.AddHs(m)
            return analyse(m2, smiles, charge)
        else:
            return analyse(m, smiles, charge)
