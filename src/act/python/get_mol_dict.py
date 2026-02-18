#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os, sys, tempfile, xmltodict
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import AllChem

debug     = False
atomtype  = "atomtype"
atomtypes = "atomtypes"
bondtype  = "bondtype"
bondtypes = "bondtypes"

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

def get_atom_bond_xml()->list:
    actdata = "ACTDATA"
    xmlname = "atom_bond.xml"
    dbname  = None
    if actdata in os.environ:
        dbname = ("%s/%s" % ( os.environ[actdata], xmlname ))
    if not dbname or not os.path.exists(dbname):
        dbname = xmlname
    if not os.path.exists(dbname):
        sys.exit("Cannot not find %s" % xmlname)

    with open(dbname) as fd:
        doc = xmltodict.parse(fd.read())
        abe = []
        for abtype in doc["atombondtypes"]['atombondtype']:
            ab = { "name": abtype["@name"],
                   "smarts": abtype["@smarts"],
                   "charge": abtype["@charge"],
                   "multiplicity": abtype["@multiplicity"],
                   atomtypes: [],
                   bondtypes: [] }
            atypes = abtype[atomtypes]
            if isinstance(atypes[atomtype], list):
                for atp in atypes[atomtype]:
                    if debug:
                        print(f"{ab['name']} {atp}")
                    ab[atomtypes].append({ "index": atp["@index"],
                                           "name": atp["@name"],
                                           "atomnumber": atp["@atomnumber"] })
            else:
                myatp = atypes[atomtype]
                ab[atomtypes].append({ "index": myatp["@index"],
                                       "name": myatp["@name"],
                                       "atomnumber": myatp["@atomnumber"] })
            if bondtypes in abtype and bondtype in abtype[bondtypes]:
                if isinstance(abtype[bondtypes][bondtype], list):
                    for bt in abtype[bondtypes][bondtype]:
                        ab[bondtypes].append({ "ai": bt["@ai"],
                                               "aj": bt["@aj"],
                                               "order": bt["@order"] })
                else:
                    mybtp = abtype[bondtypes][bondtype]
                    ab[bondtypes].append({ "ai": mybtp["@ai"],
                                           "aj": mybtp["@aj"],
                                           "order": mybtp["@order"] })

            abe.append(ab)
    return abe

class MoleculeDict:
    '''Class to hold and extract molecule information'''
    def __init__(self, verbose=False):
        self.molecule = {}
        self.atoms    = []
        self.bonds    = {}
        self.verbose  = verbose
        self.charge   = 0

    def lookUpSpecial(self, mol2, oneH:bool):
        abe = get_atom_bond_xml()

        NOTSET = -666
        # Check whether we have any of those special cases.
        for i in range(len(abe)):
            # Make RDKit matcher for our pattern
            pattern = Chem.MolFromSmarts(abe[i]["smarts"])
            if mol2.HasSubstructMatch(pattern):
                # Loop over all matches, there may be more than one, e.g.
                # a compound with two carboxylic groups.
                # https://www.rdkit.org/docs/GettingStartedInPython.html
                for matchPat in mol2.GetSubstructMatches(pattern):
                    mapAtoms = [NOTSET]*len(self.atoms)
                    if len(abe[i][atomtypes]) == 1:
                        mapAtoms[matchPat[0]] = 0
                    else:
                        if self.verbose:
                            print(f"{matchPat}")
                        for k in range(len(matchPat)):
                            mapAtoms[matchPat[k]] = k

                    # Update atom types
                    for j in range(len(self.atoms)):
                        if not mapAtoms[j] == NOTSET:
                            atype = abe[i][atomtypes][mapAtoms[j]]
                            isH   = atype["atomnumber"] == 1

                            if not (oneH and isH):
                                self.atoms[j]["obtype"] = atype["name"]

                    # Now fix the bonds
                    for b in self.bonds:
                        ai = mapAtoms[b[0]-1]
                        aj = mapAtoms[b[1]-1]
                        if not ai == NOTSET and not aj == NOTSET:
                            for ab in abe[i][bondtypes]:
                                if ((int(ab["ai"]) == ai and int(ab["aj"]) == aj) or
                                    (int(ab["ai"]) == aj and int(ab["aj"]) == ai)):
                                    self.bonds[b] = ab["order"]
                                    break

    def analyse(self, mol, molname:str, mycharge=None)->bool:
        self.title      = molname
        self.mol_weight = 0
        self.numb_atoms = len(mol.GetAtoms())
        try:
            self.formula    = Chem.rdMolDescriptors.CalcMolFormula(mol)
        except RuntimeError:
            print(f"Weird molecule {molname}")
            self.formula = "N/A"
        if None != mycharge:
            self.charge  = mycharge
            rdkit_charge = Chem.GetFormalCharge(mol)
            if self.charge != rdkit_charge:
                print(f"Warning: user passed charge {self.charge} for {molname} ({self.formula}) but RDKit says charge is {rdkit_charge}")
        # Returns a tuple and we need the first element
        myinchi         = Chem.rdinchi.MolToInchi(mol)
        if len(myinchi) > 0:
            self.inchi = myinchi[0]
        else:
            self.inchi = "N/A"
        self.mult       = 1
        coords          = mol.GetConformer().GetPositions()
        # First do the atoms
        ii = 0
        mw = 0
        for atom in mol.GetAtoms():
            if self.verbose:
                print("Mol %s Atom %d atomic number %d valence %s hybridization %s" %
                      ( molname, ii, atom.GetAtomicNum(),
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
            print(f"self.numb_atoms {self.numb_atoms} != len(self.atoms) {ii} for {molname}")
            return False

        self.mol_weight = mw

        # Now do the bonds
        for bb in mol.GetBonds():
            ai    = bb.GetBeginAtomIdx()
            aj    = bb.GetEndAtomIdx()
            order = get_order(str(bb.GetBondType()))
            if self.verbose:
                print("Bond %d from %d to %d order %s" % ( ii, ai, aj, order ) )
            self.bonds[(ai+1, aj+1)] = order

        # Finally, look for special cases
        self.lookUpSpecial(mol, False)

        return True

    def read(self, molname:str, filename:str,
             fileformat:str, mycharge:int)->bool:
        success = False
        if self.verbose:
            print("Analyzing %s for %s" % ( filename, molname ))
        m = None
        if None == fileformat:
            fileformat = filename[-3:]
        # First try and read, but do not crash
        try:
            if fileformat == "sdf":
                m = Chem.MolFromMolFile(filename, sanitize=True, removeHs=False)
            elif fileformat == "xyz":
                raw_mol = Chem.MolFromXYZFile(filename)
                m = Chem.Mol(raw_mol)
                if debug:
                    print("raw_mol #atoms %d mol #atoms %d" % ( len(raw_mol.GetAtoms()), len(m.GetAtoms())))
            elif fileformat == "pdb":
                try:
                    m = Chem.MolFromPDBFile(filename, sanitize=True, removeHs=False)
                except Chem.AtomValenceException:
                    m = Chem.MolFromPDBFile(filename, sanitize=False, removeHs=False)
        except ValueError:
            print(f"Problem reading {molname} from {filename}")
            return False
        if m == None:
            print(f"RDKit returned an empty molecule {molname}")
            return False
        # We trust the sdf files ...
        determine_bonds = fileformat != "sdf"
        # ... and if the user specified bonds in the pdb, weuse those.
        if fileformat == "pdb" and m.GetNumBonds() > 0:
            determine_bonds = False

        if determine_bonds:
            # The routine to DetermineBonds may crash for weird molecules
            # therefore it is good to try and catch exceptions.
            try:
                # Single ions are a special case
                matoms = m.GetAtoms()
                if len(matoms) == 1:
                    if mycharge != 0:
                        matoms[0].SetFormalCharge(mycharge)
                else:
                    # Ion pairs are a special case as well
                    if not (len(matoms) == 2 and
                            matoms[0].GetFormalCharge() != 0 and
                            matoms[1].GetFormalCharge() != 0):
                        rdDetermineBonds.DetermineBonds(m, useHueckel=False, charge=mycharge)
                try:
                    m.UpdatePropertyCache(strict=True)
                except Chem.AtomValenceException:
                    m.UpdatePropertyCache(strict=False)
            except ValueError:
                print(f"Problem determining bonds from {molname}")
                return False

        if m:
            # TODO remove extension in the correct way
            molname = os.path.basename(os.path.normpath(os.path.splitext(filename)[0]))
            success = self.analyse(m, molname, mycharge)

        return success

    def write_pdb(self, filenm:str, elements:list, coords:list):
        with open(filenm, "w") as outf:
            outf.write("MODEL        1\n")
            for i in range(len(elements)):
                outf.write("ATOM      1 %2s   UNK     1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % ( elements[i], coords[i][0], coords[i][1], coords[i][2], elements[i] ) )
            outf.write("ENDMDL\n")
            outf.write("TER\n")

    def from_coords_elements(self, molname:str, elements:list,
                             coords:list, mycharge:int):
        if len(coords) != len(elements):
            print("Inconsistent input: There are %d coordinates but %d elements for %s." % ( len(coords), len(elements), molname ))
            return False
        elif len(coords) == 0:
            print("No coordinates for %s" % molname)
            return False
        atomnumber = len(elements)
        success    = False
        self.atoms = []
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, f'mol-{molname}.pdb')
            # use path
            self.write_pdb(path, elements, coords)
            success = self.read(molname, path, "pdb", mycharge)
            os.unlink(path)
            os.rmdir(tmp)

        return success
        
    def from_smiles(self, smiles:str, addH:bool, mycharge:int):
        if self.verbose:
            print("Analyzing %s" % smiles)
        m  = Chem.MolFromSmiles(smiles)
        # Code taken from https://www.rdkit.org/docs/GettingStartedInPython.html
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xf00d # optional random seed for reproducibility
        if addH:
            m2 = Chem.AddHs(m)
            AllChem.EmbedMolecule(m2, params)
            return self.analyse(m2, smiles, mycharge)
        else:
            AllChem.EmbedMolecule(m, params)
            return self.analyse(m, smiles, mycharge)
