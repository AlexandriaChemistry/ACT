#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import os, sys, tempfile, xmltodict
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdDetermineBonds

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
            if hasattr(abtype, bondtypes):
                if isinstance(abtype[bondtypes], list):
                    for bt in abtype[bondtypes]['bondtype']:
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
                        print(f"{matchPat}")
                        for k in range(len(matchPat[0])):
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
                        if not (ai == NOTSET and aj == NOTSET):
                            for ab in abe[i][bondtypes]:
                                if ((ab.ai == ai and ab.aj == aj) or
                                    (ab.ai == aj and ab.aj == ai)):
                                    self.bonds[b] = ab.order
                                    break

    def analyse(self, mol, molname:str, charge=None)->bool:
        self.title      = molname
        self.mol_weight = 0
        self.numb_atoms = len(mol.GetAtoms())
        self.formula    = Chem.rdMolDescriptors.CalcMolFormula(mol)
        if None != charge:
            self.charge     = charge
        if self.charge != Chem.GetFormalCharge(mol):
            print(f"Warning: user passed charge {charge} for {molname} but RDKit says charge is {Chem.GetFormalCharge(mol)}")
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
            self.bonds[(ai+1, aj+1)] = order

        # Finally, look for special cases
        self.lookUpSpecial(mol, False)

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
            rdDetermineBonds.DetermineBonds(m)
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
