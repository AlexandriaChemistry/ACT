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

def check_bondorder(molecule_dict):
    my_pairs = { "no": "on", "cm": "om", "p5": "om", "s6": "om", "py": "om", "s3": "om", "cz": "n", "s4": "om"  }
    for mp in my_pairs.keys():
        for i in range(1, 1+len(molecule_dict["atoms"])):
            if mp == molecule_dict["atoms"][i]["atomtype"]:
                bondIndex = []
                for (ai,aj) in molecule_dict["bonds"].keys():
                    if ((i == ai and my_pairs[mp] == molecule_dict["atoms"][aj]["atomtype"]) or
                        (i == aj and my_pairs[mp] == molecule_dict["atoms"][ai]["atomtype"])):
                        bondIndex.append((ai, aj))
                if len(bondIndex) >= 2:
                    for b in bondIndex:
                        molecule_dict["bonds"][b] = 1.5

def get_mol_dict(filename, fileformat, forcefield=None):
    molecule_dict = {"molecule": {}, "atoms": {}, "bonds": {}}
    obconversion = ob.OBConversion()
    obconversion.SetInFormat(fileformat)
    obmol = ob.OBMol()
    
    notatend   = obconversion.ReadFile(obmol,filename)
    title      = obmol.GetTitle()
    mol_weight = obmol.GetMolWt()
    numb_atoms = obmol.NumAtoms()
    formula    = obmol.GetFormula()
    charge     = obmol.GetTotalCharge()
    obmol.AssignTotalChargeToAtoms(charge)
    obmol.SetAromaticPerceived(False) 
    obconversion.SetOutFormat("inchi")
    inchi      = obconversion.WriteString(obmol).strip()
    obconversion.SetOutFormat("inchikey")
    inchikey   = obconversion.WriteString(obmol).strip()

    molecule_dict["molecule"].update({"title": title, "mol_weight": mol_weight,
                                      "numb_atoms": numb_atoms, "formula": formula,
                                       "charge": charge,  "multiplicity": None,
                                        "inchi": inchi, "inchikey": inchikey})

    if forcefield:
        ff = ob.OBForceField.FindForceField(forcefield)
        if not ff:
            sys.exit("No OpenBabel support for force field %s" % forcefield)
        if not ff.Setup(obmol):
            print("Could not setup the force field %s for %s" % ( forcefield, filename))
        else:
            if not ff.GetAtomTypes(obmol):
                print("Could not get atomtypes from force field %s for %s" % ( forcefield, filename))

    g2a      = GaffToAlexandria()
    # Add the atoms
    for atom in ob.OBMolAtomIter(obmol):
        index      = atom.GetIdx()
        atomtype   = "Z"
        ffatomtype = "FFAtomType"
        if forcefield and atom.HasData(ffatomtype):
            atp = atom.GetData(ffatomtype)
            if atp:
                atomtype = g2a.rename(atp.GetValue())
        X = atom.GetX()
        Y = atom.GetY()
        Z = atom.GetZ()
        atomic_number = atom.GetAtomicNum()
        mass = atom.GetExactMass()
        molecule_dict["atoms"].update({index: {}})
        molecule_dict["atoms"][index] = {"atomic_number": atomic_number,
                                         "atomtype": atomtype, "mass": mass,
                                          "X": X, "Y": Y, "Z": Z}
    
    # Add the bonds
    for bond in ob.OBMolBondIter(obmol):
        bbb = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        molecule_dict["bonds"][bbb] = bond.GetBondOrder()
        if bond.IsAromatic() and forcefield == "alexandria":
            molecule_dict["bonds"][bbb] = 1.5
    if forcefield == "alexandria":
        check_bondorder(molecule_dict)
        
    if debug:
        print(molecule_dict)
        
    return molecule_dict

def get_info_from_coords_elements(elements, coords, forcefield="alexandria"):
  if len(coords) != len(elements):
    print("Inconsistent input: There are %d coordinates but %d elements." % ( len(coords), len(elements) ))
    return None, None, None, None, None, None
  elif len(coords) == 0:
    return None, None, None, None, None, None
  atomnumber = len(elements)
  xyzstring = ("%d\nCoordinates\n" % atomnumber)
  for i in range(len(elements)):
      xyzstring += ("%2s%22.3f%22.3f%22.3f\n" % ( elements[i], coords[i][0], coords[i][1], coords[i][2] ))

  obConversion = ob.OBConversion()
  obConversion.SetInFormat("xyz")
  obmol = ob.OBMol()

  obConversion.ReadString(obmol, xyzstring)

  numb_atoms = obmol.NumAtoms()
  weight = obmol.GetMolWt()
  formula = obmol.GetFormula()
  multiplicity = obmol.GetTotalSpinMultiplicity()

  ff = ob.OBForceField.FindForceField(forcefield)
  ff.Setup(obmol)
  ff.GetAtomTypes(obmol)

  atomtypes = []
  bonds_dict = {}
  for obatom in ob.OBMolAtomIter(obmol):
    coordinates_element = []
    type = obatom.GetData("FFAtomType")
    atomtypes.append(type.GetValue())
    index_atom = obatom.GetIdx()
    bonds_dict[index_atom]={}

    for neighbour_atom in ob.OBAtomAtomIter(obatom):
      index_neighbour = neighbour_atom.GetIdx()
      if index_neighbour > index_atom:
        bonds_dict[index_atom][index_neighbour]={}
        bond = obatom.GetBond(neighbour_atom)
        bond_order = bond.GetBondOrder()
        bonds_dict[index_atom][index_neighbour]={"bond_order": bond_order}

  for keys in sorted(bonds_dict.keys()):
    has_items =  bool(bonds_dict[keys])
    if has_items != True:
      bonds_dict.pop(keys)

  return weight, numb_atoms, formula, multiplicity, atomtypes, bonds_dict
