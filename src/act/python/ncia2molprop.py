#!/usr/bin/env python3
#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import argparse, sys
#​
# ACT python code
# sys.path.insert(1, "@ACT_PYTHON_LIB_PATH@")
#from /src/python/molprops import *

from molprops import *
from get_mol_dict import *
import openbabel as ob
import lxml
#import mol_csv_api as piskot



int_e = ""
POINT = "SP"
piskot = False
#sel_a_end = 0
#sel_b_end = 0
#charge_a = 0
#charge_b = 9
#multi
selection = ""
def parseArguments():
    desc = '''This script will read an xyz file
    and write an XML file for input to the Alexandria Chemistry Toolkit.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-i", "--infile", help="XYZ file for reading", type=str)
    outf = "test.xml"
    parser.add_argument("-o", "--outfile", help="Output xml file for writing, default "+outf, type=str, default=outf)
    molname = "chemical"
    parser.add_argument("-n", "--molname", help="Molecule name, default "+molname, type=str, default=molname)
    parser.add_argument("-basis", "--basisset", help="Basis set, will override what is in the input file", type=str, default="")
    parser.add_argument("-v", "--verbose", help="write debug statements", action="store_true")
    parser.add_argument("-int", "--int_e", help="interaction energy", type=float, default="")

    parser.add_argument("-charge_a", "--charge_a", help="charge of frag A", type=int, default="")
    parser.add_argument("-charge_b", "--charge_b", help="charge of frag B", type=int, default="")
    
#    parser.add_argument("-mult_b", "--multiplicity_b", help="multiplicity of frag A", type=int, default="")
#    parser.add_argument("-mult_a", "--multiplicity_a", help="multiplicity of frag A", type=int, default="")
    

    parser.add_argument("-sel_a_end", "--sel_a_end", help="end atom of selection A", type=int, default="")
    parser.add_argument("-sel_b_end", "--sel_b_end", help="end of atom selection B", type=int, default="")

    parser.add_argument("-FRAGA", "--FRAGA", help="fragment A xyz", type=str, default="")

    parser.add_argument("-FRAGB", "--FRAGB", help="fragment B xyz", type=str, default="")

    parser.add_argument("-MIN", "--minimum", help="is this a minimum", type=str, default="")
    parser.add_argument("-ID", "--identifier", help="this is a name of a specific geometry within the experiment", type=str, default="")
    parser.add_argument("-ref", "--reference", help="Reference", type=str, default="")
#    parser.add_argument("-ref", "--refe", help="reference", type=str, default="")
    return parser.parse_args()
#​


def read_auxA(infile:str):
    md = MoleculeDict()
    if not md.read(infile, "xyz"):
        sys.exit("Could not read file %s" % infile)
    global weight_B
    weight_B = str(md.mol_weight)
    global form_B
    form_B = str(md.formula)
    global inchiA
    inchiA = str(md.inchi)
    return mp1


def read_auxB(infile:str):
    md = MoleculeDict()
    if not md.read(infile, "xyz"):
        sys.exit("Could not read file %s" % infile)
    global weight_A
    weight_A = str(md.mol_weight)
    global form_A
    form_A = str(md.formula)
    global inchiB
    inchiB = str(md.inchi)

    return mp2





def read_xyz(infile:str) -> Molprop:
    # Create a molecule reader
    md = MoleculeDict()
    if not md.read(infile, "xyz"):
        sys.exit("Could not read file %s" % infile)
    # Now we have read the atoms, coordinates and all that
    # now the energy needs to be read from somewhere too.
    # Also the variables below need to be corrected for each
    # data set.
    reference = args.reference
    program   = "N/A"
    method    = "CCSD(T)"
    basisset  = "CBS"

    POINT = "SP"
#    print(args.minimum)
    personnummer = args.identifier
    if args.minimum.lower() == "yes":
        POINT = "Opt"
        piskot = True
    else:
        piskot = False
    print(f"IS THIS MINIMUM: {piskot} {POINT}")
    exper = Experiment("Theory", reference, program, method, basisset,
                       "minimum", POINT, personnummer, piskot)
#    print(exper)
#    mp = Molprop(md.inchi)
#    wat1 = Molprop("water1")
#    SELA = Molprop("selection_a")
#    SELB = Molprop("selection_b")
    interaction_energy = f'{args.int_e:.10f}'
#    chargeA = args.charge_a
    exper.add_energy("DeltaE0", "kJ/mol", 0.0, "gas", interaction_energy)

    rangeA=[*range(1,args.sel_a_end + 1)]
    rangeB=[*range(args.sel_a_end + 1,args.sel_b_end + 1)]
###Loop for a fragment A

    g2a      = GaffToAlexandria()
    atomtypes = []
    atomname  = []

# Auxiliary structure to rename atom types
#move above    g2a      = GaffToAlexandria()
    # Now create a molprop object
    mp.add_experiment(exper)
    mp.add_prop("mass", str(md.mol_weight))
    mp.add_prop("formula", str(md.formula))

#    print(selection)
#    if selection == "A":
#        global weight_A
#        weight_A = str(md.mol_weight)
#        global form_A
#        form_A = str(md.formula)
#        Afrag = Fragment(md.inchi, args.charge_a, 1, 1, rangeA, weight_A, form_A)
#        mp.add_fragment(Afrag)

#        global ID_A
#        ID_A = 1
#    if selection == "B":
#        global weight_B
#        weight_B = str(md.mol_weight)
#        global form_B
#        form_B = str(md.formula)
#        Bfrag = Fragment(md.inchi, args.charge_b, 1, 1, rangeB, weight_B, form_B)
#        mp.add_fragment(Bfrag)
#        global ID_B
#        ID_B = 2
#    print(str(md.formula))
#    print(str(md.mol_weight))
    for index_tuple in md.bonds:
        mp.add_bond(index_tuple[0], index_tuple[1], md.bonds[index_tuple])

    if selection == "AB":
###generate fragments 
###################

        mp.set_molname(md.inchi)
        Afrag = Fragment(inchiA, args.charge_a, 1, 1, rangeA, weight_A, form_A)
        mp.add_fragment(Afrag)
        Bfrag = Fragment(inchiB, args.charge_b, 1, 1, rangeB, weight_B, form_B)
        mp.add_fragment(Bfrag)
##Now add the atoms        
        atomtypes = []
        atomname  = []
#        print(str(md.atoms))
        for atom in md.atoms:
            atomtypes.append(md.atoms[atom]["atomtype"])
            atomname.append(md.atoms[atom]["atomtype"])
#            print(str(atomtypes))
        for i in range(len(atomtypes)):
            alextype = g2a.rename(atomtypes[i])
            exper.add_atom(atomname[i], alextype, i+1, "pm",
                       100*md.atoms[i+1]["X"],
                       100*md.atoms[i+1]["Y"],
                       100*md.atoms[i+1]["Z"],
                       "Hartree/Bohr", 0, 0, 0, None)



        del Afrag
        del Bfrag
    # Done! Return the result.
    return mp

#    Afrag = Fragment(args.charge_a, args.multiplicity_a, rangeA, weight_A, form_A)
#    Bfrag = Fragment(args.charge_b, args.multiplicity_b, rangeB, weight_B, form_B)  
#    mp.add_fragment(Afrag)
#    mp.add_fragment(Bfrag)



if __name__ == '__main__':
    args  = parseArguments()
    md = MoleculeDict()
    molprops = Molprops()
    molprops.open(args.outfile)
    mp        = None    
    mp = Molprop(args.outfile)
    mp1 = None
    mp2 = None
#    selection = "A"
    read_auxA(args.FRAGA)
    read_auxB(args.FRAGB)
#    SELA = read_xyz(args.FRAGA)
#    inchi = md.inchi
#    print(inchi)
#    if SELA:
#        SELA = read_xyz(args.FRAGA, md.inchi)
#        molprops = Molprops()
#    molprops.add_molecule(SELA, False)
#       iupac = piskot.Molecule.get_iupac(SELA)
#        SELA.add_experiment(exper)
#        SELA.add_prop("mass", str(md.mol_weight))
#        SELA.add_prop("formula", str(md.formula))

#        Amass = str(md.mol_weight)
#        Aform = str(md.formula)
#        print(f"AMASS {Amass}")

#    selection = "B"
#    SELB = read_xyz(args.FRAGB)
#    if SELB:
#        molprops = Molprops()
#    molprops.add_molecule(SELB, False)
#        SELB = read_xyz(args.FRAGB, md.inchi)
#        Bmass = str(md.mol_weight)
#        Bform = str(md.formula)
#        print(f"AMASS {Amass}")
#        print(f"BMASS {Bmass}")


    selection = "AB"
    WHOLE = read_xyz(args.infile)
    molprops.add_molecule(WHOLE, False)
#    mp = read_xyz(args.infile)

#​
#    if mp:
#        molprops = Molprops()
#        molprops.add_molecule(mp, False)
#        molprops.write(args.outfile)
#        print("Wrote %s" % args.outfile)
#        print(f"this is A weight {weight_A}")
#test subject water dimer
##    wat1 = read_xyz("wat1.xyz", "water1")
##    if wat1:
##        molprops = Molprops()
##        molprops.add_molecule(wat1

    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
