import gzip, math, os
from get_mol_dict import *
from molprops import *
from gaff_to_alexandria import *

def method_basis(ww: str):
    if len(ww) >= 3:
        elem = ww.split("/")
        if len(elem) == 2:
            return elem[0], elem[1]
    return "", ""

def interpret_gauss(content:list, infile:str,
                    molname:str, basisset:str, verbose:bool):
    '''Interpret the content of a Gaussian log file and put the
    contents in a molprop structure for later storage in an XML
    file.'''
    # Default stuff
    author       = "Spoel2022a"
    conformation = "minimum"
    jobtype      = "OEP"
    temperature  = 0
    qmtype       = "electronic"
    # New compound
    mp           = Molprop(molname)
    atomname     = []
    coordinates  = []
    exper        = None
    espfitcenter = []
    potential    = []
    # Thermochemistry variables
    tcmap        = { "CV": None, "Ezpe": None, "Hcorr": None,
                     "Gcorr": None, "Temp": None, "Method": None,
                     "E0": None, "Scomponent": [], "RotSymNum": 1 }
    for content_index in range(0, len(content)):
        line = content[content_index].strip()
        if line.find("Revision") > 0:
            program = line
        elif line.find("#P") >= 0 and not exper:
            method, basis = method_basis(line.split()[1])
            if len(basisset) > 0:
                basis = basisset
            exper = Experiment("Theory", author, program, method, basis,
                               conformation, jobtype, infile)
        elif line.find("Multiplicity") >= 0:
            words = line.split()
            if len(words) == 6:
                mp.add_prop("charge", words[2])
                mp.add_prop("multiplicity", words[5])
                c = content_index+1
                while len(content[c].strip()) > 0:
                    words = content[c].strip().split()
                    if len(words) == 4:
                        atomname.append(words[0])
                        coordinates.append([ float(words[1]), float(words[2]), float(words[3])])
                    else:
                        break
                    c += 1
        elif line.find("Dipole moment") >= 0:
            words = content[content_index+1].strip().split()
            if len(words) >= 6 and exper:
                mux     = float(words[1])
                muy     = float(words[3])
                muz     = float(words[5])
                average = math.sqrt(mux**2 + muy**2 + muz**2)
                exper.dipoles.clear()
                exper.add_dipole(qmtype, "D", temperature, average,
                                 0.0, mux, muy, muz)
        elif line.find("Traceless Quadrupole moment") >= 0:
            words = []
            for m in range(2):
                words += content[content_index+1+m].strip().split()
            exper.quadrupole.clear()
            exper.add_quadrupole(qmtype, "B", temperature,
                                 words[1], words[3], words[5],
                                 words[7], words[9], words[11])
        elif line.find("Octapole moment") >= 0:
            elem = {}
            for m in range(3):
                words = content[content_index+1+m].strip().split()
                for w in range(0, len(words), 2):
                    windex = words[w][:-1]
                    elem[windex.lower()] = words[w+1]
            exper.octupole.clear()
            exper.add_octupole(qmtype, "DAngstrom2", temperature,
                               elem["xxx"], elem["xxy"], elem["xxz"],
                               elem["xyy"], elem["xyz"], elem["xzz"],
                               elem["yyy"], elem["yyz"], elem["yzz"],
                               elem["zzz"])
        elif line.find("Hexadecapole moment") >= 0:
            elem = {}
            # This funny map is needed since Gaussian stores non-standard
            # components in the tensor.
            # https://en.wikipedia.org/wiki/Symmetric_tensor
            hdmap = { "yyyx": "xyyy", "zzzx": "xzzz", "zzzy": "yzzz",
                      "yyxz": "xyyz", "zzxy": "xyzz" }
            for m in range(4):
                words = content[content_index+1+m].strip().split()
                for w in range(0, len(words), 2):
                    windex = words[w][:-1].lower()
                    if windex in hdmap:
                        windex = hdmap[windex]
                    elem[windex] = words[w+1]
            exper.hexadecapole.clear()
            exper.add_hexadecapole(qmtype, "DAngstrom3", temperature,
                                   elem["xxxx"], elem["xxxy"], elem["xxxz"],
                                   elem["xxyy"], elem["xxyz"], elem["xxzz"],
                                   elem["xyyy"], elem["xyyz"], elem["xyzz"],
                                   elem["xzzz"], elem["yyyy"], elem["yyyz"],
                                   elem["yyzz"], elem["yzzz"], elem["zzzz"])
        elif line.find("Exact polarizability") >= 0:
            words = line.split()
            if len(words) == 8:
                average = (float(words[2])+float(words[3])+float(words[4]))/3
                exper.polarisability.clear()
                exper.add_polarisability(qmtype, "Bohr3", temperature, 
                                         average, 0.0, 
                                         words[2], words[3], words[4],
                                         words[5], words[6], words[7])
        elif line.find("Atomic Center") >= 0:
            words     = line.split()
            if len(words) < 8:
                print("Do not understand line '%s' in %s." % (line, infile))
                return None
            espindex  = int(words[2])-1
            mycoords  = [ float(words[5]), float(words[6]), float(words[7]) ]
            if espindex >= len(espfitcenter):
                espfitcenter.append(mycoords)
            else:
                espfitcenter[espindex] = mycoords
        elif line.find("ESP Fit Center") >= 0:
            words     = line.split()
            if len(words) < 9:
                print("Do not understand line '%s' in %s" % (line, infile))
                return None
            espindex  = int(words[3])-1
            try:
                mycoords  = [ float(words[6]), float(words[7]), float(words[8]) ]
            except:
                print("'%s' does not contain ESP in %s" % ( line, infile))
                return None
            if espindex >= len(espfitcenter):
                espfitcenter.append(mycoords)
            else:
                espfitcenter[espindex] = mycoords
        elif line.find("Potential          X             Y             Z") >= 0:
            nesppoints = len(espfitcenter)
            for dd in range(nesppoints):
                cindex = content_index+2+dd
                if verbose:
                    print("dd %d cindex %d" % ( dd, cindex ))
                if cindex >= len(content):
                    print("Inconsistency. Expected %d esp points but found only %d in %s" % ( nesppoints, len(content)-cindex-1, infile))
                    return None
                thisline = content[cindex].strip()
                words    = thisline.split()
                if len(words) == 3 and words[1] in [ "Atom", "Fit" ]:
                    if dd >= len(potential):
                        potential.append(float(words[2]))
                    else:
                        if verbose:
                            print(words)
                        potential[dd] = words[2]
                elif thisline.find("----------") >= 0:
                    break
        elif line.find(" This molecule is ") >= 0:
            words = content[content_index+1].strip().split()
            tcmap["RotSymNum"] = int(words[3])
        elif line.find("Zero-point correction=") >= 0:
            words = line.split()
            tcmap["Ezpe"] = float(words[2])
        elif line.find("Thermal correction to Enthalpy=") >= 0:
            words = line.split()
            tcmap["Hcorr"] = float(words[4])
        elif line.find("Thermal correction to Gibbs Free Energy=") >= 0:
            words = line.split()
            tcmap["Gcorr"] = float(words[6])
        elif line.find("CV") >= 0:
            words = content[content_index+2].strip().split()
            if len(words) == 4 and words[0] == "Total":
                tcmap["CV"] = float(words[2])
            words = content[content_index+4].strip().split()
            if len(words) == 4 and words[0] == "Translational":
                tcmap["Scomponent"].append(float(words[2]))
            words = content[content_index+5].strip().split()
            if len(words) == 4 and words[0] == "Rotational":
                tcmap["Scomponent"].append(float(words[2]))
            words = content[content_index+6].strip().split()
            if len(words) == 4 and words[0] == "Vibrational":
                tcmap["Scomponent"].append(float(words[2]))

            
        elif line.find("Temperature=") >= 0 and line.find("Pressure=") >= 0:
            words  = line.split()
            tcmap["Temp"] = float(words[1])
        else:
            # Will check for thermochemistry now
            # This has to be the last else!
            methodmap = { "CBS-QB3 (0 K)": ["CBS-QB3", 3], "G2(0 K)":["G2", 3], 
                      "G3(0 K)": ["G3", 2], "G4(0 K)":["G4", 2],
                      "W1BD (0 K)":["W1BD", 3], "W1U  (0 K)":["W1U", 3] }
            for tc in methodmap.keys():
                if line.find(tc) >= 0:
                    tcmap["Method"] = methodmap[tc][0]
                    tcmap["E0"]     = float(line.split()[methodmap[tc][1]])
    
    exper.extract_thermo(tcmap, atomname)
    weight, numb_atoms, formula, multiplicity, atomtypes, bonds_dict = get_info_from_coords_elements(atomname, coordinates)
    mp.add_prop("mass", str(weight))
    mp.add_prop("formula", formula)
    for index_atom in bonds_dict:
        for index_neighbour in bonds_dict[index_atom]:
            mp.add_bond(index_atom, index_neighbour, bonds_dict[index_atom][index_neighbour]["bond_order"])
    Angstrom = "A"
    pm       = "pm"
    g2a      = GaffToAlexandria()
    if exper:
        if len(potential) != len(espfitcenter):
            print("Found %d potentials for %d centers in %s" % ( len(potential), len(espfitcenter), infile))
            return None
        for i in range(len(potential)):
            exper.add_potential(str(i+1), "Hartree/e", pm,
                                str(100*espfitcenter[i][0]),
                                str(100*espfitcenter[i][1]),
                                str(100*espfitcenter[i][2]),
                                potential[i])
        if len(atomtypes) != len(coordinates):
            print("Found %d atomtype for %d coordinates in %s" % ( len(atomtypes), len(coordinates), infile))
            return None
        for i in range(len(atomtypes)):
            alextype = g2a.rename(atomtypes[i])
            exper.add_atom(atomname[i], alextype, i+1, pm, 
                           (100*coordinates[i][0]),
                           (100*coordinates[i][1]),
                           (100*coordinates[i][2]))
        mp.add_experiment(exper)
        return mp
    else:
        return None
    
def read_gaussian_log(infile:str, molname:str, basisset:str, verbose:bool):
    '''Read the output from a Gaussian calculation to extract
    coordinates, charge, multiplicity, multipole moments and
    more.'''
    if not os.path.exists(infile):
        print("No such file " + infile)
        return None
    try:
        with gzip.open(infile, "rt") as inf:
            content = inf.readlines()
            return interpret_gauss(content, infile, molname, basisset, verbose)
    except gzip.BadGzipFile:
        print("Something fishy with " + infile)
        return None
