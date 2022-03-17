#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import gzip, math, os
from get_mol_dict import *
from molprops import *
from gaff_to_alexandria import *
from get_csv_rows import *
from elements import AtomNumberToAtomName
from dofit import *

debug = False

def method_basis(ww: str):
    if len(ww) >= 3:
        elem = ww.split("/")
        if len(elem) == 2:
            return elem[0], elem[1], "OEP"
    else:
        return ww, ww, ww

def matrix_op(matrix: list, vector: list)->list:
    result = [ 0.0, 0.0, 0.0 ]
    for m in range(3):
        for n in range(3):
            result[m] += matrix[m][n]*vector[n]
    return result

class GaussianReader:
    '''Read the output from a Gaussian calculation to extract
    coordinates, charge, multiplicity, multipole moments and
    more.'''
    def __init__(self, molname:str, basisset:str, verbose:bool):
        # Default stuff
        self.author        = "Spoel2022a"
        self.conformation  = "minimum"
        self.jobtype       = None
        self.qmtype        = "electronic"
        # New compound
        self.mp            = Molprop(molname)
        self.atomname      = []
        self.atomtypes     = None
        self.coordinates   = []
        self.exper         = None
        self.espfitcenter  = []
        self.lastespindex  = None
        self.potential     = []
        self.qEsp          = []
        self.qCM5          = []
        self.qHirshfeld    = []
        self.qMulliken     = []
        self.atomicCenter  = 0
        self.readPotential = 0
        self.molname       = molname
        self.userbasis     = basisset
        self.verbose       = verbose
        # Thermochemistry variables
        self.tcmap         = { "CV": None, "Ezpe": None, "Hcorr": None,
                               "Gcorr": None, "Temp": None, "Method": None,
                               "E0": None, "Scomponent": [], "RotSymNum": 1 }
        
    def rotate_esp_and_add_to_exper(self, infile:str):
        # Added by MMW START
        # Get optimized coords and the last set of coords to determine the rotation matrix
        # Using the rotation matrix we will then rotate the esp grid points 
        natom     = len(self.coordinates)
        if natom < 1:
            print("No atoms")
            return
        if len(self.espfitcenter) < natom:
            print("Not enough fitting centers (%d for %d atoms)." % ( len(self.espfitcenter), natom))
            return
        refcoord  = np.zeros((natom, 3))
        testcoord = np.zeros((len(self.espfitcenter), 3))
        ref_com   = np.zeros(3)
        test_com  = np.zeros(3)
        for i in range(natom):
            refcoord[i] = self.coordinates[i]
            ref_com    += self.coordinates[i]
            if debug:
                print("ref  {}".format(refcoord[i]))
        for i in range(len(self.espfitcenter)):
            testcoord[i] = self.espfitcenter[i]
            # Only include atomic positions in the center of mass
            if i < natom:
                if debug:
                    print("test {}".format(testcoord[i]))
                test_com    += self.espfitcenter[i]
        # Divide COM by number of atoms
        ref_com  /= natom
        test_com /= natom
        # Subtract center of mass
        for i in range(natom):
            refcoord[i]  -= ref_com
            if debug:
                print("ref-com  {}".format(refcoord[i]))
        rmsbefore = 0.0
        for i in range(len(self.espfitcenter)):
            testcoord[i] -= test_com
            if i < natom:
                if debug:
                    print("test-com {}".format(testcoord[i]))
                xdiff       = testcoord[i]-refcoord[i]
                rmsbefore  += xdiff.dot(xdiff)

        # Now do the fit and get the matrix in return
        Rfit = calc_fit_R(natom, refcoord, testcoord)
        if debug:
            print("Rfit {}".format(Rfit))
        # Now rotate the esp points and add them to the experiment
        rmsafter = 0.0
        for i in range(len(self.espfitcenter)):
            newx = Rfit.dot(testcoord[i]) + test_com
            if i < natom:
                if debug:
                    print("newx {}".format(newx))
                xdiff     = newx-refcoord[i]-ref_com
                rmsafter += xdiff.dot(xdiff)
            espx = 100*newx 
            self.exper.add_potential(str(i+1), "Hartree/e", "pm",
                                     str(espx[0]),
                                     str(espx[1]),
                                     str(espx[2]),
                                     self.potential[i])
        rmsb = math.sqrt(rmsbefore/natom)
        rmsa = math.sqrt(rmsafter/natom)
        if debug or (rmsa - rmsb) > 1e-8 or rmsa > 0.1:
            print("RMSD before %g after %g Angstrom in %s" % (rmsb, rmsa, infile))

    def get_revision(self, line):
        self.program = line
        if self.program[len(self.program)-1] == ",":
            self.program = self.program[:-1]
        return 1

    def get_method(self, content, content_index:int):
        methodline = content[content_index].strip()
        while content[content_index+1].find("------") < 0:
            methodline += content[content_index+1].strip()
            content_index += 1
        self.tcmap["Method"], basis, self.jobtype = method_basis(methodline.split()[1])
        if len(self.userbasis) > 0:
            basis = self.userbasis
        return basis, content_index
        
    def get_qmult(self, line):
        words = line.split()
        if len(words) == 6:
            self.mp.add_prop("charge", words[2])
            self.mp.add_prop("multiplicity", words[5])
        return 1

    def get_standard_orientation(self, content, content_index:int) -> int:
        c = content_index+5
        self.atomname.clear()
        self.coordinates.clear()
        while content[c].strip().find("--------------------") < 0:
            words = content[c].strip().split()
            if len(words) == 6:
                atomnumber = int(words[1])
                self.atomname.append(AtomNumberToAtomName(atomnumber))
                self.coordinates.append([ float(words[3]), float(words[4]), float(words[5])])
            else:
                break
            c += 1
        return c-content_index

    def get_energy(self, line:str) -> int:
        words = line.split()
        if len(words) >= 8:
            try:
                self.tcmap["E0"] = float(words[4])
            except ValueError:
                print("Do not understand energy in line '%s'" % line)
        return 1

    def get_dipole(self, content, content_index:int) -> int:
        words = content[content_index+1].strip().split()
        if len(words) >= 6 and self.exper:
            mux     = float(words[1])
            muy     = float(words[3])
            muz     = float(words[5])
            average = math.sqrt(mux**2 + muy**2 + muz**2)
            self.exper.dipoles.clear()
            self.exper.add_dipole(self.qmtype, "D", self.tcmap["Temp"], average,
                                  0.0, mux, muy, muz)
        return 2

    def get_quadrupole(self, content, content_index:int) -> int:
        words = []
        for m in range(2):
            words += content[content_index+1+m].strip().split()
        self.exper.quadrupole.clear()
        self.exper.add_quadrupole(self.qmtype, "B", self.tcmap["Temp"],
                                  words[1], words[3], words[5],
                                  words[7], words[9], words[11])
        return 3
        
    def get_octupole(self, content, content_index:int) -> int:
        elem = {}
        for m in range(3):
            words = content[content_index+1+m].strip().split()
            for w in range(0, len(words), 2):
                windex = words[w][:-1]
                elem[windex.lower()] = words[w+1]
        self.exper.octupole.clear()
        self.exper.add_octupole(self.qmtype, "D.Angstrom2", self.tcmap["Temp"],
                                elem["xxx"], elem["xxy"], elem["xxz"],
                                elem["xyy"], elem["xyz"], elem["xzz"],
                                elem["yyy"], elem["yyz"], elem["yzz"],
                                elem["zzz"])
        return 4

    def get_hexadecapole(self, content, content_index:int) -> int:
        elem = {}
        # This funny map is needed since Gaussian stores non-standard
        # components in the tensor.
        # https://en.wikipedia.org/wiki/Symmetric_tensor
        hdmap = { "yyyx": "xyyy", "zzzx": "xzzz", "zzzy": "yzzz",
                  "yyxz": "xyyz", "zzxy": "xyzz" }
        elem  = {}
        for m in range(4):
            words = content[content_index+1+m].strip().split()
            for w in range(0, len(words), 2):
                windex = words[w][:-1].lower()
                if windex in hdmap:
                    windex = hdmap[windex]
                elem[windex] = words[w+1]
        self.exper.hexadecapole.clear()
        self.exper.add_hexadecapole(self.qmtype, "D.Angstrom3", self.tcmap["Temp"],
                                    elem["xxxx"], elem["xxxy"], elem["xxxz"],
                                    elem["xxyy"], elem["xxyz"], elem["xxzz"],
                                    elem["xyyy"], elem["xyyz"], elem["xyzz"],
                                    elem["xzzz"], elem["yyyy"], elem["yyyz"],
                                    elem["yyzz"], elem["yzzz"], elem["zzzz"])
        return 5

    def get_polarizability(self, line:str) -> int:
        words = line.split()
        if len(words) == 8:
            average = (float(words[2])+float(words[4])+float(words[7]))/3
        self.exper.polarisability.clear()
        self.exper.add_polarisability(self.qmtype, "Bohr3", self.tcmap["Temp"], 
                                      average, 0.0, 
                                      words[2], words[4], words[7],
                                      words[3], words[5], words[6])
        return 1

    def get_esp_centers(self, content, content_index:int) -> int:
        self.espfitcenter.clear()
        acenter  = "Atomic Center"
        ecenter  = "ESP Fit Center"
        espindex = 1
        c        = content_index
        while (content[c].find(acenter) >= 0 or
               content[c].find(ecenter) >= 0):
            line    = content[c].strip()
            words   = line.split()
            mywords = [ None, None, None, None ]
            if content[c].find(acenter) >= 0:
                if len(words) == 8:
                    mywords = [ words[2], words[5], words[6], words[7] ]
                elif len(words) < 8:
                    mywords[0] = line[14:18]
                    mywords[1] = line[25:35]
                    mywords[2] = line[35:45]
                    mywords[3] = line[45:55]
            elif content[c].find(ecenter) >= 0:
                if len(words) == 9:
                    mywords = [ words[3], words[6], words[7], words[8] ]
                elif len(words) < 9:
                    mywords[0] = line[14:19]
                    mywords[1] = line[26:36]
                    mywords[2] = line[36:46]
                    mywords[3] = line[46:56]
            else:
                print("Do not understand line '%s' in %s" % (line, infile))
                return 0
            if self.verbose and False:
                print("found these words {}".format(mywords))
                
            try:
                # The index field may overflow if there are more than
                # 9999 points. In that case Fortran prints ****. We
                # try to circumvent that problem here.
                if mywords[0].find("****") >= 0:
                    espindex += 1
                else:
                    if espindex != int(mywords[0]):
                        print("Expected espindex %d, found %d" % ( espindex, int(mywords[0])))
                espindex += 1
                mycoords  = [ float(mywords[1]), float(mywords[2]), float(mywords[3]) ]
                self.espfitcenter.append(mycoords)
            except ValueError:
                print("Do not understand line '%s' in %s." % (line, infile))
                return 0
            c += 1
        return c-content_index

    def get_esp_potential(self, content, content_index:int) -> int:
        nesppoints = len(self.espfitcenter)
        self.potential.clear()
        cindex = 0
        for dd in range(nesppoints):
            cindex = content_index+2+dd
            if self.verbose and False:
                print("dd %d cindex %d" % ( dd, cindex ))
            if cindex >= len(content):
                print("Inconsistency. Expected %d esp points but found only %d in %s" % ( nesppoints, len(content)-cindex-1, infile))
                return 0
            thisline = content[cindex].strip()
            words    = thisline.split()
            if len(words) == 3 and words[1] in [ "Atom", "Fit" ]:
                self.potential.append(float(words[2]))
            elif thisline.find("----------") >= 0:
                break
        return cindex-content_index
        
    def get_esp_charges(self, content, content_index:int) -> int:
        self.qEsp.clear()
        for c in range(content_index+3, content_index+3+len(self.coordinates)):
            words = content[c].strip().split()
            if len(words) == 3:
                try:
                    self.qEsp.append(words[2])
                except ValueError:
                    print("No charge on this line '%s'" % content[c].strip())
                    break
        return len(self.coordinates)+3

    def get_mulliken_charges(self, content, content_index:int) -> int:
        self.qMulliken.clear()
        for c in range(content_index+2, content_index+2+len(self.coordinates)):
            words = content[c].strip().split()
            if len(words) == 3:
                try:
                    self.qMulliken.append(words[2])
                except ValueError:
                    print("No charge on this line '%s'" % content[c].strip())
                    break
        return len(self.coordinates)+2
 
    def get_cm5_hf_charges(self, content, content_index:int) -> int:
        self.qCM5.clear()
        self.qHirshfeld.clear()
        for c in range(content_index+2, content_index+2+len(self.coordinates)):
            words = content[c].strip().split()
            if len(words) == 8:
                try:
                    self.qHirshfeld.append(words[2])
                    self.qCM5.append(words[7])
                except ValueError:
                    print("No charge on this line '%s'" % content[c].strip())
                    break
        return 2+len(self.coordinates)

    def get_entropy(self, content, content_index:int) -> int:
        words = content[content_index+2].strip().split()
        if len(words) == 4 and words[0] == "Total":
            self.tcmap["CV"] = float(words[2])
        words = content[content_index+4].strip().split()
        if len(words) == 4 and words[0] == "Translational":
            self.tcmap["Scomponent"].append(float(words[2]))
        words = content[content_index+5].strip().split()
        if len(words) == 4 and words[0] == "Rotational":
            self.tcmap["Scomponent"].append(float(words[2]))
        words = content[content_index+6].strip().split()
        if len(words) == 4 and words[0] == "Vibrational":
            self.tcmap["Scomponent"].append(float(words[2]))
        return 6
 
    def add_atoms(self, g2a) -> bool:
        if len(self.atomtypes) != len(self.coordinates):
            print("Found %d atomtype for %d coordinates in %s" % ( len(self.atomtypes), len(self.coordinates), infile))
            return False
        for i in range(len(self.atomtypes)):
            alextype = g2a.rename(self.atomtypes[i])
            qmap = {}
            if len(self.qEsp) == len(self.atomtypes):
                qmap["ESP"] = self.qEsp[i]
            if len(self.qCM5) == len(self.atomtypes):
                qmap["CM5"] = self.qCM5[i]
            if len(self.qHirshfeld) == len(self.atomtypes):
                qmap["Hirshfeld"] = self.qHirshfeld[i]
            if len(self.qMulliken) == len(self.atomtypes):
                qmap["Mulliken"] = self.qMulliken[i]
            # Convert Angstrom to pm
            myx = 100*self.coordinates[i]
            self.exper.add_atom(self.atomname[i], alextype, i+1, "pm", 
                                myx[0], myx[1], myx[2], qmap)
        return True

    def interpret_gauss(self, content:list, infile:str) -> Molprop:
        '''Interpret the content of a Gaussian log file and put the
        contents in a molprop structure for later storage in an XML
        file.'''

        self.tcmap["Temp"] = 0
        content_index      = 0
        atomicCenter       = 0
        readPotential      = 0
        while content_index < len(content):
            line = content[content_index].strip()
            if line.find("Revision") > 0:
                content_index += self.get_revision(line)
                
            elif line.find("#P") >= 0 and None == self.exper:
                basis, content_index = self.get_method(content, content_index)
                self.exper = Experiment("Theory", self.author, self.program, self.tcmap["Method"], basis,
                                        self.conformation, self.jobtype, infile)
                                   
            elif line.find("Multiplicity") >= 0:
                content_index += self.get_qmult(line) 
                
            elif (line.find("Standard orientation") >= 0
                  and line.find("Standard orientation") < 23):
                content_index += self.get_standard_orientation(content, content_index)
                
            elif line.find("SCF Done:") >= 0:
                content_index += self.get_energy(line)
                
            elif line.find("Dipole moment") >= 0:
                content_index += self.get_dipole(content, content_index)
                
            elif line.find("Traceless Quadrupole moment") >= 0:
                content_index += self.get_quadrupole(content, content_index)
                
            elif line.find("Octapole moment") >= 0:
                content_index += self.get_octupole(content, content_index)
                
            elif line.find("Hexadecapole moment") >= 0:
                content_index += self.get_hexadecapole(content, content_index)
                
            elif line.find("Exact polarizability") >= 0:
                content_index += self.get_polarizability(line)
                
            elif line.find("Atomic Center") >= 0 and atomicCenter < 2:
                delta_c = self.get_esp_centers(content, content_index)
                if delta_c == 0:
                    return None
                atomicCenter += 1
                if delta_c != len(self.espfitcenter):
                    print("delta_c = %d len(espfitcenter) = %d" % ( delta_c, len(self.espfitcenter)))
                content_index += delta_c

            elif line.find("Potential          X             Y             Z") >= 0 and readPotential < 2:
                content_index += self.get_esp_potential(content, content_index)
                readPotential += 1
                
            elif line.find("Charges from ESP fit") >= 0:
                content_index += self.get_esp_charges(content, content_index)
                
            elif line.find("Mulliken charges:") >= 0:
                content_index += self.get_mulliken_charges(content, content_index)
        
            elif (line.find("Hirshfeld charges") >= 0 and 
                  line.find("CM5 charges") >= 0):
                content_index += self.get_cm5_hf_charges(content, content_index)
              
            elif line.find(" This molecule is ") >= 0:
                words = content[content_index+1].strip().split()
                tcmap["RotSymNum"] = int(words[3])
                content_index += 1

            elif line.find("Zero-point correction=") >= 0:
                words = line.split()
                self.tcmap["Ezpe"] = float(words[2])
                content_index += 1
                
            elif line.find("Thermal correction to Enthalpy=") >= 0:
                words = line.split()
                self.tcmap["Hcorr"] = float(words[4])
                content_index += 1
                
            elif line.find("Thermal correction to Gibbs Free Energy=") >= 0:
                words = line.split()
                self.tcmap["Gcorr"] = float(words[6])
                content_index += 1
                
            elif line.find("CV") >= 0:
                content_index += self.get_entropy(content, content_index)
                
            elif line.find("Temperature=") >= 0 and line.find("Pressure=") >= 0:
                words  = line.split()
                self.tcmap["Temp"] = float(words[1])
                content_index += 1
                
            else:
                # Will check for thermochemistry now
                # This has to be the last else!
                methodmap = { "CBS-QB3 (0 K)": ["CBS-QB3", 3], "G2(0 K)":["G2", 3], 
                              "G3(0 K)": ["G3", 2], "G4(0 K)":["G4", 2],
                              "W1BD (0 K)":["W1BD", 3], "W1U  (0 K)":["W1U", 3] }
                for tc in methodmap.keys():
                    if line.find(tc) >= 0:
                        self.tcmap["Method"] = methodmap[tc][0]
                        self.tcmap["E0"]     = float(line.split()[methodmap[tc][1]])
                content_index += 1
                
        # End of the big while loop!
        
        # See if we gathered any thermochemistry data
        if self.verbose:
            print("tcmap {}".format(self.tcmap))
        if None != self.tcmap["Temp"] and None != self.tcmap["Method"]:
            leveloftheory =  self.tcmap["Method"] + "/" + basis
            ahof = AtomicHOF(leveloftheory, self.tcmap["Temp"], self.verbose)
            self.exper.extract_thermo(self.tcmap, self.atomname, ahof)
            weight, numb_atoms, formula, multiplicity, self.atomtypes, bonds_dict = get_info_from_coords_elements(self.atomname, self.coordinates)
            if None != weight and None != formula:
                g2a      = GaffToAlexandria()
            if None != self.tcmap["E0"] and None != self.tcmap["Ezpe"]:
                self.tcmap["Temp"]   = 0
                ahof = AtomicHOF(leveloftheory, self.tcmap["Temp"], self.verbose)
                eHF  = compute_dhform(self.tcmap["E0"], self.atomtypes, g2a, ahof,
                                      leveloftheory, self.tcmap["Temp"])
                if debug:
                    print("Computed eHF = %g" % eHF)
                self.exper.add_energy("DeltaE0", "Hartree", 0.0, "gas", eHF)
                if not self.add_atoms(g2a):
                    return None
                self.mp.add_prop("mass", str(weight))
                self.mp.add_prop("formula", formula)
                for index_atom in bonds_dict:
                    for index_neighbour in bonds_dict[index_atom]:
                        self.mp.add_bond(index_atom, index_neighbour, bonds_dict[index_atom][index_neighbour]["bond_order"])
                Angstrom = "A"
                pm       = "pm"
                if self.exper:
                    if len(self.potential) != len(self.espfitcenter):
                        print("Found %d potentials for %d centers in %s" % ( len(self.potential), len(self.espfitcenter), infile))
                        return None
                # Now we have to rotate the electrostatic potential grid points from 
                # the standard orientation to the input orientation.
                if len(self.espfitcenter) > 0:
                    self.rotate_esp_and_add_to_exper(infile)
                self.mp.add_experiment(self.exper)
                return self.mp
            else:
                print("tcmap2 {}".format(tcmap))
        return None
    
    def read(self, infile:str) -> Molprop:
        if not os.path.exists(infile):
            print("No such file " + infile)
            return None
        try:
            with gzip.open(infile, "rt") as inf:
                return self.interpret_gauss(inf.readlines(), infile)
        except gzip.BadGzipFile:
            try:
                with open(infile, "r") as inf:
                    return self.interpret_gauss(inf.readlines(), infile)
            except:
                print("Something fishy with " + infile)
        return None
        
def read_gaussian_log(infile:str, molname:str, basisset:str, verbose:bool) -> Molprop:
    gr = GaussianReader(molname, basisset, verbose)
    return gr.read(infile)
