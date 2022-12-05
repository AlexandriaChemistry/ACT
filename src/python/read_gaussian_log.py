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
            return elem[0], elem[1], "Opt"
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
    def __init__(self, molname:str, basisset:str, verbose:bool, coordinate_set:int):
        # Default stuff
        self.author        = "Spoel2022a"
        self.conformation  = "minimum"
        self.jobtype       = None
        self.qmtype        = "electronic"
        self.coordset      = coordinate_set
        # New compound
        self.mp            = Molprop(molname)
        self.charge        = 0
        self.multiplicity  = 1
        self.atomname      = []
        self.atomtypes     = None
        self.coord_archive = []
        self.forces        = []
        self.exper         = None
        self.espfitcenter  = []
        self.lastespindex  = None
        self.potential     = []
        self.qEsp          = []
        self.qCM5          = []
        self.qHirshfeld    = []
        self.qMulliken     = []
        self.frequencies   = []
        self.intensities   = []
        self.atomicCenter  = 0
        self.readPotential = 0
        self.molname       = molname
        self.userbasis     = basisset
        self.verbose       = verbose
        # Thermochemistry variables
        self.tcmap         = { "CV": None, "Ezpe": None, "Hcorr": None,
                               "Gcorr": None, "Temp": None, "Method": None,
                               "E0": 0, "Scomponent": [], "RotSymNum": 1 }
        self.testvalue = 0
        self.oldvalue = 0
        
    def coordinates(self):
        myindex = -1
        if self.coordset > 0 and self.coordset < len(self.coord_archive):
            myindex  = self.coordset
        return self.coord_archive[myindex]

    def rotate_esp_and_add_to_exper(self, infile:str):
        # Added by MMW START
        # Get optimized coords and the last set of coords to determine the rotation matrix
        # Using the rotation matrix we will then rotate the esp grid points
        mycoords = self.coordinates()
        natom    = len(mycoords)
        if natom < 1:
            print("No atoms")
            return
        if len(self.espfitcenter) < natom:
            print("Not enough fitting centers (%d for %d atoms)." % ( len(self.espfitcenter), natom))
            return
        if len(self.espfitcenter) != len(self.potential):
            print("Inconsistency in %s. %d ESP fit centers and %d potentials" %
                  ( infile, len(self.espfitcenter), len(self.potential)))
            return
        refcoord  = np.zeros((natom, 3))
        testcoord = np.zeros((len(self.espfitcenter), 3))
        ref_com   = np.zeros(3)
        test_com  = np.zeros(3)
        for i in range(natom):
            refcoord[i] = mycoords[i]
            ref_com    += mycoords[i]
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
        Rfit = None
        if natom > 2:
            Rfit = calc_fit_R(natom, refcoord, testcoord)
        if debug or self.verbose:
            print("Rfit {}".format(Rfit))
        # Now rotate the esp points and add them to the experiment
        rmsafter = 0.0
        for i in range(len(self.espfitcenter)):
            oldx = testcoord[i]
            if natom > 2:
                newx = Rfit.dot(oldx) + ref_com #test_com
            else:
                newx = oldx + ref_com
            if i < natom:
                if debug:
                    print("refx {}".format(mycoords[i]))
                    print("oldx {}".format(oldx))
                    print("newx {}".format(newx))
                xdiff     = newx - (refcoord[i]+ref_com)
                rmsafter += xdiff.dot(xdiff)
            espx = 100*newx 
            self.exper.add_potential(str(i+1), "Hartree/e", "pm",
                                     str(espx[0]),
                                     str(espx[1]),
                                     str(espx[2]),
                                     self.potential[i])
        rmsb = math.sqrt(rmsbefore/natom)
        rmsa = math.sqrt(rmsafter/natom)
        if debug or (rmsa - rmsb) > 1e-8 or rmsa > 0.01:
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
            self.charge       = int(words[2])
            self.multiplicity = int(words[5])
        return 1

    def get_forces(self, content, content_index:int) -> int:
        c         = content_index+3
        newforces = []
        while content[c].strip().find("--------------------") < 0:
            words = content[c].strip().split()
            if len(words) == 5:
                atomnumber = int(words[1])
                newforces.append([ float(words[2]), float(words[3]), float(words[4])])
                c += 1
        if len(newforces) > 0:
            self.forces = newforces
        return c-content_index
        
    def get_frequencies(self, content, content_index:int) -> int:
        myc = content_index
        self.frequencies = []
        self.intensities = []
        while len(content[myc].strip()) > 0:
            line = content[myc].strip()
            words = line.split()
            if line.find("Frequencies --") >= 0:
                for fff in range(2, len(words)):
                    self.frequencies.append(words[fff])
            elif line.find("IR Inten    --") >= 0:
                for iii in range(3, len(words)):
                    self.intensities.append(words[iii])
            myc += 1
        return myc-content_index
        
    def get_standard_orientation(self, content, content_index:int) -> int:
        c = content_index+5
        newatomname    = []
        newcoordinates = []
        while content[c].strip().find("--------------------") < 0:
            words = content[c].strip().split()
            if len(words) == 6:
                atomnumber = int(words[1])
                newatomname.append(AtomNumberToAtomName(atomnumber))
                newcoordinates.append([ float(words[3]), float(words[4]), float(words[5])])
            elif len(words) == 5:
                atomnumber = int(words[1])
                newatomname.append(AtomNumberToAtomName(atomnumber))
                newcoordinates.append([ float(words[2]), float(words[3]), float(words[4])])
            else:
                break
            c += 1
        if debug:
            print("Found %d new atoms" % len(newatomname))
        if len(newatomname) > 0:
            self.coord_archive.append(newcoordinates)
            self.atomname    = newatomname
        return c-content_index

    def get_energy(self, line:str) -> int:
        words = line.split()
        if len(words) >= 8:
            try:
                self.testvalue = float(words[4])
                self.oldvalue = float(self.tcmap["E0"])
                if self.testvalue < self.oldvalue:
                    self.tcmap["E0"] = self.testvalue
#                self.tcmap["E0"] = float(words[4])
            except ValueError:
                print("Do not understand energy in line '%s'" % line)
        if line.find("CCSD(T)=") >= 0:
            self.testvalue = float(line.partition("CCSD(T)=")[2].partition("\\")[0])
            self.oldvalue = float(self.tcmap["E0"])
        if self.testvalue < self.oldvalue:
            self.tcmap["E0"] = self.testvalue
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
                if None != mywords[0] and None != mywords[1] and None != mywords[2] and None != mywords[3]:
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

    def natom(self):
        return len(self.coord_archive[-1])
        
    def get_esp_charges(self, content, content_index:int) -> int:
        self.qEsp.clear()
        
        for c in range(content_index+3, content_index+3+self.natom()):
            words = content[c].strip().split()
            if len(words) == 3:
                try:
                    self.qEsp.append(words[2])
                except ValueError:
                    print("No charge on this line '%s'" % content[c].strip())
                    break
        return self.natom()+3

    def get_mulliken_charges(self, content, content_index:int) -> int:
        self.qMulliken.clear()
        for c in range(content_index+2, content_index+2+self.natom()):
            words = content[c].strip().split()
            if len(words) == 3:
                try:
                    self.qMulliken.append(words[2])
                except ValueError:
                    print("No charge on this line '%s'" % content[c].strip())
                    break
        return self.natom()+2
 
    def get_cm5_hf_charges(self, content, content_index:int) -> int:
        self.qCM5.clear()
        self.qHirshfeld.clear()
        for c in range(content_index+2, content_index+2+self.natom()):
            words = content[c].strip().split()
            if len(words) == 8:
                try:
                    self.qHirshfeld.append(words[2])
                    self.qCM5.append(words[7])
                except ValueError:
                    print("No charge on this line '%s'" % content[c].strip())
                    break
        return 2+self.natom()

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
        if len(self.atomtypes) != self.natom():
            print("Found %d atomtype for %d coordinates in %s" % ( len(self.atomtypes), self.natom(), infile))
            return False
        for i in range(len(self.atomtypes)):
            alextype = g2a.rename(self.atomtypes[i])
            qmap = {}
            if len(self.qEsp) == len(self.atomtypes):
                qmap["ESP"] = self.qEsp[i]
            if False:
                if len(self.qCM5) == len(self.atomtypes):
                    qmap["CM5"] = self.qCM5[i]
                if len(self.qHirshfeld) == len(self.atomtypes):
                    qmap["Hirshfeld"] = self.qHirshfeld[i]
                if len(self.qMulliken) == len(self.atomtypes):
                    qmap["Mulliken"] = self.qMulliken[i]
            # If no forces are found, the fc_unit will be None
            # and no forces will be written.
            ff      = [ 0.0, 0.0, 0.0 ]
            fc_unit = None
            if len(self.forces) == len(self.atomtypes):
                ff      = self.forces[i]
                fc_unit = "Hartree/Bohr"
            # Convert Angstrom to pm
            mycoords = self.coordinates()
            self.exper.add_atom(self.atomname[i], alextype, i+1, "pm", 
                                100*mycoords[i][0], 
                                100*mycoords[i][1], 
                                100*mycoords[i][2], 
                                fc_unit, ff[0], ff[1], ff[2], qmap)
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
                myForce = self.jobtype == "Opt"
                self.exper = Experiment("Theory", self.author, self.program, self.tcmap["Method"], basis,
                                        self.conformation, self.jobtype, infile, useForces=myForce)
                content_index += 1
            elif line.find("Multiplicity") >= 0:
                content_index += self.get_qmult(line) 
                
            elif line.find("Standard orientation") >= 0:
                if debug:
                    print("Found new coords at line %d" % content_index)
                content_index += self.get_standard_orientation(content, content_index)
            
            elif line.find("Forces (Hartrees/Bohr)") >= 0:
                content_index += self.get_forces(content, content_index)

            elif line.find("SCF Done:") >= 0:
                content_index += self.get_energy(line)
            elif line.find("\CCSD=") >= 0 or line.find("\CCSD(T)=") >= 0:
                nextline = content[content_index + 1].strip()
                stringline = line + nextline
                content_index += self.get_energy(stringline)
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
                    print("Cannot find ESP centers")
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

            elif line.find("Harmonic frequencies") >= 0:
                content_index += self.get_frequencies(content, content_index)
                
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
            print("tcmap1 {}".format(self.tcmap))
        if None != self.tcmap["Temp"] and None != self.tcmap["Method"]:
            leveloftheory =  self.tcmap["Method"] + "/" + basis
            ahof = AtomicHOF(leveloftheory, self.tcmap["Temp"], self.verbose)
            self.exper.extract_thermo(self.tcmap, self.atomname, ahof)
            
            md = MoleculeDict()
            if not md.from_coords_elements(self.atomname, self.coordinates(), "alexandria"):
                print("Cannot deduce weight or formula from %s" % infile)
                return None
            self.atomtypes = []
            for atom in md.atoms:
                self.atomtypes.append(md.atoms[atom]["atomtype"])
            g2a      = GaffToAlexandria()
            if None != self.tcmap["E0"] and None != self.tcmap["Ezpe"]:
                self.tcmap["Temp"]   = 0
                ahof = AtomicHOF(leveloftheory, self.tcmap["Temp"], self.verbose)
                lots = []
                charges = []
                for i in range(len(self.atomname)):
                    lots.append(leveloftheory)
                    charges.append(0.0)
                eHF  = compute_dhform(self.tcmap["E0"], self.atomname,
                                      g2a, ahof, lots, charges,
                                      self.tcmap["Temp"])
                if debug:
                    print("Computed eHF = %g" % eHF)
                self.exper.add_energy("DeltaE0", "Hartree", 0.0, "gas", eHF)
                # Now we have to rotate the electrostatic potential grid points 
                # from the standards orientation to the input orientation.
                if len(self.espfitcenter) > 0:
                    self.rotate_esp_and_add_to_exper(infile)
                # Now add frequencies and intensities if they look OK.
                if len(self.frequencies) > 0 and len(self.frequencies) == len(self.intensities):
                    self.exper.frequencies = self.frequencies
                    self.exper.intensities = self.intensities
                # We can only add the atoms after selecting the coordinates.
                # This is likely a peculiarity of the Alexandria library.
                if not self.add_atoms(g2a):
                    print("Cannot add the atoms or atomtypes")
                    return None
                if len(md.inchi) > 0:
                    print("Empy InChi")
                    return None
                frag = Fragment(md.inchi, self.charge, self.multiplicity, 1, range(1,1+len(self.atomtypes)), md.mol_weight, md.formula)
                self.mp.add_fragment(frag)
                for index_tuple in md.bonds:
                    self.mp.add_bond(index_tuple[0], index_tuple[1], md.bonds[index_tuple])
                Angstrom = "A"
                pm       = "pm"
                if self.exper:
                    if len(self.potential) != len(self.espfitcenter):
                        print("Found %d potentials for %d centers in %s" % ( len(self.potential), len(self.espfitcenter), infile))
                        return None
                self.mp.add_experiment(self.exper)
                return self.mp
            else:
                print("tcmap2 {}".format(tcmap))
        else:
            print("No temperature or method? tcmap3 {}".format(tcmap))
        return None
    
    def read(self, infile:str) -> Molprop:
        if not os.path.exists(infile):
            print("No such file " + infile)
            return None
        try:
            with gzip.open(infile, "rt") as inf:
                mp = self.interpret_gauss(inf.readlines(), infile)
                if None == mp:
                    sys.exit("Could not read %s" % infile)
                else:
                    return mp
        except gzip.BadGzipFile:
            try:
                with open(infile, "r") as inf:
                    return self.interpret_gauss(inf.readlines(), infile)
            except:
                print("Something fishy with " + infile)
        return None
        
def read_gaussian_log(infile:str, molname:str, basisset:str, verbose:bool, coordinate_set:int) -> Molprop:
    gr = GaussianReader(molname, basisset, verbose, coordinate_set)
    mp = gr.read(infile)
    if None == mp:
        sys.exit("Could not read %s" % infile)
    else:
        return mp
