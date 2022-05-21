#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/dspoel/ACT
#
import xml.etree.ElementTree as ET
from xml.dom import minidom
from atomic_heat_of_formation import *

debug = False

class Molprops:
    '''
    Class to write molprop files (molecule properties) from the
    Alexandria Chemistry Toolkit
    '''

    def __init__(self):
        self.molecules = ET.Element("molecules")
 
    def add_molecule(self, molecule):
        if len(molecule.fragments) == 0:
            molname = ""
            if "molname" in molecule.properties:
                molname = molecule.properties["molname"]
            print("There are no fragments. Ignoring molecule %s" % molname)
            return
        lastmol = ET.SubElement(self.molecules, "molecule")
        for prop in molecule.properties.keys():
            lastmol.set(prop, molecule.properties[prop])
        if len(molecule.bonds) > 0:
            for bond in molecule.bonds:
                newbond = ET.SubElement(lastmol, "bond")
                newbond.set("ai", bond[0])
                newbond.set("aj", bond[1])
                newbond.set("bondorder", bond[2])
        fragments = ET.SubElement(lastmol, "fragments")
        for frag in molecule.fragments:
            xmlfrag = ET.SubElement(fragments, "fragment")
            frag.add_xml(xmlfrag)

        for ep in molecule.experiments:
            exper = ET.SubElement(lastmol, "experiment")
            fc_unit = None
            if ep.useForce:
                fc_unit = "force_unit"
            for prop in ep.properties:
                exper.set(prop, ep.properties[prop])
            for ene in ep.energies:
                myene = ET.SubElement(exper, "energy")
                for prop in ene.keys():
                    aver = "average"
                    if prop == aver:
                        ccc = ET.SubElement(myene, aver)
                        ccc.text = ene[prop]
                    else:
                        myene.set(prop, ene[prop])
            for dip in ep.dipoles:
                mydip = ET.SubElement(exper, "dipole")
                for myprop in [ "type", "unit", "temperature" ]:
                    if myprop in dip:
                        mydip.set(myprop, dip[myprop])
                for value in [ "average", "error", "x", "y", "z"]:
                    if value in dip:
                       myvalue = ET.SubElement(mydip, value)
                       myvalue.text = dip[value]
            if len(ep.frequencies) > 0:
                myharm = ET.SubElement(exper, "harmonics")
                myharm.set("frequency_unit", "cm^-1")
                myharm.set("intensity_unit", "km/mol")
                for k in range(len(ep.frequencies)):
                    harm = ET.SubElement(myharm, "harmonic")
                    harm.set("frequency", ep.frequencies[k])
                    harm.set("intensity", ep.intensities[k])
            for atom in ep.atoms:
                myatm   = ET.SubElement(exper, "atom")
                testprops = [ "name", "obtype", "atomid", "coord_unit" ]
                if None != fc_unit:
                    testprops.append(fc_unit)
                for myprop in testprops:
                    if myprop in atom and None != atom[myprop]:
                        myatm.set(myprop, atom[myprop])
                myelem = [ "x", "y", "z" ]
                if None != fc_unit and fc_unit in atom and None != atom[fc_unit]:
                    myelem += [ "fx", "fy", "fz" ]
                for coord in myelem:
                    if coord in atom:
                        myx = ET.SubElement(myatm, coord)
                        myx.text = atom[coord]
                qmap = "qmap"
                if qmap in atom:
                    for q in atom[qmap].keys():
                        myq = ET.SubElement(myatm, "q")
                        myq.set("type", q)
                        myq.text = atom[qmap][q]
                        
            for pol in ep.polarisability:
                mypol = ET.SubElement(exper, "polarizability")
                for myprop in ["type", "unit", "temperature"]:
                    if myprop in pol:
                        mypol.set(myprop, pol[myprop]) 
                for value in [ "average", "error", "xx", "yy", "zz", "xy", "xz", "yz"]:
                    if value in pol:
                       myvalue = ET.SubElement(mypol, value)
                       myvalue.text = pol[value]

            for quad in ep.quadrupole:
                myquad = ET.SubElement(exper, "quadrupole")
                for myprop in ["type", "unit", "temperature"]:
                    if myprop in quad:
                        myquad.set(myprop, quad[myprop]) 
                for value in ["xx", "yy", "zz", "xy", "xz", "yz"]:
                    if value in quad:
                       myvalue = ET.SubElement(myquad, value)
                       myvalue.text = quad[value]

            for octu in ep.octupole:
                myoctu = ET.SubElement(exper, "octupole")
                for myprop in ["type", "unit", "temperature"]:
                    if myprop in octu:
                        myoctu.set(myprop, octu[myprop]) 
                for value in ["xxx", "xxy", "xxz", "xyy", "xyz", "xzz", "yyy", "yyz", "yzz", "zzz"]:
                    if value in octu:
                       myvalue = ET.SubElement(myoctu, value)
                       myvalue.text = octu[value]

            for hexadeca in ep.hexadecapole:
                myhexadeca = ET.SubElement(exper, "hexadecapole")
                for myprop in ["type", "unit", "temperature"]:
                    if myprop in hexadeca:
                        myhexadeca.set(myprop, hexadeca[myprop]) 
                for value in ["xxxx", "xxxy", "xxxz", "xxyy", "xxyz", "xxzz", "xyyy", "xyyz", "xyzz", "xzzz", "yyyy", "yyyz", "yyzz", "yzzz", "zzzz"]:
                    if value in hexadeca:
                       myvalue = ET.SubElement(myhexadeca, value)
                       myvalue.text = hexadeca[value]

            for esp in ep.potential:
                myesp = ET.SubElement(exper, "potential")
                for myprop in ["espid", "coord_unit", "potential_unit"]:
                    if myprop in esp:
                        myesp.set(myprop, esp[myprop]) 
                for value in [ "x", "y", "z", "V"]:
                    if value in esp:
                       myvalue = ET.SubElement(myesp, value)
                       myvalue.text = esp[value]

    def write(self, outfile):
        xmlstr      = minidom.parseString(ET.tostring(self.molecules)).toprettyxml(indent="  ")
        with open(outfile, "w") as f:
            f.write(xmlstr)

class Fragment:
    '''
    Class to store molecule "fragments" which in this context 
    means compounds in e.g. a dimer.
    '''
    
    def __init__(self, identifier:str, charge:int, multiplicity:int, symmetry_number: int, atoms:list, mass:float, formula: str):
        self.identifier      = identifier
        self.charge          = charge
        self.multiplicity    = multiplicity
        self.symmetry_number = symmetry_number
        self.mass            = mass
        self.formula         = formula
        if len(atoms) == 0:
            sys.exit("Trying to create a fragment without atoms")
        self.atoms = ""
        for a in atoms:
            self.atoms += " " + str(a)
        self.atoms = self.atoms.strip()
        
    def add_xml(self, xml):
        xml.set("identifier", self.identifier)
        xml.set("charge", str(self.charge))
        xml.set("multiplicity", str(self.multiplicity))
        xml.set("symmetry_number", str(self.symmetry_number))
        xml.set("mass", str(self.mass))
        xml.set("formula", self.formula)
        xml.text = self.atoms

class Molprop:
    '''
    Class to store a single molprop (molecule properties) from the
    Alexandria Chemistry Toolkit
    '''

    def __init__(self, molname):
        self.properties = {}
        self.properties["molname"] = molname
        self.bonds = []
        self.fragments = []
        self.experiments = []
        self.compounds = []
        
    def prop(self, propname):
        return self.properties[propname]

    def add_prop(self, propname, value):
        self.properties[propname] = str(value)

    def add_bond(self, atom1, atom2, bondorder):
        self.bonds.append([str(atom1), str(atom2), str(bondorder)])
    
    def add_experiment(self, exper):
        self.experiments.append(exper)
        
    def add_fragment(self, f:Fragment):
        self.fragments.append(f)

    def add_compound(self, comp): 
        self.compounds.append(comp)   
        
class Experiment:
    '''
    Class to hold Experiment data for a molecule
    '''
    def __init__(self, datasource, reference, program, method, basisset, conformation, jobtype, datafile, useForces):
        self.properties = {}
        self.properties["datasource"] = datasource
        self.properties["reference"]  = reference
        self.properties["program"]    = program
        self.properties["method"]     = method
        self.properties["basisset"]   = basisset
        self.properties["conformation"] = conformation
        self.properties["jobtype"]  = jobtype
        self.properties["datafile"] = datafile
        self.energies   = []
        self.dipoles    = []
        self.potential  = []
        self.quadrupole = []
        self.octupole   = []
        self.hexadecapole   = []
        self.polarisability = []
        self.frequencies = []
        self.intensities = []
        self.atoms = []
        self.tcmap = None
        self.useForce = useForces
        
    def add_prop(self, propname, value):
        self.properties[propname] = str(value)
        
    def add_energy(self, type, unit, temperature, phase, value):
        self.energies.append({ "type": type, "unit": unit, "temperature": str(temperature), "phase": phase, "average": str(value)})
    
    def add_dipole(self, type, unit, temperature, average, error, x, y, z): # phase,value
        self.dipoles.append({ "type": type, "unit": unit, "temperature": str(temperature), "average": str(average), "error": str(error), "x": str(x), "y": str(y), "z": str(z)})   #"phase": phase, "value": str(value)
    
    def add_potential(self, espid, potential_unit, coord_unit, x, y, z, V):
        self.potential.append({"espid": espid, "potential_unit": potential_unit, "coord_unit": coord_unit,"x": str(x), "y": str(y), "z": str(z), "V": str(V)})    

    def add_quadrupole(self, type, unit, temperature, xx, yy, zz, xy, xz, yz):
        self.quadrupole.append({"type":type, "unit":unit, "temperature": str(temperature), "xx":xx, "yy":yy, "zz":zz, "xy":xy, "xz":xz, "yz":yz})

    def add_octupole(self, type, unit, temperature, xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz):
        self.octupole.append({"type":type, "unit":unit, "temperature": str(temperature), "xxx":xxx, "xxy":xxy, "xxz":xxz, "xyy":xyy, "xyz":xyz, "xzz":xzz, "yyy":yyy, "yyz":yyz, "yzz":yzz, "zzz":zzz})

    def add_hexadecapole(self, type, unit, temperature, xxxx, xxxy, xxxz,
                         xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz,
                         yyyy, yyyz, yyzz, yzzz, zzzz):
        self.hexadecapole.append({"type":type, "unit":unit, "temperature": str(temperature), "xxxx":xxxx, "xxxy":xxxy, "xxxz":xxxz, "xxyy":xxyy, "xxyz":xxyz, "xxzz":xxzz, "xyyy":xyyy, "xyyz":xyyz, "xyzz":xyzz, "xzzz":xzzz, "yyyy":yyyy, "yyyz":yyyz, "yyzz":yyzz, "yzzz":yzzz, "zzzz":zzzz})

    def add_polarisability(self, type, unit, temperature,
                           average, error, xx, yy, zz, xy, xz, yz):    
        self.polarisability.append({"type":type, "unit": unit, "temperature": str(temperature), "average": str(average), "error":str(error), "xx":xx, "yy":yy, "zz":zz, "xy":xy, "xz":xz, "yz":yz})

    def add_atom(self, name, obtype, atomid, coord_unit, x, y, z,
                 force_unit, fx, fy, fz, qmap=None):
        newatom = { "name": name, "obtype": obtype, "atomid": str(atomid),
                    "coord_unit": coord_unit, 
                    "x": str(x), "y": str(y), "z": str(z),
                    "force_unit": force_unit,
                    "fx": str(fx), "fy": str(fy), "fz": str(fz) }
        if qmap:
            newatom["qmap"] = {}
            for q in qmap.keys():
                newatom["qmap"][q] = qmap[q]
        self.atoms.append(newatom)
        
    def extract_thermo(self, tcmap, atomname, ahof, verbose=False):
        if (not tcmap["Ezpe"] or not tcmap["Hcorr"] or not tcmap["Gcorr"] or
            not tcmap["E0"] or not tcmap["CV"] or not tcmap["Method"]):
            return
        eFactor   = UnitToConversionFactor("Hartree")
        Rgas      = 1.9872041*UnitToConversionFactor("kcal/mol")
        S0MT      = 0
        kilo      = 1000
        if tcmap["Temp"] > 0:
            S0MT += kilo*eFactor*(tcmap["Hcorr"]-tcmap["Gcorr"])/tcmap["Temp"]
        Srot    = -Rgas*math.log(float(tcmap["RotSymNum"]))
        if tcmap["RotSymNum"] > 1:
            Srot = 0
        S0MT += Srot
        DeltaSMT = S0MT
        dhofM0   = tcmap["E0"]*eFactor
        dhofMT   = dhofM0+(tcmap["Hcorr"]-tcmap["Ezpe"])*eFactor
        # Now fetch the atomic contributions to the heat of formation.
        # These depend on atomtype, temperature and method used as
        # well as the charge.
        charge = 0
        atomid = 0
        foundThermo = True
        for a in atomname:
            dhfx0, dhfxT, S0xT = ahof.get(a, tcmap["Temp"], charge)
            if dhfx0 and dhfxT and S0xT:
                dhofM0   += dhfx0
                dhofMT   += dhfxT
                DeltaSMT += S0xT
                atomid   += 1
            else:
                foundThermo = False
        # Final results, please check values and units!
        myT = float(tcmap["Temp"])
        self.add_energy("ZPE", "Hartree", 0, "gas", tcmap["Ezpe"])
        self.add_energy("CV", "J/mol K", 0, "gas", tcmap["CV"])
        self.add_energy("CP", "J/mol K", 0, "gas", Rgas+tcmap["CV"])
        if not foundThermo:
            if debug:
                print("Could not get atomization energies for %s at T = %g" % ( a, tcmap["Temp"]))
        else:
            self.add_energy("DeltaHform", "kJ/mol", myT, "gas", dhofMT)
            self.add_energy("DeltaGform", "kJ/mol", myT, "gas", dhofMT-myT*DeltaSMT/kilo)
            self.add_energy("DeltaHform", "kJ/mol", 0, "gas", dhofM0)
            self.add_energy("DeltaSform", "J/mol K", myT, "gas", DeltaSMT)
        Scomps = [ "Strans", "Srot", "Svib" ]
        if len(tcmap["Scomponent"]) == len(Scomps):
            for i in range(len(Scomps)):
                self.add_energy(Scomps[i], "J/mol K", myT, "gas", tcmap["Scomponent"][i])


    
def test_molprops():
    molprops = Molprops()
    
    mp1 = Molprop("water", 0, 0, 1, 1, "compound_1", "compound_2")
    mp1.add_prop("formula", "H2O")
    mp1.add_prop("mass", "18")
    mp1.add_bond(1, 2, 1)
    mp1.add_bond(1, 3, 1)
    
    exper = Experiment("Theory", "Author2021a")
    exper.add_energy("HF", "Hartree", 0.0, "G", -76.2)
    exper.add_energy("MP2", "Hartree", 0.0, "G", -76.3)
    exper.add_atom("O", "ow", 1, 0.0, 0.0, 0.117, "Angstrom")
    exper.add_atom("H", "hw", 2, 0.0, 0.7634, -0.4681, "Angstrom")
    exper.add_atom("H", "hw", 3, 0.0, -0.7634, -0.4681, "Angstrom")
    
    mp1.add_experiment(exper)
    molprops.add_molecule(mp1)
    molprops.write("test.xml")
    
#test_molprops()
