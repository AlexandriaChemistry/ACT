#!/usr/bin/env python3

import os, sys, subprocess, argparse, csv
from enum    import Enum
#from dbutils import *
#from legacy_models import *
    
BohrInverse2nmInverse = 18.897268

delim = "|"

def get_csv_rows(csv_file, minimum_nr_columns, enc='utf-8'):
    inputfile = open(csv_file, "r", encoding=enc)
    csv.register_dialect('pipes', delimiter=delim)

    rows = []
    try:
        reader = csv.reader(inputfile, dialect='pipes')

        for row in reader:
            if (len(row) >= minimum_nr_columns) and (row[0].find("#") < 0):
                rows.append(row)
    finally:
        inputfile.close()
    print("Read %d rows from %s" % ( len(rows), csv_file ))
    return rows

def library_filename(filename):
  if not os.path.exists(filename):
    actdata = "ACTDATA"
    if actdata in os.environ:
      filename = ("%s/forcefield/%s" % ( os.environ[actdata], filename ))
    else:
      sys.exit("Cannot find %s" % filename)
  return filename

def string2bool(string):
  return string == "True" or string == "true" or string == "yes" or string == "Yes"

def param_header(out, identifier):
  out.write("    <parameterlist identifier=\"%s\">\n" % ( identifier ) )

def param_footer(out):
  out.write("    </parameterlist>\n")
  
def format_params(function, types, units, minval, maxval, sigma,
                  mutability, ntrain, nonnegative):
  try:
    value = str((float(minval)+float(maxval))/2)
  except ValueError:
    return ""
  if not sigma:
    sigma = ""
  nnstr = "no"
  if nonnegative:
    nnstr = "yes"
  return ("      <parameter type=\"%s\" unit=\"%s\" value=\"%s\" uncertainty=\"%s\" minimum=\"%s\" maximum=\"%s\" ntrain=\"%d\" mutability=\"%s\" nonnegative=\"%s\"/>" % ( types, units, value, sigma, minval, maxval, ntrain, mutability, nnstr))

class AtomTypes:
  """Class to handle atomtypes for force field files"""
  def __init__(self):
    self.atp = {}
    
  def read4(self, filename, aprops, polarizable):
    rows   = get_csv_rows(library_filename(filename), 35)
    header = rows[0]
    # Expecting these elements in the file
    elements = [ "atomtype", "element", "q_min", "q_max", "charge_mutability", "bondtype", "poltype", "zetatype", "acmtype", "epsilon_min", "epsilon_max", "epsilon_mutability", "gamma_min", "gamma_max", "gamma_mutability", "sigma_min", "sigma_max", "sigma_mutability", "radius", "row", "alpha_min", "alpha_max", "alpha_mutability", "zeta_min", "zeta_max", "zeta_mutability", "jaa_min", "jaa_max", "jaa_mutability", "chi_min", "chi_max", "chi_mutability", "comment", "ref_enthalpy", "reference" ]
    # Simple check for completeness and correctness
    if len(header) != len(elements):
      sys.exit("Found %d elements on the header line in %s, expected %d" % ( len(header), filename, len(elements)))
    for k in range(len(header)):
      if header[k] != elements[k]:
        sys.exit("Found header[%d] = %s, expected %s" % ( k, header[k], elements[k] ))
    
    # Now fill the atp dictionnary
    n = int(len(self.atp.keys()))
    for rr in range(1,len(rows)):
      row = rows[rr]
      key = row[0]
      if not polarizable and (key[:-2] == "_s" or row[1] == 'X'):
        sys.exit("Atomtypes file contains shells. Please select another one.")
      if key in self.atp:
        sys.exit("Duplicate atomtype %s in %s" % ( key, filename ))
      self.atp[key] = {}
      self.atp[key]["index"] = n
      for kk in range(1,len(row)):
        if len(row[kk]) > 0:
          self.atp[key][elements[kk]] = row[kk]
      elem = self.atp[key]["element"]
      self.atp[key]["mass"]       = aprops.mass(elem)
      self.atp[key]["atomnumber"] = aprops.atomnumber(elem)
      self.atp[key]["vdwtype"]    = key
      
  def write_parm(self, out, mytype, myunit, mymut, minval, maxval, nonnegative):
    nnstr = "no"
    if nonnegative:
      nnstr = "yes"
    value = str((float(minval)+float(maxval))/2)
    out.write("      <parameter type=\"%s\" unit=\"%s\" value=\"%s\" uncertainty=\"0\" minimum=\"%s\" maximum=\"%s\" ntrain=\"0\" mutability=\"%s\" nonnegative=\"%s\"/>\n" % ( mytype, myunit, value, minval, maxval, mymut, nnstr ))
    
  def write_xml(self, myff, out, nexcl):
    out.write("  <version checksum=\"\" timestamp=\"\"/>\n")
    out.write("  <particletypes nexclusions=\"%s\" epsilonr=\"1\">\n" % (nexcl))
    shell_written = {}
    for key in sorted(list(self.atp.keys()), key=lambda y: (self.atp[y]["index"])):
      mymass = self.atp[key]["mass"]
      if not "comment" in self.atp[key]:
        sys.exit("No comment to describe atomtype %s" % key)
      mytype = "Atom"
      if self.atp[key]["element"] == "X":
        mytype = "Shell"
      out.write("    <particletype identifier=\"%s\" type=\"%s\" description=\"%s\">\n" % (key, mytype, self.atp[key]["comment"]))
      for opt in [ "acmtype", "atomnumber", "bondtype", "element", "poltype", "row", "vdwtype", "zetatype" ]:
        if (opt != "poltype" or myff.polarizable):
          if opt in [ "row", "atomnumber" ] or (opt in self.atp[key] and len (self.atp[key][opt]) > 0):
            out.write("      <option key=\"%s\" value=\"%s\"/>\n" % ( opt, self.atp[key][opt]))
      self.write_parm(out, "charge", "e", self.atp[key]["charge_mutability"], self.atp[key]["q_min"], 
                      self.atp[key]["q_max"], False)
      self.write_parm(out, "mass", "Da", "Fixed", mymass, mymass, True)
      self.write_parm(out, "ref_enthalpy", "kJ/mol", "Fixed", self.atp[key]["ref_enthalpy"],
                      self.atp[key]["ref_enthalpy"], True)
      out.write("    </particletype>\n")
    out.write("  </particletypes>\n")

  def keys(self):
    return sorted(list(self.atp.keys()), key=lambda y: (self.atp[y]["index"]))
    
  def subtype(self, key, subtype):
    if subtype in self.atp[key]:
      return self.atp[key][subtype]
    return None
  
  def set_subtype(self, key, subtype, value):
    self.atp[key][subtype] = value

  def delete(self, key):
    del self.atp[key]
    
  def set_vdwparam(self, atomtype, sigma, epsilon, gamma):
    if not atomtype in self.atp:
      sys.exit("No such atomtype %s" % atomtype)
    self.atp[atomtype]["sigma"]   = sigma
    self.atp[atomtype]["gamma"]   = gamma
    self.atp[atomtype]["epsilon"] = epsilon
    
  def set_reference(self, atomtype, reference):
    if not atomtype in self.atp:
      sys.exit("Unknown atomtype %s when adding reference %s" % ( atomtype, reference ))
    self.atp[atomtype]["reference"] = reference
    
  def print_vdwparams(self, outf):
    vdwparamlist = []
    units = { "epsilon": "kJ/mol", "gamma": "", "sigma": "nm" }
    for key in sorted(list(self.atp.keys()), key=lambda y: (self.atp[y]["index"])):
      param_header(outf, key)
      for u in units.keys():
        outf.write("%s\n" % format_params(key, u, units[u], 
                                          self.atp[key][u+"_min"], self.atp[key][u+"_max"],
                                          None, self.atp[key][u+"_mutability"], 0, True))
      param_footer(outf)
      
class ForceField:
  """Class to hold generic FF information"""
  def __init__(self, qdist, filename="ffnames.csv"):
    rows  = get_csv_rows(library_filename(filename), 7)
    found = False
    for row in rows:
      if row[0] == qdist:
        self.qdist = qdist
        self.chargetype = row[1]
        self.chargealgorithm = row[2]
        self.polarizable = string2bool(row[3])
        self.vanderwaals = row[4]
        self.useIons = string2bool(row[5])
        self.reference = row[6]
        found = True
        print("chargetype  = {}".format(self.chargetype))
        print("vanderwaals = {}".format(self.vanderwaals))
        print("polarizable = {}".format(self.polarizable))
    if not found:
      print("Unsupported force field %s" % qdist)
      exit(1)

class AtomProps:
  """Basic atom properties"""
  def __init__(self, filename="atomprops.csv"):
    self.atompropsdb = {}
    for row in get_csv_rows(library_filename(filename), 4):
      if len(row) == 4:
        self.atompropsdb[row[0]] = { "atomnumber": int(row[2]), "mass": float(row[3]) }
    
  def mass(self, elem):
    return self.atompropsdb[elem]["mass"]
  
  def atomnumber(self, elem):
    return self.atompropsdb[elem]["atomnumber"]

def write_syms(out):
  syms = [[ "C", "H", "3" ],[ "N", "H", "3" ],[ "N", "H", "2" ],
          [ "N", "O", "2" ],[ "O", "H", "2" ],[ "C", "F", "3" ],
          [ "C", "Cl", "3" ],[ "C", "Br", "3" ],[ "C", "I", "3" ],
          [ "Si", "H", "3" ]]
  out.write("  <symmetric_charges>\n")
  for sym in range(len(syms)):
     out.write("    <sym_charge central=\"%s\" attached=\"%s\" numattach=\"%s\"/>\n" % (syms[sym][0], syms[sym][1], syms[sym][2]))
  out.write("  </symmetric_charges>\n")

def write_bonds(out):
  with open(library_filename("gt_bonds.xml"), "r") as inf:
    for line in inf.readlines():
      out.write(line)

def write_interactions(out, fftype, function, canswap, terminate, options=None):
  if canswap:
    swap = "true"
  else:
    swap = "false"
  out.write("  <interaction type=\"%s\" function=\"%s\" canswap=\"%s\">\n" %
            ( fftype, function, swap ))
  if options:
    for opt in options:
      if len(opt) == 2:
        out.write("    <option key=\"%s\" value=\"%s\"/>\n" % ( opt[0], opt[1] ))
      else:
        print("Options should contain a type and a value, ignoring {}".format(opt))
  if terminate:
    out.write("  </interaction>\n")

def print_one(out, mytype: str, unit: str, minval: str, maxval: str, mutability: str, nonnegative: bool):
  value = str((float(minval)+float(maxval))/2)
  nneg = "no"
  if nonnegative:
    nneg = "yes"
  out.write("      <parameter type=\"%s\" unit=\"%s\" value=\"%s\" minimum=\"%s\" maximum=\"%s\" uncertainty=\"0\" mutability=\"%s\" ntrain=\"1\" nonnegative=\"%s\"/>\n" %
            ( mytype, unit, value, minval, maxval, mutability, nneg ) )

class BondVsites:
  """Class to handle bond shells for force field files"""
  def __init__(self):
    self.bs = []
    
  def read(self, filenm):
    self.bs = get_csv_rows(filenm, 4)
    
  def print_vsites(self, outf):
    outf.write("  <interaction type=\"VSITE2\" function=\"vsite2\" canswap=\"false\">\n")
    for bs in self.bs:
      myid = ("%s~%s~%s" % ( bs[0], bs[1], bs[2]) )
      param_header(outf, myid)
      print_one(outf, "v2_a", "", bs[3], bs[3], "Fixed", True)
      param_footer(outf)
      if bs[0] != bs[1]:
        myid = ("%s~%s~%s" % ( bs[1], bs[0], bs[2]) )
        param_header(outf, myid)
        aa = str(1.0-float(bs[3]))
        print_one(outf, "v2_a", "", aa, aa, "Fixed", True)
        param_footer(outf)
    outf.write("  </interaction>\n")

def write_charge_dist(out, myatp):
  out.write("  <interaction type=\"CHARGEDISTRIBUTION\" function=\"\" canswap=\"false\">\n")
  out.write("    <option key=\"chargetype\" value=\"%s\"/>\n" % (myff.chargetype))
  out.write("    <option key=\"reference\" value=\"%s\"/>\n" % (myff.reference))
  mycd = {}
  for elem in sorted(myatp.keys()):
    zt = "zetatype"
    if not zt in myatp.atp[elem]:
      sys.exit("No %s for %s" % ( zt, elem ))
    zetatype = myatp.atp[elem][zt]
    if not zetatype in mycd:
      mycd[zetatype] = True
      # Check mutability
      mutab   = myatp.atp[elem]["zeta_mutability"]
      param_header(out, myatp.atp[elem][zt])
      print_one(out, "zeta","1/nm", myatp.atp[elem]["zeta_min"], myatp.atp[elem]["zeta_max"], mutab, True)
      param_footer(out)
  out.write("  </interaction>\n")

def write_eem(out, myatp, reference):
  out.write("  <interaction type=\"ELECTRONEGATIVITYEQUALIZATION\" function=\"\" canswap=\"false\">\n")
  out.write("    <option key=\"reference\" value=\"%s\"/>\n" % (reference))
  myeem = {}
  for elem in sorted(myatp.keys()):
    if len(elem) > 0:
      if elem in myeem:
        continue
      myeem[elem] = True
      if "chi_min" in myatp.atp[elem] and "chi_max" in myatp.atp[elem] and "jaa_min" in myatp.atp[elem] and "jaa_max" in myatp.atp[elem] and "acmtype" in myatp.atp[elem]:
        if len(myatp.atp[elem]["acmtype"]) > 0:
          out.write("    <parameterlist identifier=\"%s\">\n" % myatp.atp[elem]["acmtype"])
          print_one(out, "chi", "eV",    myatp.atp[elem]["chi_min"], myatp.atp[elem]["chi_max"],
                    myatp.atp[elem]["chi_mutability"], True)
          print_one(out, "jaa", "eV/e",  myatp.atp[elem]["jaa_min"], myatp.atp[elem]["jaa_max"],
                    myatp.atp[elem]["jaa_mutability"], True)
          out.write("    </parameterlist>\n")
  out.write("  </interaction>\n")

def write_polarization(out, myff, myatp):
  write_interactions(out, "POLARIZATION", "Polarization", False, False)
  if myff.polarizable:
    out.write("    <option key=\"reference\" value=\"Molina2011a,Ghahremanpour2018b\"/>\n")
    for key in myatp.keys():
      alpha_min = myatp.subtype(key, "alpha_min")
      alpha_max = myatp.subtype(key, "alpha_max")
      if alpha_min and alpha_max:
        mutab  = myatp.subtype(key, "alpha_mutability")
        param_header(out, key)
        out.write("%s\n" % format_params("Polarization", "alpha", 
                                         "Angstrom3", alpha_min, alpha_max, "0", mutab, 
                                         0, True))
        param_footer(out)
  out.write("  </interaction>\n")

def print_gentop(output, myff, myatp, bs,
                 combRule, nexcl, pdihs, 
                 vsite, slaterMax, reference):
  out = open(output, "w")
  out.write("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n")
  out.write("<!DOCTYPE gentop.dtd PUBLIC \"gentop.dtd\" \"gentop.dtd\">\n")
  out.write("<gentop>\n")

  myatp.write_xml(myff, out, nexcl)
  
  write_polarization(out, myff, myatp)
  
  if vsite:
    out.write("  <gt_vsites angle_unit=\"degree\" length_unit=\"pm\">\n")
    for key in sorted(vsite.keys()):
      out.write("    <gt_vsite atype=\"%s\" vtype=\"%s\" number=\"%s\" distance=\"%0.0f\" angle=\"%s\" ncontrolatoms=\"%s\"/>\n" %
             (key, vsite[key]['vtype'],
             vsite[key]['number'], 1000*float(vsite[key]['distance']),
             vsite[key]['angle'], vsite[key]['ncontrolatoms']))
    out.write("  </gt_vsites>\n")
  write_interactions(out, "VANDERWAALS", myff.vanderwaals,
                    True, False,
                     [ [ "combination_rule", combRule ] ] )
  myatp.print_vdwparams(out)
  out.write("  </interaction>\n")
  if True:
    write_bonds(out)
  else:
    write_interactions(out, "BONDS", "MORSE", True, True)
    write_interactions(out, "ANGLES", "UREY_BRADLEY", True, True)
    write_interactions(out, "LINEAR_ANGLES", "LINEAR_ANGLES", False, True)
    write_interactions(out, "PROPER_DIHEDRALS", pdihs, True, True)
    write_interactions(out, "IMPROPER_DIHEDRALS", "IDIHS", True, True)

  write_syms(out)  
  write_charge_dist(out, myatp)
  bs.print_vsites(out)
  write_eem(out, myatp, reference)
  out.write("</gentop>\n")
  out.close()

  print("Generated %s\n" % output)

if __name__ == '__main__':

    def parseArguments():
        parser = argparse.ArgumentParser()
        parser.add_argument("-o",        "--output",     help="Output file name, will be generated if empty",     type=str, default="")
        parser.add_argument("-comb",     "--combRule",   help="combination rule:  arithmetic, geometric, kong.",  type=str,   default = "Geometric")
        parser.add_argument("-g",        "--gamma",      help="gamma value for wbk potential.",                   type=float, default =  13.0)
        parser.add_argument("-nexcl",    "--nexcl",      help="number of exclusions.",                            type=int,   default =  0)
        scaleFactor = 0.8
        parser.add_argument("-sf",       "--scaleFactor",help="Scaling factor to set minimum and maximum allowed values of parameters, according to min = factor*value, max = value/factor. Number should be > 0 and <= 1", type=float, default=scaleFactor)
        parser.add_argument("-dih",      "--dihedral",   help="Which GROMACS function to use for dihedrals, PDIHS or FOURDIHS", type=str, default="FOURDIHS")
        parser.add_argument("-ff",    "--force_field",   help="The charge distribution scheme for which eemprops will be written", type=str, default="ACM-g")
        parser.add_argument("-smax",     "--slater_max", type=int, default=3)
        parser.add_argument("-bv",       "--bondvsites", help="Added vsites to selected bonds and for polarizable models a shell to the vsite", type=str, default=None)
        defaprops = "atomtypes.csv"
        parser.add_argument("-atypes", "--atomtypes", help=("Definition to use when generating a force field file (default %s)" % defaprops), type=str, default=defaprops)
        args = parser.parse_args()
        return args

    arg    = parseArguments()
    myff   = ForceField(arg.force_field)
    aprops = AtomProps()
    myatp  = AtomTypes()
    myatp.read4(arg.atomtypes, aprops, myff.polarizable)
      
    # Polarizability check
    bs = BondVsites()
    if arg.bondvsites:
      bs.read(arg.bondvsites)
    output = arg.output
    if len(output) == 0:
      output = arg.force_field
      if arg.bondvsites:
        output += "-bv"
      output += ".xml"

    print_gentop(output, myff, myatp, bs,
                 arg.combRule, arg.nexcl, arg.dihedral,
                 None, arg.slater_max, 
                 "Ghahremanpour2022a")
