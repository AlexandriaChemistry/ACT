#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, shutil, argparse

from pathlib     import Path

def get_environ(myvar: str):
  if myvar in os.environ:
    return os.environ[myvar]
  else:
    print("Please set the %s environment variable" % myvar)
    exit(1)
    
AlexandriaLib = get_environ("AlexandriaLib")

def parseArguments():
  parser = argparse.ArgumentParser()
  maxpot = 100
  parser.add_argument("-maxpot", "--maxpot",       help="Percent of the total number of ESP, default "+str(maxpot), type=int, default=maxpot)
  sel = "../SELECTIONS/alcohol.dat"
  parser.add_argument("-sel", "--selection",       help="Selection file, default "+sel, type=str, default=sel)
  output = "alcohol.dat"
  parser.add_argument("-o", "--output",            help="Output file, default "+output, type=str, default=output)
  gentop = "../ACS-pg.dat"
  parser.add_argument("-d", "--gentop",            help="Force field input file, default "+gentop, type=str, default=gentop)
  args = parser.parse_args()
  return args

def XML_directory():
  return ALEXANDRIA + "XML"

def name_of_data_log_file(molname: str, lot: str):
  logfile = ( "%s/compounds/%s/%s-%s.log.gz" % ( AlexandriaLib, molname, molname, lot ))
  if os.path.exists(logfile):
    datafile = ( "%s-3.dat" % molname )
    return datafile, logfile
  else:
    logfile = ( "%s/compounds/%s/%s-%s-oep.log.gz" % ( AlexandriaLib, molname, molname, lot ))
    if os.path.exists(logfile):
      datafile = ( "%s-3.dat" % molname )
      return datafile, logfile

  return None, None

def get_alexandria():
  for name in ["alexandria", "alexandria_d", "alexandria_mpi", "alexandria_mpi_d"]:
    alexandria = shutil.which(name)
    if alexandria:
      return alexandria
  return None
  
def write_makefile(selection, gentop, maxpot, method, basisset, qm2molprop, merge_mp, output):
  cwd = os.getcwd()
  makefile = open("Makefile", "w")
  makefile.write("all: %s\n\n" % output)
  datafiles = []
  lot = ( "%s-%s" % ( method, basisset ))
  with open(selection, "r") as inf:
    for line in inf:
      words = line.split("|")
      if len(words) != 2:
        continue
      molname = words[0]
      datafile, logfile = name_of_data_log_file(molname, lot)
      if not datafile or not logfile:
        print("Missing logfile for target %s" % ( datafile ))
        continue
      datafiles.append(datafile)
      makefile.write("# Rule for %s.\n" % datafile)
      makefile.write("%s: %s\n" % (datafile, logfile))
      makefile.write("\t\t%s -g03 %s -o %s -d %s -maxpot %d -molnm %s -iupac %s -basis %s -jobtype Opt\n\n" %
                     (qm2molprop, logfile, datafile, gentop, maxpot, molname, molname, basisset))
  files = ' '.join(map(str, datafiles))
  makefile.write("\n\n")
  makefile.write("\n%s: %s\n" % ( output, files ) )
  makefile.write("\t\t%s -di %s -o %s -maxwarn 0 -f %s\n" % (merge_mp, gentop, output, files))
  makefile.write("\nclean:\n")
  makefile.write("\trm -f *.dat *.sdf *.debug \\#* Makefile\n")
  makefile.write("\techo 'run ./extract_xml.py to recreate the Makefile'\n")
  makefile.close()      
  os.chdir(cwd)

if __name__ == "__main__":

  args        = parseArguments()
  alexandria  = get_alexandria()
  if not alexandria:
    print("Can not find the alexandria executable. Stopping script.")
    exit(1)
  selection     = args.selection
  gentop        = args.gentop
  qm2molprop    = alexandria + " qm2molprop"
  merge_mp      = alexandria + " merge_mp"

  write_makefile(args.selection, gentop, min(100, max(20, args.maxpot)), 
                 "B3LYP", "aug-cc-pVTZ", qm2molprop, merge_mp, args.output)

  print("Makefile generated. Now please run make -j 8 %s" % args.output)
