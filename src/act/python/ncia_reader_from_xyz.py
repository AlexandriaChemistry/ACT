#!/usr/bin/env python3
import os
from weakref import ref
import yaml
import json
import numpy as np
#from xyz2molprop import *
#import re
from itertools import islice
import subprocess
import time
import argparse
listinchi = []
listmons = []
if os.path.exists("sel_s66.dat"):
  os.remove("sel_s66.dat")
if os.path.exists("sel_s66_mons.dat"):
  os.remove("sel_s66_mons.dat")
if os.path.exists("inchi.dat"):
  os.remove("inchi.dat")
selfile = open("sel_s66.dat", "a+")
MIN = "No"
def parseArguments():
    desc = '''This script will read xyz files in a NCIA directory and pass it on the conversion script. Provided that they have a header with information like molecular segments, etc...
    the level of theory is CCSD(T)/CBS'''
    parser = argparse.ArgumentParser(description=desc)
    defdir = "S66x8_geom"
    parser.add_argument("-dir", "--directory", help="The directory with xyz files"+defdir, type=str, default=defdir)
    deffor = "S66"
    parser.add_argument("-frm", "--format", help="Either NCIA, or S66."+deffor, type=str, default=deffor)

#    parser.add_argument("-ref", "--refe", help="reference", type=str, default="")
    return parser.parse_args()
args  = parseArguments()
listofxyz = list(os.listdir(args.directory))
for line in listofxyz:
    
        nameinc = line.strip(".xyz\n")
        name = nameinc.strip()
        print(f"Name: {name}")
        selfile.write(f"{name}|Train\n")
        Aname = name + "_A"
        Bname = name + "_B"
        if args.format == "S66":
            if "." not in nameinc:
                print("This is probably version with only equilibrium geometries.")
                MIN = "Yes"
                scaling = 1.0
            else:
                scaling = nameinc.partition("_")[2].partition("_")[2]
                if "_" in scaling:
                    scaling = scaling.partition("_")[2]
                print(f"Scaling {scaling}" )
                if float(scaling) == 1.0:
                    print("EQVILIBRIUM")
                    MIN = "Yes"
                else: 
                    MIN = "No"
        else: 
            scaling = 0.0

###list of monomers
        if Aname not in listmons:
            listmons.append(Aname)
        if Bname not in listmons:
            listmons.append(Bname)
        
        geo = open(f"{args.directory}/{line}", "r")
        for header in islice(geo, 1, 2):
            #count += 1
            benchmark = 0.0
            unit = "empty"

            print(header.strip())
            if "-" in header.partition("selection_a=")[2].partition("selection_b")[0]:
             #   SEL_A_1 = header.partition("selection_a=")[2].partition("selection_b")[0].split("-")[0]
                SEL_A_1 = header.partition("selection_a=")[2].split("-")[0]
            else:
            #    SEL_A_1 = header.partition("selection_a=")[2].partition("selection_b")[0]
                SEL_A_1 = header.partition("selection_a")[2].partition("selection_b")[0].split("=")[1]
                print("A1, no - CCCCCCCCCCCCCCCCCCCCCCCCC")
            if "-" in header.partition("selection_a=")[2].partition("selection_b")[0]:
                SEL_A_2 = header.partition("selection_a=")[2].partition("selection_b")[0].split("-")[1]
            else:
                SEL_A_2 = header.partition("selection_a")[2].partition("selection_b")[0].split("=")[1]
                print("A1, no - DDDDDDDDDDDDDDDDDDDDDDDDDD")
         #   if "-" in  
           # if "-" in header.partition("selection_b=")[2]
            if args.format == "S66": ## no monoatomics in S66
                SEL_B_1 = header.partition("selection_b=")[2].split("-")[0]
                SEL_B_2 = header.partition("selection_b=")[2].split("-")[1]
            if args.format == "NCIA":
                if "-" in header.partition("selection_b=")[2].partition("scaling")[0]:
                    SEL_B_1 = header.partition("selection_b=")[2].partition("scaling")[0].split("-")[0]
                else: 
                    SEL_B_1 = header.partition("selection_b")[2].partition("scaling")[0].split("=")[1]
                    print("B1, no - AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")

            CHAR_A = header.partition("charge_a=")[2].partition("charge_b")[0]
            CHAR_B= header.partition("charge_b=")[2].partition("selection_a")[0]
            print(f"charge B: {CHAR_B}")
            print(f"charge A: {CHAR_A}")
            if args.format == "NCIA":
                benchmark = float(header.partition("benchmark_Eint=")[2].partition("benchmark_unit")[0])
                if "-" in header.partition("selection_b=")[2].partition("scaling")[0]:
                    SEL_B_2 = header.partition("selection_b=")[2].partition("scaling")[0].split("-")[1]
                else:
                    SEL_B_2 = header.partition("selection_b")[2].partition("scaling")[0].split("=")[1]
                    print("B2, no - XXXXXXXXXXXXXXXXXXXXXXXXXX")
                unit = header.partition("benchmark_unit=")[2].partition("group")[0]
                scaling = float(header.partition("scaling=")[2].partition("benchmark_Eint")[0])
            if args.format == "S66":
                benchmark = float(header.partition("reference_value=")[2].partition("charge")[0])
                SEL_B_2 = header.partition("selection_b=")[2].split("-")[1]
                unit = header.partition("reference_unit=")[2].partition("reference_value")[0]
            if benchmark == 0.0:
                print("Warning, benchamrk value 0.0, probably not recognized **************************************")
            if unit == "empty":
                print("Warning, benchmark unit not recognized ******************************************")
            print(f"Unit: {unit}")
            if "kcal" in unit and "mol" in unit:
                print("Converting kcal/mol to kj/mol")
                benchmark = benchmark * 4.183995
                print(f"this is benchmark {benchmark} kj/mol")
            else:
                print(f"this is benchmark {benchmark} {unit}")
            if args.format == "NCIA":
                if scaling == 1.00 or scaling == 1.0:
                    MIN = "Yes"
                    print("Equilibrium")
                else:
                    MIN = "No"
                    print("Not equilibrium")
            print(f"Scaling: {scaling}")

            print(f"SEL_B_1: {SEL_B_1}")
            print(f"SEL_B_2: {SEL_B_2}")
            print(f"SEL_A_1: {SEL_A_1}")
            print(f"SEL_A_2: {SEL_A_2}")
            if os.path.exists("SEL_A.xyz"):
                os.remove("SEL_A.xyz")

            if os.path.exists("SEL_B.xyz"):
                os.remove("SEL_B.xyz")
            sub = int(SEL_B_2) - int(SEL_A_2)
            print(f"subtract {sub}")
            SEL_A = open("SEL_A.xyz", "a+")
            SEL_A.write(f"{SEL_A_2}\n") 
            SEL_A.write(f"\n") 
            SEL_B = open("SEL_B.xyz", "a+")
            SEL_B.write(f"{sub}\n") 
            SEL_B.write(f"\n") 
###temporary selection files
            for lines in range(int(SEL_B_1) - 1):
                ln = geo.readline()
                #print(f"{lines}th {ln}")
                SEL_A.write(f"{ln}") 
            for lines in range(int(SEL_A_2), int(SEL_B_2)):
                ln = geo.readline()
                #print(f"{lines}th {ln}")
                SEL_B.write(f"{ln}") 
            SEL_A.close()
            SEL_B.close()
            ###INCHI PART
            subprocess.run("obabel -ixyz SEL_A.xyz -oinchi -OSELAinchi".format().split())
            subprocess.run("obabel -ixyz SEL_B.xyz -oinchi -OSELBinchi".format().split())
            Binchi = open(f"SELBinchi", "r")
            Binchi_read = str(Binchi.read()).split("=")[1].strip()
            Ainchi = open(f"SELAinchi", "r")
            Ainchi_read = str(Ainchi.read()).split("=")[1].strip()
            print(Binchi_read)
            print(Ainchi_read)
            if Ainchi_read not in listinchi:
                listinchi.append(f"{Ainchi_read}")
            if Binchi_read not in listinchi:
                listinchi.append(Binchi_read)
###
        ##    try:
            subprocess.run("./ncia2molprop.py -i {}/{}.xyz -o {}.xml -n {} -int {} -charge_a {} -charge_b {} -sel_a_end {} -sel_b_end {} -MIN {} -ID {} -FRAGA SEL_A.xyz -FRAGB SEL_B.xyz -ref {}".format(args.directory,name, name,name,benchmark,CHAR_A,CHAR_B,SEL_A_2,SEL_B_2, MIN, name, args.directory).split())
          ##  except:
            #    print("error occoured in the molprop convertor")

#            time.sleep(1)
selmons = open("sel_s66_mons.dat", "a+")
inchies = open("inchi.dat", "a+")
for inch in listinchi:
    print([inch])
    inchies.write(f"{inch}|Train\n")
for mon in listmons:
    selmons.write(f"{mon}|Train\n")









print("Taking kcal/mol as reference unit")
