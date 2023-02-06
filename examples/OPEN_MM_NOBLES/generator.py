#!/usr/bin/env python3
import random
import numpy as np
import argparse, sys

def parseArguments():
    desc = '''This script will generate a box filled with n atoms (recommended for noble gases)'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-num", "--num", help="how many atoms are to be generated", type=int)
    outf = "out.pdb"
    parser.add_argument("-out", "--outfile", help="Output pdb file for writing, default "+outf, type=str, default=outf)
    parser.add_argument("-type", "--type", help="Atom type ", type=str, default="Ar")
    return parser.parse_args()
args  = parseArguments()

atoms = args.num
typ = args.type

out = open(args.outfile, "w")
out.write(f"REMARK    THIS IS A SIMULATION BOX\n")
out.write(f"CRYST1   28.200   28.200   28.200  90.00  90.00  90.00 P 1           1\n")
out.write(f"MODEL        1\n")

Xlist = list(range(0, 70))
Ylist = list(range(0, 70))
Zlist = list(range(0, 70))
for i in range(atoms):
    print(i)
    PICKX = random.choice(Xlist) + 10
    PICKZ = random.choice(Xlist) + 10
    PICKY = random.choice(Xlist) + 10
    if i < 9:
       out.write(f"ATOM      {i+1} {typ}   {typ}      {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (PICKX, PICKY, PICKZ))
    if i == 9: 
     out.write(f"ATOM     {i+1} {typ}   {typ}     {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (PICKX, PICKY, PICKZ))
    if i < 99 and i > 9:
       out.write(f"ATOM     {i+1} {typ}   {typ}     {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (PICKX, PICKY, PICKZ))
    if i > 98 and i < 999: 
        out.write(f"ATOM    {i+1} {typ}   {typ}    {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (PICKX, PICKY, PICKZ))
    if i > 998:
        out.write(f"ATOM   {i+1} {typ}   {typ}   {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (PICKX, PICKY, PICKZ))

out.write(f"TER\n")
out.write(f"ENDMDL\n")
