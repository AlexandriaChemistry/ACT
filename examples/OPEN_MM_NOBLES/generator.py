#!/usr/bin/env python3
import random
import numpy as np
import argparse, sys
import time
coordlist = []
coords = [0, 0, 0]
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
out.write(f"CRYST1   100.20   100.20   100.20  90.00  90.00  90.00 P 1           1\n")
out.write(f"MODEL        1\n")

def generate():
    Xcoord = int(random.randint(5, 49)*2)
    Ycoord = int(random.randint(5, 49)*2)
    Zcoord = int(random.randint(5, 49)*2)
    coords = [Xcoord, Ycoord, Zcoord]
    return coords


def checker(loop):
    global coords
    coords = generate()
    if coords not in coordlist:
        print(f"unique{coords}")
        coordlist.append(coords)
        return(coords)
        loop = 0
    else:
        print("not unique")
        loop += 1
        print(loop)
        if loop > 100:
            print("Not found in 100 iterations")
            time.sleep(0.5)
            sys.exit()
            time.sleep(0.5)
        checker(loop)


#generate
#print(generate())


for i in range(atoms):
    print(i)

    checker(0)
    print(f" zde {coords}")
#    PICKX = random.choice(Xlist) + 10
#    PICKZ = random.choice(Xlist) + 10
#    PICKY = random.choice(Xlist) + 10
    if i < 9:
        out.write(f"ATOM      {i+1} {typ}   {typ}      {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (coords[0], coords[1], coords[2]))
    if i == 9: 
        out.write(f"ATOM     {i+1} {typ}   {typ}     {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (coords[0], coords[1], coords[2]))
    if i < 99 and i > 9:
        out.write(f"ATOM     {i+1} {typ}   {typ}     {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (coords[0], coords[1], coords[2]))
    if i > 98 and i < 999: 
        out.write(f"ATOM    {i+1} {typ}   {typ}    {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (coords[0], coords[1], coords[2]))
    if i > 998:
        out.write(f"ATOM   {i+1} {typ}   {typ}   {i+1}      %0.3f  %0.3f  %0.3f  1.00  0.00          {typ}\n"% (coords[0], coords[1], coords[2]))
















out.write(f"TER\n")
out.write(f"ENDMDL\n")
