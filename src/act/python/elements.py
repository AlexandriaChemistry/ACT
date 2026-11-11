#
# This source file is part of the Alexandria Chemistry Toolkit.
#
# Copyright (C) 2026
#
# Developers:
#             Mohammad Mehdi Ghahremanpour,
#             Paul J. van Maaren,
#             David van der Spoel (Project leader)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.
#

import os, sys
from actutils     import *
from get_csv_rows import *

atomprops = {}
def get_atomprops():
    adata = act_library_filename("atomprops.csv")
    for row in get_csv_rows(adata, 7, delim=","):
        try:
            atomprops[row[0]] = { "symbol": row[1], "name": row[2], "atomicnumber": int(row[3]), "mass": float(row[4]),
                                  "charge": int(row[5]), "mult": int(row[6]) }
        except ValueError:
            print("Could not understand line {}".format(row))

def AtomNumberToAtomName(atomnumber:int) -> str:
    if len(atomprops) == 0:
        get_atomprops()
    for k in atomprops.keys():
        if atomprops[k]["atomicnumber"] == atomnumber and atomprops[k]["charge"] == 0:
            return k
    sys.exit("Invalid atomnumber %d" % atomnumber)

def AtomNameToAtomNumber(atomname:str) -> int:
    if len(atomprops) == 0:
        get_atomprops()
    for k in atomprops.keys():
        if atomname.upper() == k.upper():
            return atomprops[k]["atomicnumber"]
    sys.exit("Invalid atomname %s" % atomname)

def ElementName(atomname:str) -> str:
    if len(atomprops) == 0:
        get_atomprops()
    for k in atomprops.keys():
        if atomname.upper() == k.upper():
            return AtomNumberToAtomName(AtomNameToAtomNumber(atomname))
    sys.exit("Invalid element %s" % atomname)

def StringIsElement(atomname:str) -> bool:
    return atomname in atomprops

