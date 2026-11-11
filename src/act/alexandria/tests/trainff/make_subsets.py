#!/usr/bin/env python3
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

"""Make a small molprop file with a few conformations per molecule.

The molprop files for the train_ff tests were made from the benchmarks in
ACTdata/Performance with

  make_subsets.py booze/meth-eth-water.xml booze-inter.xml 8 'water#water' 'methanol#water' 'methanol#methanol'
  make_subsets.py booze/alc-mono.xml booze-mono.xml 8 methanol ethanol
  make_subsets.py hf-hf/sapt2+aug-cc-pvtz-0.03H-rmin0.9-HF_mp2-hh.xml hf-hf.xml 10 'hydrogen-fluoride' 'hydrogen-fluoride#hydrogen-fluoride'

The force fields lj-pc.xml, halgren-PG.xml, 4S-lj14-7-pg-hvdw.xml and
alc-pol-wbh.xml are copied from the same benchmarks.

usage: make_subsets.py input.xml output.xml nconf molname [molname ...]
"""
import sys
import xml.etree.ElementTree as ET

infile, outfile, nconf = sys.argv[1], sys.argv[2], int(sys.argv[3])
names = sys.argv[4:]
tree  = ET.parse(infile)
root  = tree.getroot()
for mol in list(root):
    if mol.get('molname') not in names:
        root.remove(mol)
        continue
    exps = mol.findall('experiment')
    # Evenly spaced conformations, to span the energy range of a scan
    keep = set(round(i*(len(exps)-1)/max(1, nconf-1)) for i in range(min(nconf, len(exps))))
    for i, e in enumerate(exps):
        if i not in keep:
            mol.remove(e)
        else:
            # The electrostatic potentials are not used in these tests
            for p in e.findall('potential'):
                e.remove(p)
ET.indent(tree, space='  ')
tree.write(outfile, encoding='utf-8', xml_declaration=True)
print(outfile, [ (m.get('molname'), len(m.findall('experiment'))) for m in root ])
