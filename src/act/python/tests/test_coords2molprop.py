#!/usr/bin/env python3
#
# This file is part of the Alexandria Chemistry Toolkit
# https://github.com/AlexandriaChemistry/ACT
#
"""
Tests of coords2molprop, that the compounds A and B of a dimer are found
from the residue numbers in pdb files and from the bonds in xyz files, see
https://github.com/AlexandriaChemistry/ACT/issues/691
They are skipped if RDKit or xmltodict cannot be imported.
"""

import os, shutil, subprocess, sys, tempfile, unittest
import xml.etree.ElementTree as ET

try:
    import rdkit, xmltodict
    available = True
except ImportError:
    available = False

testDir = os.path.dirname(os.path.abspath(__file__))
script  = os.path.join(os.path.dirname(testDir), "coords2molprop")

waterDimer = [ ( "O", 1, -0.067,  0.000,  1.494 ),
               ( "H", 1,  0.816,  0.000,  1.866 ),
               ( "H", 1,  0.069,  0.000,  0.539 ),
               ( "O", 2,  0.063,  0.000, -1.423 ),
               ( "H", 2, -0.407, -0.760, -1.772 ),
               ( "H", 2, -0.407,  0.760, -1.772 ) ]

def writePdb(fileName:str, atoms:list, charges:dict={}):
    """
    Write atoms as tuples of element, residue number and coordinates, with
    the formal charges of the ions, by element, in columns 79-80
    """
    with open(fileName, "w") as pdb:
        pdb.write("REMARK   Epot -20.5 kJ/mol\n")
        for i, ( elem, resnr, x, y, z ) in enumerate(atoms):
            charge = ""
            if elem in charges:
                charge = "%d%s" % ( abs(charges[elem]), "+" if charges[elem] > 0 else "-" )
            pdb.write("HETATM%5d %-4s MOL  %4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s%2s\n" %
                      ( i+1, elem, resnr, x, y, z, elem, charge ))
        pdb.write("END\n")

@unittest.skipUnless(available, "RDKit or xmltodict not available")
class Coords2MolpropTest(unittest.TestCase):

    def setUp(self):
        self.workDir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.workDir, ignore_errors=True)

    def run_script(self, inFile:str, chargeA:int=0):
        """ Return the exit code and the formula and atoms of the fragments """
        outFile = os.path.join(self.workDir, "out.xml")
        result  = subprocess.run([ sys.executable, script, "-i", inFile, "-nA", "a", "-nB", "b",
                                   "-qa", str(chargeA), "-o", outFile ],
                                 cwd=self.workDir, capture_output=True, text=True,
                                 env=dict(os.environ, PYTHONPATH=os.path.dirname(script)))
        fragments = []
        if result.returncode == 0:
            for frag in ET.parse(outFile).getroot().iter("fragment"):
                fragments.append(( frag.get("formula"), frag.text.split() ))
        return result.returncode, fragments

    def test_water_dimer(self):
        pdb = os.path.join(self.workDir, "dimer.pdb")
        writePdb(pdb, waterDimer)
        self.assertEqual((0, [ ( "H2O", [ "1", "2", "3" ] ), ( "H2O", [ "4", "5", "6" ] ) ]),
                         self.run_script(pdb))

    def test_interleaved(self):
        """ The atoms of the residues do not have to be consecutive """
        pdb = os.path.join(self.workDir, "interleaved.pdb")
        writePdb(pdb, [ waterDimer[i] for i in [ 0, 3, 1, 2, 4, 5 ] ])
        self.assertEqual((0, [ ( "H2O", [ "1", "3", "4" ] ), ( "H2O", [ "2", "5", "6" ] ) ]),
                         self.run_script(pdb))

    def bonds(self):
        """ Return the bonds in the output """
        outFile = os.path.join(self.workDir, "out.xml")
        return sorted(( int(b.get("ai")), int(b.get("aj")) ) for b in ET.parse(outFile).getroot().iter("bond"))

    def test_ion_pair(self):
        """
        A compound without bonds between its atoms, NaCl with water. The
        ions have formal charges and are not bonded
        """
        pdb = os.path.join(self.workDir, "nacl_water.pdb")
        writePdb(pdb, [ ( "Na", 1,  0.0,  0.00, 0.0 ), ( "Cl", 1,  2.4,  0.00, 0.0 ),
                        ( "O",  2, -2.3,  0.00, 0.0 ), ( "H",  2, -2.9,  0.76, 0.0 ),
                        ( "H",  2, -2.9, -0.76, 0.0 ) ],
                 charges={ "Na": 1, "Cl": -1 })
        self.assertEqual((0, [ ( "ClNa", [ "1", "2" ] ), ( "H2O", [ "3", "4", "5" ] ) ]),
                         self.run_script(pdb))
        self.assertEqual([ ( 4, 3 ), ( 5, 3 ) ], self.bonds())

    def test_charged_ion(self):
        """ A sodium ion with water """
        pdb = os.path.join(self.workDir, "na_water.pdb")
        writePdb(pdb, [ ( "Na", 1,  0.0,  0.00, 0.0 ), ( "O",  2, -2.3,  0.00, 0.0 ),
                        ( "H",  2, -2.9,  0.76, 0.0 ), ( "H",  2, -2.9, -0.76, 0.0 ) ],
                 charges={ "Na": 1 })
        self.assertEqual((0, [ ( "Na+", [ "1" ] ), ( "H2O", [ "2", "3", "4" ] ) ]),
                         self.run_script(pdb, chargeA=1))
        self.assertEqual([ ( 3, 2 ), ( 4, 2 ) ], self.bonds())

    def test_xyz(self):
        """ Without residues the compounds are found from the bonds """
        xyz = os.path.join(self.workDir, "dimer.xyz")
        with open(xyz, "w") as outf:
            outf.write("%d\nEpot -20.5 kJ/mol\n" % len(waterDimer))
            for elem, _, x, y, z in waterDimer:
                outf.write("%s %g %g %g\n" % ( elem, x, y, z ))
        self.assertEqual((0, [ ( "H2O", [ "1", "2", "3" ] ), ( "H2O", [ "4", "5", "6" ] ) ]),
                         self.run_script(xyz))

    def test_three_residues(self):
        """ A dimer should have two compounds """
        pdb = os.path.join(self.workDir, "trimer.pdb")
        writePdb(pdb, waterDimer + [ ( elem, 3, x + 3, y, z ) for elem, _, x, y, z in waterDimer[:3] ])
        self.assertNotEqual(0, self.run_script(pdb)[0])

if __name__ == "__main__":
    unittest.main()
