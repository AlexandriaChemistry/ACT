##################################################################
############ ACT reference data calculation using psi4 ###########
##################################################################

### Files ###
grid.dat
input.py
water_Li.xyz
write_xml.py
#############


The psi4 calculation is done be running the input.py script (psi4 has to be installed, e.g. by using anaconda, https://anaconda.org/psi4/psi4).

For the specified structure in the input.py file (in the section psi4.geometry) and level of theory, the electrostatic potential (ESP, calculated for the grid points specified in grid.dat), multipole moments (dipole, quadrupole, octupole, and hexadecapole) and polarisability will be calculated.

In the input file (input.py) you will find an example for a Li+ - water complex.
Certain information can be specified in the input.py file that will be written in the .xml file (e.g. level of theory).

The calculated reference data will be written into the .xml file that can be used in ACT as a reference data input file. The name of the .xml file can be specified in the input.py file.

The write_xml_ACT.py script will read the grid.dat and the .xyz file. The .xyz file is processed by a modified openbabel version (obenbabel-alexandria) that will extract information such as atomtypes and bonds that is also written in the .xml file. The ESP points, the multipole moments and the polarisability matrix are passed on directly to write_xml_ACT.py using a python dictionary. In this dictionary, reference data for several compounds can be added which can then be written to the same .xml file.   

Feel free to change the input.py file.

