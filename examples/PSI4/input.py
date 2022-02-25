#! Reference calculation of electrostatic potential, multipole moments and polarisability .

import psi4
import write_xml_ACT as ACT

psi4.set_memory('1 GB')
psi4.set_num_threads(2)
psi4.core.set_output_file('output.dat', True)



molecule_name = "water_Li"

water_Li = psi4.geometry("""
 noreorient
 nocom
0 1
O            0.0     0.0         0.117
H            0.0     0.7634     -0.4681
H            0.0     -0.7634     -0.4681
--
1 1
Li            0.0     0.0     2.197
""")


psi4.basis_helper("""
assign aug-cc-pvtz
""")

psi4.set_options({"PROPERTIES_ORIGIN":["NUCLEAR_CHARGE"]}) 

########## ESP at Grid ############### ########## polarizability #############   
E, wfn = psi4.prop("b3lyp", properties=["GRID_ESP", "DIPOLE_POLARIZABILITIES"], return_wfn=True)
ESP_grid = wfn.oeprop.Vvals()

######### Dipoles, Quadrupoles, Octupoles, Hexadecapoles ##############
psi4.oeprop(wfn, 'MULTIPOLE(4)',  title='multipoles') 


Multipoles, Multipoles_ACT, Polarizability, Polarizability_ACT = ACT.multipoles_pol_dict()

for pol in range(len(Polarizability)):
    Polarizability_ACT.append(psi4.variable(Polarizability[pol]))

for multipole in range(len(Multipoles)):
    for property in Multipoles[multipole]:
        Multipoles_ACT[multipole].append(wfn.variable(property))


######### Information for .xml file ##############
datasource = "Theory"
reference  = "N.A."
program    = "Psi4/1.3.2"
method     = "B3LYP"
basisset   = "aug-cc-pVTZ"
conformation = "excited" 
jobtype    = "SP"
    
tot_charge = "1"
multiplicity = "1"

gridfile    = "grid.dat"
esp_data    = "grid_esp.dat"
filename    = "water_Li.xyz"
output_name = "ACT_ref_data.xml"

ESP_results_dict = {}
ESP_results_dict[filename] = {"ESP_values": ESP_grid, "Multipoles": Multipoles_ACT, "Pols": Polarizability_ACT}

### Write .xml file for ACT
ACT.write_molprops(output_name, ESP_results_dict, tot_charge, multiplicity, filename, fileformat, gridfile, datasource, reference, program, method, basisset, conformation, jobtype, esp_data)

