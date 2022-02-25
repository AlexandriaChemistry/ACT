#!/usr/bin/env python3

import openbabel as ob

import subprocess
import os,sys, glob
sys.path.append('/Users/walz/GG/ACT/src/python')
import get_mol_dict as gmd
import molprops as mp
import numpy as np

def multipoles_pol_dict():
    Dipoles       = ["Dipole X", "Dipole Y", "Dipole Z"] 
    Quadrupoles   = ["Quadrupole XX", "Quadrupole XY", "Quadrupole XZ", "Quadrupole YY", "Quadrupole YZ", "Quadrupole ZZ"]        
    Octupoles     = ["Octupole XXX", "Octupole XXY", "Octupole XXZ", "Octupole XYY", "Octupole XYZ", "Octupole XZZ", "Octupole YYY", "Octupole YYZ", "Octupole YZZ", "Octupole ZZZ"]
    Hexadecapoles = ["Hexadecapole XXXX", "Hexadecapole XXXY", "Hexadecapole XXXZ", "Hexadecapole XXYY", "Hexadecapole XXYZ", "Hexadecapole XXZZ", "Hexadecapole XYYY", "Hexadecapole XYYZ", "Hexadecapole XYZZ", "Hexadecapole XZZZ", "Hexadecapole YYYY", "Hexadecapole YYYZ", "Hexadecapole YYZZ", "Hexadecapole YZZZ", "Hexadecapole ZZZZ"]   

    Multipoles = [Dipoles, Quadrupoles, Octupoles, Hexadecapoles]
    Multipoles_ACT = [[], [], [], []] # The Multipoles_ACT are in the same order as in Multipoles

    Polarizability = ['DIPOLE POLARIZABILITY XX', 'DIPOLE POLARIZABILITY XY', 'DIPOLE POLARIZABILITY XZ', 'DIPOLE POLARIZABILITY YY', 'DIPOLE POLARIZABILITY YZ', 'DIPOLE POLARIZABILITY ZZ']
    Polarizability_ACT = []
    
    return Multipoles, Multipoles_ACT, Polarizability, Polarizability_ACT  

def extract_grid_points(grid_file):
    x_values = []
    y_values = []
    z_values = []

    if os.path.isfile("{0}".format(grid_file)):     
        f = open(grid_file, 'r')
        for line in f:
            if (len(line.split()) == 3):
                x, y, z = line.split()
                x_values.append(float(x))
                y_values.append(float(y))
                z_values.append(float(z))
        f.close()         
    else:
        print("{0} is not available.".format(grid_file)) 

    return x_values, y_values, z_values



def write_molprops(output_name, ESP_results_dict, charge, multiplicity, filename, fileformat, gridfile, datasource, reference, program, method, basisset, conformation, jobtype, esp_data):

    molprops = mp.Molprops()
    for filename in ESP_results_dict:
        molecule_dict = gmd.get_mol_dict(filename, fileformat, forcefield="alexandria")
        bond_dict = molecule_dict["bonds"]

        x,y,z                    = extract_grid_points(gridfile)
        GRID = [x,y,z]
        
        ESP_values         = ESP_results_dict[filename]["ESP_values"]
        Dipoles            = ESP_results_dict[filename]["Multipoles"][0]
        Quads              = ESP_results_dict[filename]["Multipoles"][1]
        Octus              = ESP_results_dict[filename]["Multipoles"][2]
        Hexadecas          = ESP_results_dict[filename]["Multipoles"][3]
        Pols               = ESP_results_dict[filename]["Pols"]
        
        title   = str(filename[:-4])
        formula = molecule_dict["molecule"]["formula"]
        weight  ="{:.5f}".format(molecule_dict["molecule"]["mol_weight"])
           
        
        mp1 = mp.Molprop(str(title)) 
        mp1.add_prop("formula", str(formula))
        mp1.add_prop("mass", str(weight))
        mp1.add_prop("charge", str(charge))
        mp1.add_prop("multiplicity", str(multiplicity))
        

        for keys in bond_dict:
            (atom_1,atom_2) = keys
            bond_order = bond_dict[(atom_1, atom_2)]
            mp1.add_bond(str(atom_1), str(atom_2), str(bond_order))            


        exper = mp.Experiment(datasource, reference, program, method, basisset, conformation, jobtype, esp_data)

        unit_str_potential = "Hartree/e"
        unit_str_coord     = "Angstrom"
        temp = "0"

        for i in range(len(ESP_values)):
            esp_value = "{0:.5f}".format(ESP_values[i])
            x = "{0:.3f}".format(GRID[0][i])
            y = "{0:.3f}".format(GRID[1][i])
            z = "{0:.3f}".format(GRID[2][i])
            exper.add_potential(str(i), unit_str_potential, unit_str_coord, x, y, z, esp_value)

        for atom in molecule_dict["atoms"]: 
            element = molecule_dict["atoms"][atom]["atomic_number"]
            atomtype = molecule_dict["atoms"][atom]["atomtype"]
            
            x ="{:.5f}".format(molecule_dict["atoms"][atom]["X"])
            y ="{:.5f}".format(molecule_dict["atoms"][atom]["Y"])
            z ="{:.5f}".format(molecule_dict["atoms"][atom]["Z"])
            #print("{0} {1} {2} {3} {4} {5} {6}".format(str(element), atomtype, str(atom), str(x), str(y), str(z), "Angstrom"))
            exper.add_atom(str(element), atomtype, str(atom), "Angstrom", str(x), str(y), str(z))


        ###################################
        ######## Multipole moments ######## 
        ###################################

        # dipoles ### unit Debeye
        #Conversion of dipole unit:
        #L = 1.  Multiply by 2.5417464519 to convert to Debye
        conv_fact_dipole = 2.5417464519
        dipole_moment = "{:.5f}".format(np.sqrt((conv_fact_dipole * Dipoles[0])**2 + (conv_fact_dipole * Dipoles[1])**2 + (conv_fact_dipole * Dipoles[2])**2))
        dipole_x = "{:.5f}".format(conv_fact_dipole * Dipoles[0])
        dipole_y = "{:.5f}".format(conv_fact_dipole * Dipoles[1])
        dipole_z = "{:.5f}".format(conv_fact_dipole * Dipoles[2])
        exper.add_dipole("electronic", "Debye", str(temp), str(dipole_moment), str(0.000),str(dipole_x),str(dipole_y),str(dipole_z))

       
        # Traceless Quadrupoles ### unit B (D A)
        ### Conversion of quads
        #L = 2.  Multiply by 1.3450342976 to convert to Debye.ang
        conv_fact_quad = 1.3450342976
        
        # for the traceless quadrupole only the diagonal elements have to be modified: q_ii' = q_ii - 1/3 trace
        quad_xx = conv_fact_quad * Quads[0]
        quad_yy = conv_fact_quad * Quads[3]
        quad_zz = conv_fact_quad * Quads[5]
        trace = quad_xx + quad_yy + quad_zz

        quad_xx_traceless = "{:.5f}".format(quad_xx - (1/3)*trace)
        quad_yy_traceless = "{:.5f}".format(quad_yy - (1/3)*trace)
        quad_zz_traceless = "{:.5f}".format(quad_zz - (1/3)*trace)

        quad_xy = "{:.5f}".format(conv_fact_quad * Quads[1])
        quad_xz = "{:.5f}".format(conv_fact_quad * Quads[2])
        quad_yz = "{:.5f}".format(conv_fact_quad * Quads[4])
        exper.add_quadrupole("electronic","Buckingham",str(temp),str(quad_xx_traceless),str(quad_yy_traceless),str(quad_zz_traceless),str(quad_xy),str(quad_xz),str(quad_yz))


        ##########################################################################################################
        ### In psi4 and gaussian, both the octupole and the hexadecapole are given in there non-traceless form ###
        ##########################################################################################################
        ### Conversion of octs
        # L = 3.  Multiply by 0.7117614979 to convert to Debye.ang^2
        conv_fact_octs = 0.7117614979
        octus_xxx = "{:.5f}".format(conv_fact_octs * Octus[0])
        octus_xxy = "{:.5f}".format(conv_fact_octs * Octus[1])
        octus_xxz = "{:.5f}".format(conv_fact_octs * Octus[2])
        octus_xyy = "{:.5f}".format(conv_fact_octs * Octus[3])
        octus_xyz = "{:.5f}".format(conv_fact_octs * Octus[4])
        octus_xzz = "{:.5f}".format(conv_fact_octs * Octus[5])
        octus_yyy = "{:.5f}".format(conv_fact_octs * Octus[6])
        octus_yyz = "{:.5f}".format(conv_fact_octs * Octus[7])
        octus_yzz = "{:.5f}".format(conv_fact_octs * Octus[8])
        octus_zzz = "{:.5f}".format(conv_fact_octs * Octus[9])
        exper.add_octupole("electronic","D.Angstrom2",str(temp),str(octus_xxx),str(octus_xxy),str(octus_xxz),str(octus_xyy),str(octus_xyz),str(octus_xzz),str(octus_yyy),str(octus_yyz),str(octus_yzz),str(octus_zzz))


        ### Conversion of hexadecas
        # L = 4.  Multiply by 0.3766479641 to convert to Debye.ang^3 
        conv_fact_hexadeca = 0.3766479641

        hexadecas_xxxx = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[0])
        hexadecas_xxxy = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[1])
        hexadecas_xxxz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[2])
        hexadecas_xxyy = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[3])
        hexadecas_xxyz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[4])
        hexadecas_xxzz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[5])
        hexadecas_xyyy = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[6])
        hexadecas_xyyz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[7])
        hexadecas_xyzz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[8])
        hexadecas_xzzz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[9])
        hexadecas_yyyy = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[10])
        hexadecas_yyyz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[11])
        hexadecas_yyzz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[12])
        hexadecas_yzzz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[13])
        hexadecas_zzzz = "{:.5f}".format(conv_fact_hexadeca * Hexadecas[14])
        exper.add_hexadecapole("electronic","D.Angstrom3",str(temp),str(hexadecas_xxxx),str(hexadecas_xxxy),str(hexadecas_xxxz),str(hexadecas_xxyy),str(hexadecas_xxyz),str(hexadecas_xxzz),str(hexadecas_xyyy),str(hexadecas_xyyz),str(hexadecas_xyzz),str(hexadecas_xzzz),str(hexadecas_yyyy),str(hexadecas_yyyz),str(hexadecas_yyzz),str(hexadecas_yzzz),str(hexadecas_zzzz))

        # polarizabilites for dimer in Å3 (index 0)  0.148035889 : Polarizability (alpha): Cubic atomic units and cubic Angstroms.  1 au3 = (0.529)3 Angstroms3
        conversion = 0.148035889
        pol_xx = "{:.5f}".format(Pols[0]*conversion)
        pol_xy = "{:.5f}".format(Pols[1]*conversion)
        pol_xz = "{:.5f}".format(Pols[2]*conversion)
        pol_yy = "{:.5f}".format(Pols[3]*conversion)
        pol_yz = "{:.5f}".format(Pols[4]*conversion)
        pol_zz = "{:.5f}".format(Pols[5]*conversion)
        pol_average = "{:.5f}".format((Pols[0]*conversion + Pols[3]*conversion + Pols[5]*conversion)/3)
        exper.add_polarisability("electronic","Angstrom3",str(temp),str(pol_average),str(0.000), str(pol_xx),str(pol_yy),str(pol_zz),str(pol_xy),str(pol_xz),str(pol_yz))
 
        mp1.add_experiment(exper)
     
        molprops.add_molecule(mp1)
    molprops.write(output_name)


if __name__ == "__main__":

    write_molprops(output_name, ESP_results_dict, charge, multiplicity, filename, fileformat, gridfile, datasource, reference, program, method, basisset, conformation, jobtype, esp_data)

    