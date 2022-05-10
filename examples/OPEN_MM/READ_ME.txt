In this folder you find a python script that can be used to run OpenMM (once you have installed it). 
Use the python script with the following flags:

python run_openmm_simulation.py -pdb input_structures/your_coordinate_file.pdb -dat your_simulation_settings_file.dat -xml parameter_xml_files/your_force_field_parameter_file.xml 

As examples the following files are provided:
input_structures/your_coordinate_file.pdb:               MX_1000.pdb (a NaCl crystal lattice)
your_simulation_settings_file.dat:                       simulation.dat (a file specifying simulation settings)
parameter_xml_files/your_force_field_parameter_file.xml  openmm_example.xml (an example of an OpenMM compatible force field file that can be generated using alexandria gentop)

The output of the MD run is saved into the folder "output".

Proceed at your own risk and modify the scripts to your needs.


