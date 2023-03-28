Running Alexandria force field simulations with OpenMM
------------------------------------------------------
In this folder you find a python script that can be used to run OpenMM
(once you have installed it). 

Use the python script with the following flags:

```
./run.py -pdb your_coordinate_file.pdb -dat your_simulation_settings_file.dat -xml your_force_field_parameter_file.xml 
```

As examples the following files are provided:
+ your_coordinate_file.pdb: MX_1000.pdb (a NaCl crystal lattice)
+ your_simulation_settings_file.dat: simulation.dat (a file specifying simulation settings)
+ your_force_field_parameter_file.xml  openmm_ff.xml (an example of an OpenMM compatible force field file that can be generated using alexandria gentop)

The output of the MD run is saved into the folder "output".

To run with the above inputs, please use:
```
./run.py -pdb MX_1000.pdb -xml openmm_ff.xml -dat simulation.dat
```
The final output line should be similar to
```
potential energy = -389092 kJ/mol
```
(some numerical differences may occur).

Proceed at your own risk and modify the scripts to your needs.

[OpenMM documentation](http://docs.openmm.org/development/userguide/).

[Generic information on OpenMM](https://openmm.org).



