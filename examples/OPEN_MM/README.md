Running Alexandria force field simulations with OpenMM
------------------------------------------------------
In this folder you find a python script that can be used to run OpenMM
(once you have installed it). The system is a sodium chloride crystal using
the force field due to Walz et al.
https://doi.org/10.1021/acs.jctc.8b00507

Use the python script with the following flags:

```
./run.py
```

As examples the following files are provided:
+ your_coordinate_file.pdb: MX_1000.pdb (a NaCl crystal lattice)
+ your_simulation_settings_file.dat: simulation.dat (a file specifying simulation settings)
+ your_force_field_parameter_file.xml  openmm_ff.xml (an example of an OpenMM compatible force field file that can be generated using alexandria gentop)

The final output line should be similar to
```
potential energy = -389092 kJ/mol
```
(some numerical differences may occur).

Proceed at your own risk and modify the scripts to your needs.

[OpenMM documentation](http://docs.openmm.org/development/userguide/).

[Generic information on OpenMM](https://openmm.org).



