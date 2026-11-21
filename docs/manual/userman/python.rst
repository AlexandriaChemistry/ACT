.. _sec-python:

**********************
Using ACT from Python
**********************

The Python package *act* runs the commands of the *alexandria* program from Python: the options of a command are keyword arguments, and the results come back as Python objects and files. It is installed together with ACT, next to the program it starts, and needs no other package than the Python standard library. The graphical user interface, see :ref:`sec-gui`, is built on it.

The package starts the program as a separate process for each command. A second layer, the classes of ACT in Python, is being built: see :ref:`sec-python-classes` for what there is.

==================
Running a command
==================
Each command of *alexandria* is a function in *act.commands* with the name of the command:

.. code-block:: python

   import act

   res = act.commands.minimize(ff="ACS-pg1.xml", smiles="CC(=O)C")
   print(res.data["rmsd"])        # {'value': 0.0019, 'unit': 'nm'}
   print(res.file("c"))           # the file with the minimized structure

   print(act.commands.available())    # the commands
   help(act.commands.train_ff)        # the options of a command

The keyword arguments are the options of the command without the dash, so *alexandria minimize -ff ACS-pg1.xml -maxiter 500* becomes *minimize(ff="ACS-pg1.xml", maxiter=500)*. The rules are:

* An option that is a flag takes *True* or *False*: *optimize=False* gives *-nooptimize*.
* An option that takes several values or files takes a list.
* An option for an output file takes a file name, or *True* for the name that the program uses by default. Output files that the program always writes, such as the log file, need not be asked for.
* File names are relative to the directory from which Python runs, also when the command runs in another directory.

The options are checked before the program starts, against the description that the program gives of itself (*alexandria help -export json*). A name that is not an option of the command, a value of the wrong type, a text that starts with a dash where a value is expected, and an input file that does not exist raise an error at once.

*act.run("minimize", ff=..., smiles=...)* does the same as *act.commands.minimize(...)*, with the name of the command as a text.

===========
The result
===========
A command returns an object of the class *act.Result*:

* *data* has what the program wrote in JSON format: numbers are numbers, and a quantity with a unit is a dict with *value* and *unit*. Commands that can write JSON are asked to do so, unless *json=False* is given.
* *files* has the output files that exist, by the name of their option, each as a list; *file("c")* gives the first file of the option *-c*, or *None*.
* *log* has the text of the log file, *stdout* and *stderr* what the program printed, *args* the command line that was run, and *returncode* the exit status.
* *workdir* is the directory in which the program ran.

By default a command runs in a temporary directory, which is removed when the result is no longer used. Copy the files that should be kept, or give a directory, which is made if needed and kept:

.. code-block:: python

   res = act.commands.nma(ff="ACS-pg1.xml", smiles="O", ir=True, workdir="water-nma")
   print(res.data["Frequencies"])
   print(res.file("ir"))          # water-nma/IRspectrum.xvg

A result can be used in a *with* statement, which removes a temporary directory at the end of the block.

=======
Errors
=======
All errors are subclasses of *act.ActError*:

* *act.ActOptionError*: an unknown command or option, or a value that the option cannot take. Nothing was started.
* *act.ActNotFoundError*: the program, the launcher for several processes or an input file was not found. Nothing was started.
* *act.ActTimeoutError* and *act.ActCancelledError*: the program was stopped because it ran too long, or because the caller asked for it.
* *act.ActError* itself: the program ran and failed. The message has the errors that the program reported, and the exception carries *returncode*, *stdout*, *stderr*, *log* and *command_line*.

.. code-block:: python

   try:
       act.commands.minimize(ff="ACS-pg1.xml", smiles="not a molecule")
   except act.ActError as error:
       print(error)               # ... Cannot interpret the SMILES string ...
       print(error.returncode)    # 1

With *check=False* a program that fails raises no error: the result is returned and its *returncode* tells what happened.

=========================================
Long runs: progress, time limit, stopping
=========================================
*simulate* and *train_ff* can report how far they are while they run (their option *-progress*). The argument *on_progress* takes a function that is called with each report:

.. code-block:: python

   def report(item):
       if item["event"] == "generation":
           print(item["generation"], item["best"]["train"])

   res = act.commands.train_ff(ff="start.xml", mp="dimers.xml", sel="dimers.dat",
                               optimizer="HYBRID", pop_size=8, max_generations=50,
                               fit="sigma epsilon", fc_inter=1.0,
                               processes=8, on_progress=report)
   print(res.file("o"))           # the trained force field

For *simulate* a report has the step, the time, the temperature and the energies. For *train_ff* the key *event* tells what a report is: *start*, with the optimizer and the numbers of generations, iterations, individuals and parameters; *generation*, at the end of each generation, with the lowest :math:`\chi^2` so far in *best*; and *step*, a few times per second at most, for a step of the Monte Carlo chain of the master process.

A training that may take longer than the time it has, for instance in a queue, is told to stop in time with the option *maxh* of *train_ff*, the number of hours since the program started. It then ends as a complete training does, with the force field and, for the genetic algorithm, the gene pool from which a following training continues (*gpin*); *res.data["time_limit"]["reached"]* tells whether that happened. This is not the argument *timeout* below, which stops the program without results.

Other arguments of every command, which are not options of the program:

* *processes*: the number of processes. More than one starts the program with the MPI launcher, see below.
* *timeout*: the number of seconds after which the program is stopped.
* *stop*: an object with a method *is_set()*, such as a *threading.Event*, with which another thread stops the program.
* *env*: environment variables for the program, on top of those of Python.

==========================================
Which program, and how it is started
==========================================
The package starts the *alexandria* it was installed with. Otherwise it takes *alexandria*, or *alexandria_d* of a build in double precision, from the *PATH*. Another program, or a command that starts it such as that of a container, is chosen with the environment variable *ACT_ALEXANDRIA* or in Python:

.. code-block:: python

   act.config.executable = "/opt/act-test/bin/alexandria"
   act.config.mpi_launcher = ["srun", "--mpi=pmix"]
   act.config.timeout = 3600

For several processes the package looks for *mpirun*, *mpiexec* or *srun*; the environment variable *ACT_MPIEXEC* or *act.config.mpi_launcher* names another launcher, and *act.config.mpi_processes_flag* the option with which it takes the number of processes, *-n* by default.

==========================================
The options of the program as data
==========================================
*act.program_info()* gives the description of all commands that the package itself uses: for each command its options with their type, default value, allowed values and help text. A program that builds forms or command lines can use it:

.. code-block:: python

   command = act.program_info().find_command("minimize")
   for option in command.options.values():
       print(option.name, option.type, option.default, option.description)

*act.build_arguments(command, options, workdir)* makes the list of arguments for a command without running it, with the output files that the command would write. This is how the graphical user interface makes the command of a job that runs on a cluster.

.. _sec-python-classes:

========================
The classes of ACT
========================
When ACT is built with *-DACT_PYTHON=ON*, which needs `pybind11`_, the package also has classes that use the code of ACT directly, without a program in between. *act.classes_available()* tells whether an installation has them; the conda package *alexandria-python* does. There are the force field, the molecule with its energy, its minimum, its vibrations and simulations of it, and the compounds of a molprop file. All functions and classes of the package are described in :ref:`sec-python-api`.

.. code-block:: python

   import act

   ff = act.ForceField.read("ACS-pg1.xml")
   print(ff.particle_types[:5], ff.interactions, ff.polarizable)

   for p in ff.parameters("VANDERWAALS")[:3]:
       print(p.identifier, p.name, p.value, p.unit, p.minimum, p.maximum, p.mutability)

   sigma = ff.parameter("VANDERWAALS", "Ar", "sigma")
   ff.set_parameter("VANDERWAALS", "Ar", "sigma", value=0.3, mutability="Fixed")
   ff.set_particle_parameter("c3_s", "charge", minimum=-3.0)
   print(ff.particle_type("c3_s").parameters)
   print(ff.validate())           # the problems that edit_ff -validate lists
   ff.write("changed.xml")

A parameter that is asked for is a copy of what the force field has at that moment; it is changed through the force field, which checks the change first: the value stays within the bounds, the minimum below the maximum, and a parameter that cannot be negative stays so. A change that is refused changes nothing.

A molecule is made for a force field, from a SMILES string or a structure file, with its topology and charges as the programs of ACT make them:

.. code-block:: python

   mol = act.Molecule.from_smiles("CC(=O)C", ff)      # or act.Molecule.from_file("acetone.sdf", ff)
   print(mol.formula, mol.charge, mol.mass)
   for atom in mol.atoms[:4]:
       print(atom["name"], atom["type"], atom["particletype"], atom["charge"])
   print(mol.bonds[0])            # {'atoms': [0, 2], 'order': 1}
   print(mol.coordinates[0])      # x, y and z in nm
   energy = mol.energy()          # by term, in kJ/mol
   print(energy["EPOT"], energy["COULOMB"])
   forces = mol.forces()          # on each particle, in kJ/mol nm

The energy is minimized and the vibrations are computed as *alexandria minimize* and *alexandria nma* do it:

.. code-block:: python

   res = mol.minimize()           # the molecule gets the structure of the minimum
   print(res.converged, res.status, res.rmsd, res.energy_after["EPOT"])
   res = mol.minimize(algorithm="Newton", max_iter=500, tolerance=1e-8, freeze=[0, 2])

   modes = mol.normal_modes()     # minimizes first, unless minimize=False
   print(modes.frequencies)       # in 1/cm
   print(modes.intensities, modes.intensity_unit)
   print(modes.modes[0]["reduced_mass"], modes.modes[0]["displacements"][0])
   print(modes.thermochemistry["T=298.15"]["Standard entropy"]["Total"])
   x, y = modes.spectrum(linewidth=24)

A molecular dynamics simulation, as *alexandria simulate* does it, starts from the structure of the molecule and leaves it as it is:

.. code-block:: python

   def follow(item):
       print(item["step"], item["time"], item["temperature"], item["energies"]["EPOT"])
       return item["temperature"] < 2000      # False ends the simulation

   traj = mol.simulate(steps=10000, timestep=0.0002, temperature=300, seed=1993,
                       energy_interval=100, frame_interval=500, on_energies=follow)
   print(traj.energy_times, traj.energies["ETOT"], traj.temperature)
   print(traj.times, traj.frames[-1][0])      # the structures, in nm

The times are in ps. Without a temperature the atoms start at rest; with one they get velocities for it, the same each time for the same seed. Other threads of Python run while the simulation does, and Ctrl-C ends it the next time the energies are due.

The charges of ACT depend on the structure. They are made when the molecule is made, and again on request, for the structure the molecule has then; they can also be set by hand. The files for a simulation with OpenMM, see :ref:`sec-simulations`, are those of *alexandria gentop -openmm*:

.. code-block:: python

   mol.minimize()
   mol.generate_charges()                     # for the minimized structure
   mol.set_charges([0.0] * len(mol.atoms))    # or one charge for each particle
   mol.generate_charges("SQE")

   mol.write_openmm("acetone-openmm.xml", structure="acetone.pdb")

A minimization that does not converge is told in the result, with the structure of the lowest energy that was found; *update=False* leaves the structure of the molecule as it is. Normal modes of a structure that is no minimum make no sense, so there a minimization that does not converge raises an error.

The particles are the atoms, shells and virtual sites of the molecule, in the order of the coordinates. The numbers are those that *alexandria gentop*, *minimize* and *nma* report for the same molecule and force field. *act.data_directory()* tells where the classes look for the data files of ACT, after the directory in the environment variable *ACTDATA*. A molecule keeps the force field it was made for, with the parameters that the force field had at that moment.

The compounds of a molprop file, see :ref:`sec-files`, with their reference data:

.. code-block:: python

   mps = act.read_molprop("alcohols.xml")
   print(len(mps), mps.names, mps.summary["properties"])
   water = mps["water"]                       # or mps[0]: a compound by name or by place
   print(water["formula"], water["charge"], water["fragments"], water["properties"])
   for calc in mps.calculations("water"):
       print(calc["method"], calc["basisset"], calc["jobtype"], calc["coordinates"][0])
       for prop in calc["properties"]:
           print(prop["observable"], prop["type"], prop["value"], prop["unit"])
   print(mps.validate())                      # the check of edit_mp -validate
   mps.write("copy.xml")

The coordinates are in nm and the values in the units that ACT computes with, which each property names. A property is a number (*value*) or a vector or tensor (*vector*); what it is not, is *None*.

Errors are exceptions, and Python lives on: *act.ActFileError* (also an *OSError*) for a file that is missing, is not what it should be or cannot be written, *act.ActInputError* (also a *ValueError*) for input that is not accepted, such as a parameter that does not exist, a SMILES string that is no molecule or a molecule that the force field has no types for, and *act.ActInternalError* for what should not happen inside ACT. All are subclasses of *act.ActError*.

The warnings and the information of ACT, such as an atom type that the force field does not have, go to the logger *act* of the module *logging*. Nothing is shown unless the program asks for it:

.. code-block:: python

   import logging
   logging.basicConfig(level=logging.INFO)   # or logging.WARNING for the warnings only

.. _pybind11: https://pybind11.readthedocs.io

====================================
The class ACT of earlier versions
====================================
Until ACT 1.9b the module *act.py* of ACT-utilities had a class *ACT* that wrote command lines for *alexandria*. It is now *act.legacy.ACT*, and *from act import ACT* still works, with a warning. New scripts should use the functions above.
