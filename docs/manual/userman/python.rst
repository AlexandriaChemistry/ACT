.. _sec-python:

**********************
Using ACT from Python
**********************

The Python package *act* runs the commands of the *alexandria* program from Python: the options of a command are keyword arguments, and the results come back as Python objects and files. It is installed together with ACT, next to the program it starts, and needs no other package than the Python standard library. The graphical user interface, see :ref:`sec-gui`, is built on it.

The package starts the program as a separate process for each command. It does not give access to the data structures of ACT inside Python; that is planned as a second layer.

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

====================================
The class ACT of earlier versions
====================================
Until ACT 1.9b the module *act.py* of ACT-utilities had a class *ACT* that wrote command lines for *alexandria*. It is now *act.legacy.ACT*, and *from act import ACT* still works, with a warning. New scripts should use the functions above.
