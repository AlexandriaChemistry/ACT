.. _sec-gui:

*****************************
The graphical user interface
*****************************

ACT-GUI is a user interface for ACT that runs in a web browser. A small server on your own computer runs the commands of *alexandria* through the Python package of :ref:`sec-python`, and a page shows the results: molecules in three dimensions with their charges, spectra, force fields, and how well a force field reproduces reference data. Long calculations can be sent from the page to a cluster on which you have an account. The interface is a repository of its own, `ACT-GUI`_, which takes over from the program *actgui* of ACT-utilities.

Everything the interface does, it does by running *alexandria* with options that it shows or that this manual describes, so a calculation that was tried in the interface can be repeated on the command line or in a script.

.. _ACT-GUI: https://github.com/AlexandriaChemistry/ACT-GUI

========================
Installing and starting
========================
The interface needs ACT with its Python package, and Python 3.9 or newer with the packages *fastapi* and *uvicorn*. To save plots as PDF files it uses `plotXVG`_ when that is installed. From a copy of the repository:

.. code-block:: bash

   pip install ".[plots]"
   act-gui

*act-gui*, or *python -m act_gui*, prints an address and opens it in the browser. The server listens on your own computer only, and the address has a key without which it does not answer. Uploaded files and results are kept in a temporary directory that is removed when the program ends; *--data DIR* keeps them in a directory. *--no-browser* does not open a browser and *--port* chooses the port.

The interface finds *alexandria* as the Python package does: the program the package was installed with, or the one in the environment variable *ACT_ALEXANDRIA*. Calculations on more than one process need *mpirun*, *mpiexec* or *srun* in the *PATH*, or the launcher in *ACT_MPIEXEC*.

.. _plotXVG: https://github.com/AlexandriaChemistry/plotXVG

======================
The molecule workspace
======================
A molecule is given as a SMILES string or as a file, and is read with the default force field, the force field of the force field workspace or a force field file. The viewer shows it coloured by element or by charge, and the tabs run the calculations for one molecule (:numref:`fig-gui-molecule`):

* *Charges*: the charges and the topology (*alexandria gentop*).
* *Minimize*: an energy minimization (*minimize*), with the minimized structure on top of the original one, the RMSD and the energy terms.
* *Normal modes*: the vibrations (*nma*) with the infrared spectrum, the animated modes and the thermochemistry.
* *Simulate*: a short molecular dynamics simulation (*simulate*), with the energies plotted while it runs and the trajectory animated afterwards.
* *OpenMM*: an archive with which the molecule is simulated with OpenMM and the same force field, see :ref:`sec-simulations`.

.. figure:: ../images/gui-molecule.png
   :name: fig-gui-molecule
   :width: 100%
   :align: center

   The molecule workspace after a normal mode analysis of acetone. The atoms are coloured by charge; a click on a line under the spectrum or on a row of the table animates that vibration.

=========================
The force field workspace
=========================
This workspace has one force field at a time, which is made there or opened from a file. Each step writes a new file that can be saved and that the other workspaces can use.

* *Generate* makes a new force field with *alexandria gen_ff* from a few choices: the charge model, the potentials with their combination rules, the atom types and which of them get a shell, and virtual sites.
* *Virtual sites* defines sites by picking two, three or four atoms of a molecule in the viewer and choosing the kind of site; the site is shown where it would be.
* *Geometry* takes bond lengths, angles and dihedral angles from the compounds of a molprop file (*geometry_ff*): a file of its own, or the file of the molprop workspace with the compounds that are ticked there.
* *Edit and validate* lists the parameters with their value, bounds and whether they are trained (:numref:`fig-gui-forcefield`), changes them with *edit_ff*, and shows the problems that *edit_ff -validate* finds.
* *Compare* lists what differs between the force field of the workspace and another one: parameters with another value, with other bounds or another way of training, and those that only one of the two has.

.. figure:: ../images/gui-forcefield.png
   :name: fig-gui-forcefield
   :width: 100%
   :align: center

   The parameters of a force field. The bar shows where a value is between its bounds; a click on a row opens the form to change the value, the bounds or whether the parameter is trained.

=======================
The molprop workspace
=======================
A molprop file with compounds and their reference data is opened and checked (*edit_mp -validate*). The page lists the compounds, the number of values of each property and the problems that were found. The table of compounds also says what each compound is used for: training, testing or nothing. A click changes that for one compound, and the buttons above the table for all compounds that the table shows.

**Selection files.** The choices of the table make the selection file that *train_ff* reads with *-sel*, see :ref:`sec-files`: a line for each compound that is used, with *Train* or *Test*. The box in front of a compound says whether it is used, and the box in the header does so for all compounds that the table shows, which can be limited to those that have values of one property. For a file with dimers the card *Selection file* lists the monomers, and ticking monomers chooses the dimers they are in: with the dimers A#B, A#C and B#C, ticking only A gives A#B and A#C but not B#C. The other setting uses a dimer only when all its monomers are ticked. *Draw at random* uses a percentage of the ticked compounds for testing. *Save* writes the file, named after the molprop file with the extension *.dat*, and *Open* reads a selection file that exists already: the table then says what the file says, and compounds of the molprop file that the selection does not have are not used. A comparison or a training, on your own computer or on a cluster, uses what the table says at that moment.

**Comparing.** A force field is compared with the data of the file by *train_ff -nooptimize*, for the properties that are chosen. The page shows the statistics for each property and data set and a plot of the values of the force field against the reference values (:numref:`fig-gui-comparison`), which can be saved as a PDF file. For a compound with several structures, such as a dimer, a histogram shows the energies of the reference and of the force field.

.. figure:: ../images/gui-comparison.png
   :name: fig-gui-comparison
   :width: 100%
   :align: center

   A force field compared with the interaction energies of dimers: each point is a structure, and the dashed line is where the force field equals the reference.

**Training.** With *Train* instead of *Compare* the force field is trained on the data (*train_ff*). The force field offers the kinds of parameters that are not fixed, each property has a weight, and the optimizer is Monte Carlo or the genetic algorithm with or without Monte Carlo, see :ref:`sec-performance` for what they need. *train_ff* has one set of parameters for each process, so the number of processes is the size of the population. The compounds for testing are those that the table says, or a percentage that is drawn at random with a fixed seed. Note that *train_ff* starts from random values of the parameters unless *From random values* is switched off.

The page follows a training with a plot of :math:`\chi^2` (:numref:`fig-gui-training`) and can stop it. Afterwards it lists the trained parameters with the result of the sensitivity analysis, compares the trained force field with the data, and shows what the training changed. The trained force field can be saved or opened in the force field workspace. A calculation on your own computer is stopped after ten minutes; longer ones belong on a cluster.

.. figure:: ../images/gui-training.png
   :name: fig-gui-training
   :width: 100%
   :align: center

   A short training: the lowest :math:`\chi^2` so far, and the parameters of the trained force field.

.. _sec-gui-cluster:

==================
Jobs on a cluster
==================
Long calculations, a comparison with a large file or a training, can run on a Slurm cluster on which you have an account. In the tab *Cluster* you give the name of the server, your user name there and a working directory; how ACT is made available in a job (a module, a file that the shell reads, a conda environment, or nothing); the partitions for jobs with and without a GPU; and the name of the program there, *alexandria* or *alexandria_d*. The settings are kept in the file *act-gui/settings.json* in your configuration directory.

The interface logs in with the *ssh* program of your computer, which uses your own keys and settings. It never asks for a password and keeps none, so the login has to work without a question: set up a key for the server, or open one connection in a terminal with connection sharing if the server asks for a password or a code each time. *Test the login* logs in once and tells whether the directory is there and the server has Slurm.

A cluster may block an account that logs in often, so the interface logs in as seldom as it can:

* Before anything is sent, the page shows the job: the script for Slurm, the command and the files that go along. The files get a name after their option, such as *ff.xml*, so no name of a file on your computer goes to the cluster.
* Sending a job is one login. The files go to a directory of the job below the working directory, where *sbatch* is started.
* One login asks for the state of all jobs together, not more often than the settings say (every two minutes by default), and only while the page is open. A training also tells how far it is (:numref:`fig-gui-cluster`).
* When a job has ended, one more login fetches its results, which are then shown as those of a calculation on your own computer.
* After a login that did not work the cluster is not asked again until you press the button.

What is known about the jobs is kept on your computer, so they are still listed when the interface is started again. *Stop* cancels a job with *scancel*. Nothing is ever removed on the cluster: the directory of a job stays there until you remove it.

The program on the cluster has to know the options that the *alexandria* of your computer uses, so the two should be of the same version; a training needs *train_ff -progress*.

.. figure:: ../images/gui-cluster.png
   :name: fig-gui-cluster
   :width: 100%
   :align: center

   The jobs on the cluster, with how far a training was when the cluster was last asked, and the settings below.
