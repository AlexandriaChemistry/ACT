.. _sec-simulations:
        
*******************************************
MD Simulations with Alexandria Force Fields
*******************************************

.. index::
   single: molecular dynamics
   single: MD simulation

==================
Using ACT simulate
==================

.. index::
   single: energy minimization

The command *alexandria simulate* will perform a simulation in the gas phase, i.e. without periodic boundary conditions. This obviously limits the usefulness but the OpenMM software can be used for condensed-phase simulations using ACT force fields :cite:p:`Kriz2024a,Spoel2025b`.
The *alexandria simulate* utility is aimed for validation of the potentials, for instance by evaluating whether the total energy is conserved, it can be verified that the force is the derivative of the potential. The integration algorithm and the relaxation of shells are described in :ref:`sec-integration`.

The energy of a molecule or small cluster can be minimized using *alexandria minimize*::

    alexandria minimize -ff actff.xml -f coords.pdb 
    -charges mp2.xml -c afterem.pdb -g minimize.log

where *actff.xml* is the force field file, *coords.pdb* is a structure file (*xyz* and *sdf*
are supported as well). The *-charges mp2.xml* option indicates a molprop file (database) with monomeric (optimized) structures used for generating charges. The files *afterem.pdb* and *minimize.log* are output files. The minimization algorithm is selected with *-minalg*, atoms can be frozen using an index file with *-freeze* and bond lengths or angles can be constrained with *-constrain*, see :ref:`sec-constraints`. The same minimization can be done using *alexandria simulate -minimize*, which can in addition run a simulation starting from the minimized structure.

.. _sec-simulate-beads:

---------------------------------------------------------
Simulations at constant temperature and zero point energy
---------------------------------------------------------

.. index::
   single: path integral molecular dynamics
   single: zero point energy
   single: binding energy
   single: thermostat

By default *alexandria simulate* conserves the energy. With the option *-nbeads* the simulation is at the constant temperature *-temp* instead. With one bead it is classical molecular dynamics with a thermostat. With more beads it is path integral molecular dynamics, in which each atom is a ring polymer of beads, and which gives the nuclear quantum effects on the energies and the structure, such as the zero point energy, for any model, also with shells and virtual sites. The algorithm is described in :ref:`sec-pimd`, the same ring polymers are used for the second virial coefficient of flexible molecules (:ref:`sec-b2-flexible`). The charges are those of the input structure. For example, a water dimer at 100 K::

    alexandria simulate -ff actff.xml -f dimer.pdb -charges mp2.xml
    -minimize -nbeads 64 -temp 100 -nsteps 400000 -deltat 0.00025
    -nstener 20 -nstxout 1000 -monomers

The number of beads that is needed is a few times :math:`\hbar\omega/k_BT` for the highest frequency :math:`\omega` in the system: for water 32 beads at room temperature, and proportionally more at lower temperatures. Energies converge more slowly with the number of beads than the second virial coefficient does, please check the convergence. The time step of 0.25 fs suits flexible water, the default of 0.2 fs is safe.

The log file, and the JSON file if *-json* is given, have the averages of the potential, kinetic and total energy, of the energy terms, and of the spread of the beads, over the steps after the first *-equil* steps, a tenth of the simulation by default. The standard errors are from eight consecutive blocks of the simulation. They are too small if the simulation is not much longer than the slowest motion, for instance compounds that bind and separate. The kinetic energy includes the translation of the center of mass, see Eqn. :eq:`centroidvirial`, such that the total energy is that of the compound or cluster in the gas phase. The energy file, the progress file and the trajectory, which has the centroids of the beads, are written as without beads. Since the thermostat changes the motion of the atoms, the time dependence has no meaning, only the averages have.

Two or more compounds in vacuum that move apart do not come back. The program warns when that happens. A wall on the distance between the centers of mass of the compounds keeps them together, with the radius *-rwall* (nm) and the force constant *-kwall*, see Eqn. :eq:`fragmentwall`. The wall acts on the centers of mass only. For two compounds the option *-hist* writes the distribution of that distance, and the potential of mean force along it, :math:`W(R) = -k_BT\ln\left[p(R)/R^2\right]`, shifted to zero at its minimum. The potential of mean force at all distances, without a simulation of the dimer, is given by *alexandria b2 -pmf*, see :ref:`sec-b2-flexible`. At low temperature the simulation is the better source close to the minimum, since in *alexandria b2* few orientations then dominate the average.

With *-monomers* each compound is simulated alone as well, with the same settings, and the program prints the energy and the enthalpy of binding of :math:`n` compounds

.. math:: \Delta U(T) ~=~ \left<E\right>_{all} - \sum_{c=1}^{n}\left<E\right>_c, \qquad \Delta H(T) ~=~ \Delta U(T) - (n-1)\,k_BT,
   :label: bindingenergy

where :math:`E` is the total energy and the last term is the change of :math:`pV` of ideal gases. At low temperature :math:`\Delta U` becomes :math:`-D_0`, the binding energy with zero point motion, which can be compared to the harmonic estimate from the energies of the minimized structures and the zero point energies of *alexandria nma* (:ref:`sec-nma`). The standard error is that of the simulations together. These numbers are those of a bound cluster only as long as the compounds stay together by themselves. A water dimer at room temperature does not: the compounds then are apart most of the time, the energies depend on the volume that the wall leaves them, and the second virial coefficient is the quantity to compute. The program reports the fraction of the steps in which the compounds were at the wall, and warns if they reached it.

========================
The ACT-OpenMM interface
========================

.. index::
   single: OpenMM

The OpenMM :cite:p:`Eastman2023a` software in its native form is controlled from Python scripts. This allows users great control over the simulations. OpenMM allows to specify user-defined energy functions, which we use to implement both Gaussian-distributed charges :cite:p:`Ghahremanpour2018b,Walz2018a` and many Van der Waals potentials :cite:p:`Kriz2024a` (see Section :ref:`sec-energy`).
To facilitate this, a special python code was implemented that makes it relatively easy for the user to run simulations and minimization. The first step is to convert an ACT force field file (*actff.xml}) to one compatible with OpenMM (*openmmff.xml*)::

    alexandria gentop -ff actff.xml -openmm openmmff.xml 
    -charges mp2.xml -db "water methanol"

with two additional arguments. First, flag *-charges mp2.xml* indicates a molprop file (database) with monomeric (optimized) structures used for generating charges and, second, the flag *-db "water methanol"* instruct the code to produce an OpenMM topology for those compounds. If one or more of the compounds is not present in the database, a warning will be issued.

Each compound becomes a residue in the OpenMM force field file. With the flag *-c conf.pdb*, gentop writes a PDB file with the same residue and atom names, such as C1, H2 and so on. The residue names can be given with *-resnames "water=SOH methanol=MOH"*, or *-resnames "SOH MOH"* for the compounds in the order of the input. A name has at most four letters or digits, and cannot be a residue that OpenMM knows itself, such as HOH or SOL. Compounds without a name get one made from their IUPAC name. A simulation box made from these PDB files, for instance with packmol, should keep the residue and atom names. Then the ACT python interface tells OpenMM the bonds of the residues from the force field file, and the PDB file does not need CONECT records.

---------------------------
MD simulations using OpenMM
---------------------------

.. index::
   single: OpenMM; MD simulation
   single: ActOpenMMSim

It is easy to perform simulations based on an Alexandria force field. For this, we rely on the `OpenMM`_ software that you need to install separately. Then you can use the ACT python interface to OpenMM, as in this simple script, that we will call run.py::

  #!/usr/bin/env python3

  from act_openmm import ActOpenMMSim

  sim = ActOpenMMSim(pdbfile="file.pdb", 
                     datfile="dimer.dat",
                     xmlfile="ff.xml",
                     txtfile="output2.txt",
                     verbose=True)
  sim.run()
  sim.log_to_xvg("energy.xvg", [ "Potential Energy (kJ/mole)"])
  sim.log_to_xvg("temperature.xvg", [ "Temperature (K)" ])
  sim.log_to_xvg("density.xvg", [ "Density (g/mL)" ])

.. _OpenMM: https://openmm.org/

Here, we first instantiate an ActOpenMMSim object, and pass it a structure file *file.pdb*, a simulation parameter file *dimer.dat* and a force field file *ff.xml*. We also instruct the code that output should be written to *output2.txt* and that we want a lot of information in that file.
Then we can run it, that's all! The last three lines are for convenience of plotting the results. We can run this script using::

  python ./run.py    

or submitted to a cluster, preferably with GPUs available.

--------------------------------
Energy Minimization using OpenMM
--------------------------------

.. index::
   single: OpenMM; energy minimization

An energy minimization of the input structure can be done using this script::

  #!/usr/bin/env python3

  from act_openmm import ActOpenMMSim

  sim = ActOpenMMSim(pdbfile="file.pdb", 
                     datfile="dimer.dat",
                     xmlfile="ff.xml",
                     txtfile="output2.txt",
                     verbose=True)
  sim.setup()
  sim.minimize()
  sim.write_coordinates("final.pdb")

which is run here using the same flags as the simulation. Note that details about simulations and minimization can be specified in the *dimer.dat* file.
