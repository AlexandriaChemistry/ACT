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

.. _sec-simulate-thermo:

-------------------------------------------
Thermochemistry from simulations with beads
-------------------------------------------

.. index::
   single: thermochemistry; path integral
   single: entropy; path integral
   single: heat capacity; path integral
   single: anharmonicity

The normal mode analysis (:ref:`sec-nma`) gives the thermochemistry of a rigid rotor with harmonic oscillators. With the option *-thermo* the simulations with beads give the internal energy, the enthalpy, the heat capacity, the standard entropy and the Gibbs energy of one compound in the gas phase at the temperature *-temp* and 1 bar, without these approximations: with the anharmonicity of the model, the coupling between rotation and vibration, and all conformations that the simulation visits. The compound is minimized first, and the results are printed next to those of the normal mode analysis of the minimum. The energies are relative to that minimum. For example::

  alexandria simulate -ff ff.xml -f hf.sdf -nbeads 32 -temp 298.15 -nsteps 4000000 -tau 0.01 -thermo

The internal energy and the heat capacity are those of the simulation at that temperature, the heat capacity from the fluctuations of the energy. The entropy cannot be had from one simulation. It is computed by thermodynamic integration from the classical rigid rotor with harmonic oscillators, see :ref:`sec-pimd-thermo`, which needs 2 x *-nquad* more simulations, eight by default: classical ones at lower temperatures, which are longer, and ones with fewer beads in which the constant of Planck is scaled down. Together they take two to three times as long as the simulation itself. With *-nquad 0* they are left out and there is no entropy. With one bead the thermochemistry is that of the classical compound.

A finite number of beads gives too low an energy and too high a heat capacity: for water at room temperature with 32 beads 1.8 kJ/mol and 11 J/mol K. The values of the simulation therefore include the difference between harmonic oscillators with infinitely many beads and with the beads of the simulation, which the program prints. What remains is the effect of the number of beads on the anharmonicity, which is small if that correction is small compared to the zero point energy.

:numref:`tab-thermo-pimd` shows the result for hydrogen fluoride with an anharmonic bond, for which the exact values follow from the energy levels of the vibration and rotation in the same potential. The heat capacity needs the longest simulations, since it is a fluctuation.

.. table:: Thermochemistry of hydrogen fluoride with the potential of Hua for the bond at 298.15 K, relative to the minimum: the normal mode analysis, simulations of 4 000 000 steps with 32 beads, and the energy levels of the vibration and rotation.
   :name: tab-thermo-pimd

   +-----------------------------+----------+--------------------+--------------+
   | Property                    | Harmonic | Simulation         | Exact levels |
   +=============================+==========+====================+==============+
   | Internal energy (kJ/mol)    | 31.118   | 30.786 +/- 0.014   | 30.778       |
   +-----------------------------+----------+--------------------+--------------+
   | :math:`C_V` (J/mol K)       | 20.79    | 20.91 +/- 0.32     | 20.82        |
   +-----------------------------+----------+--------------------+--------------+
   | :math:`S^o` (J/mol K)       | 173.64   | 173.77 +/- 0.06    | 173.83       |
   +-----------------------------+----------+--------------------+--------------+
   | Gibbs energy (kJ/mol)       | -18.173  | -18.545 +/- 0.010  | -18.570      |
   +-----------------------------+----------+--------------------+--------------+

The identical nuclei of a compound are not exchanged in the simulation. As in the normal mode analysis, the rotational symmetry number accounts for them, which is correct unless the temperature is close to the rotational temperature, as for hydrogen below room temperature. The entropy is that of the conformations that the classical simulations reach from the minimum: for a compound with high barriers between conformations they have to be long enough to cross them.

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

-------------------
GROMACS topologies
-------------------

.. index::
   single: GROMACS; topology
   single: gentop; GROMACS topology

For force fields that GROMACS can represent, *alexandria gentop* writes a GROMACS topology as well::

    alexandria gentop -ff actff.xml -f dimer.pdb -p topol.top -c conf.pdb -ws

GROMACS has the force field as it is if the charges are point charges, the Van der Waals potential is Lennard-Jones 12-6 and the bonded potentials are ones that GROMACS has (harmonic and Morse bonds, harmonic, Urey-Bradley, quartic and linear angles, proper, Fourier and improper dihedrals). The comment at the top of the topology says what is missing otherwise: Gaussian or Slater charges become point charges, other Van der Waals potentials are left out, and a bond potential that GROMACS lacks is written without parameters, such that GROMACS does not accept the topology. The Lennard-Jones parameters of all pairs of particle types are given explicitly, which reproduces any combination rule.

Shells are written with their polarization, and virtual sites in the definitions of GROMACS, where the site comes first. A virtual site on four atoms is written as a center of weights, which GROMACS has for weights that are not negative only. Each compound in the input is a molecule in the topology, and compounds that are the same, with the same charges, share a molecule type. The flag *-ws* puts the shells and virtual sites in the structure file, which GROMACS needs, with the same names and in the same order as in the topology.

The energy of GROMACS differs by a constant from that of ACT if the bonds have an energy at their minimum, which is noted in the topology. Apart from that, GROMACS gives the energies and forces of ACT.

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
