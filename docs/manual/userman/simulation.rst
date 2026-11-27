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

The angle potentials that depend on the distances between the atoms, the linear angles (Eqn. :eq:`vlin`) and the Urey-Bradley angles, are written as a *CustomAngleForce* with their energy expression and the parameters of each angle. Since the *CustomAngleForce* of OpenMM can evaluate the angle only, the file has a script for these that OpenMM runs when it makes the system. It replaces each such force by a *CustomCompoundBondForce* with the same expression, parameters and angles, named after the potential, for instance LINEAR_ANGLES. Nothing has to be done for this by the user. Proper dihedrals become a *PeriodicTorsionForce*, or an *RBTorsionForce* for the Fourier dihedrals, and the harmonic improper dihedrals a *CustomTorsionForce* with the central atom first. OpenMM applies an improper dihedral to every atom with three neighbours of the classes in the file, whereas ACT only makes one for an atom that is planar in the input structure.

-------------------
GROMACS topologies
-------------------

.. index::
   single: GROMACS; topology
   single: gentop; GROMACS topology

For force fields that GROMACS can represent, *alexandria gentop* writes a GROMACS topology as well::

    alexandria gentop -ff actff.xml -f dimer.pdb -p topol.top -c conf.pdb -ws

GROMACS has the force field as it is if the charges are point charges, the Van der Waals potential is Lennard-Jones 12-6 and the bonded potentials are ones that GROMACS has (harmonic and Morse bonds, harmonic, Urey-Bradley and linear angles, proper, Fourier and improper dihedrals). The comment at the top of the topology says what is missing otherwise: Gaussian or Slater charges become point charges, other Van der Waals potentials are left out, and a bond potential that GROMACS lacks is written without parameters, such that GROMACS does not accept the topology. The Lennard-Jones parameters of all pairs of particle types are given explicitly, which reproduces any combination rule.

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
