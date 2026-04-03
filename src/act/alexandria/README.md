# alexandria

This directory contains the high-level application layer of ACT: molecule handling, topology construction, force-field training, and the CLI sub-commands exposed by the `alexandria` executable.

## Contents

### Molecular representation and topology
- **actmol.h / actmol.cpp** — `ACTMol`: the central molecule object that combines a `MolProp` with a `Topology`, charge generation, and force-computation capabilities.
- **actmol_low.h / actmol_low.cpp** — Lower-level helpers shared by `ACTMol`.
- **topology.h / topology.cpp** — `Topology`: holds the atom list, bonded and non-bonded parameter tables, shells, and virtual sites for a single molecule.
- **fragmenthandler.h / fragmenthandler.cpp** — Manages fragmentation of a molecule for per-fragment charge generation.
- **allbondeds.h / allbondeds.cpp** — Collects and classifies all bonded interactions (bonds, angles, dihedrals, impropers) found in a topology.
- **dissociation_energy.h / dissociation_energy.cpp** — Reads and stores bond dissociation energies.

### Force-field training and optimisation
- **train_ff.h / train_ff.cpp** — Top-level routine for the `train_ff` CLI command; orchestrates genetic-algorithm or Bayesian optimisation of force-field parameters.
- **train_utility.h / train_utility.cpp** — Shared helpers for training (I/O, parameter management, convergence checks).
- **acm_ga.h / acm_ga.cpp** — ACM-specific specialisation of the genetic algorithm (`AcmGa`).
- **acmfitnesscomputer.h / acmfitnesscomputer.cpp** — Fitness function evaluator for ACM parameter optimisation.
- **acmindividual.h / acmindividual.cpp** — Genome individual carrying a set of force-field parameters.
- **acminitializer.h / acminitializer.cpp** — Population initialiser for ACM genetic algorithm runs.
- **bayes.h / bayes.cpp** — Bayesian (MCMC) parameter estimation for force-field training.
- **mcmcmutator.h / mcmcmutator.cpp** — MCMC-based mutation operator.
- **percentmutator.h / percentmutator.cpp** — Percentage-based mutation operator for the GA.
- **devcomputer.h / devcomputer.cpp** — Computes deviations between computed and reference molecular properties during training.
- **optimizationindex.h / optimizationindex.cpp** — Maps force-field parameters to their positions in the optimisation vector.
- **staticindividualinfo.h / staticindividualinfo.cpp** — Immutable per-run metadata shared across all GA individuals.
- **confighandler.h / confighandler.cpp** — Reads and validates training configuration options.

### CLI commands
- **alexandria.cpp** — Entry point; registers and dispatches all `alexandria` sub-commands.
- **alex_modules.h / alex_modules.cpp** — Registration of sub-commands (`gentop`, `train_ff`, `analyze`, `simulate`, …).
- **gentop.cpp** — `gentop` command: generates a topology from a structure or Gaussian output file.
- **analyze.cpp** — `analyze` command: analyses molecular and force-field properties and produces tables.
- **geometry_ff.cpp** — `geometry_ff` command: derives bond/angle/dihedral equilibrium distributions.
- **molgen.h / molgen.cpp** — Molecule generator used by several CLI commands.
- **molhandler.h / molhandler.cpp** — Reads molecule sets and prepares them for training or analysis.
- **molselect.h / molselect.cpp** — Selects a subset of molecules based on user-supplied criteria.
- **acthelper.h / acthelper.cpp** — Common CLI setup helpers (option parsing, MsgHandler initialisation).
- **actmiddleman.h / actmiddleman.cpp** — Mediates data flow between molecule objects and the training machinery.

### Output and export
- **pdbwriter.h / pdbwriter.cpp** — Writes molecular structures in PDB format.
- **gromacs_top.h / gromacs_top.cpp** — Exports topology in GROMACS `.top` / `.itp` format.
- **openmm_xml.h / openmm_xml.cpp** — Exports topology and force-field data as OpenMM XML.
- **symmetrize_charges.h / symmetrize_charges.cpp** — Enforces charge symmetry constraints on equivalent atoms.
- **thermochemistry.h / thermochemistry.cpp** — Computes thermochemical properties (ZPE, enthalpy, entropy) from normal modes.
- **princ.h / princ.cpp** — Principal-axes analysis of molecular geometry.
