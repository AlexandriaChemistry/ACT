# python

This directory contains Python scripts and utilities that complement the C++ `alexandria` executable, providing convenience wrappers, analysis tools, and integration with third-party packages.

## Contents

### Core libraries
- **act.py** — High-level Python API wrapping common `alexandria` command invocations.
- **actutils.py** — General utility functions shared by the other scripts (file handling, unit conversion, logging).
- **molutils.py** — Molecular structure helpers (geometry manipulation, formula parsing).
- **molprops.py** — Python representation of molecular properties; parses ACT XML/CSV output.
- **elements.py** — Periodic table data (atomic numbers, masses, covalent radii).
- **mol_csv_api.py** — API for reading and writing molecule data from CSV files.
- **get_mol_dict.py** — Builds a dictionary of molecules from a directory of structure files.
- **get_csv_rows.py** — Extracts rows matching given criteria from ACT CSV output files.
- **jacobi.py** — Jacobi eigenvalue solver used for inertia-tensor diagonalisation.
- **atomic_heat_of_formation.py** — Calculates heats of formation from atomisation energies and QM total energies.
- **read_gaussian_log.py** — Parses Gaussian 16/09 log files to extract energies, geometries, frequencies, and multipoles.

### OpenMM integration
- **act_openmm.py** — Runs OpenMM simulations using force-field parameters exported by `alexandria`.
- **act_gct.py** — Grand-canonical titration utility using OpenMM and ACT charge models.

### Converters and importers
- **gauss2molprop** (`gauss2molprop.cmakein`) — Converts Gaussian log files to ACT molprop XML.
- **ncia2molprop** — Converts NCIA (non-covalent interaction atlas) data to molprop XML.
- **donchev2molprop** — Converts the Donchev et al. dataset to molprop XML format.
- **coords2molprop** — Converts bare coordinate files to molprop XML.
- **generate_mp** (`generate_mp.cmakein`) — Batch-generates molprop XML files for a set of molecules.
- **ncia_reader_from_xyz.py** — Reads NCIA geometry/energy data from XYZ files.

### Analysis and visualisation
- **dofit.py** — Fits model curves to property data (e.g. force-field vs. QM energies).
- **actgui** — Graphical interface for visualising ACT training results.
- **view_fitness** — Plots fitness convergence curves from GA training output.
- **plot_convergence** — Plots parameter convergence during training.
- **dimer_scan** — Scans dimer interaction energies as a function of separation.

### Selection and configuration
- **molselect** (`molselect.cmakein`) — Command-line tool to select molecule subsets for training.
- **reshuffle_selection** (`reshuffle_selection.cmakein`) — Randomly reshuffles a molecule selection file.

Scripts ending in `.cmakein` are processed by CMake at build time to inject the correct installation paths.
