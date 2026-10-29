# forces

This directory implements the computation of energies and forces from a `Topology` and a set of atomic coordinates. All calculations use the SIMD routines in `forces_simd`, which compute many conformations of one topology together.

## Contents

- **forcecomputer.h / forcecomputer.cpp** — `ForceComputer`: the public interface for energy and force evaluation, including the relaxation of shells and the treatment of virtual sites.
  - `computeTotal` computes the total energies and forces of a `ConformationSet`, the conformations of one topology. `ShellOptions` sets which shells are relaxed, whether they are reset first and an optional electric field.
  - `computeInteraction` computes the interaction energies and forces of many conformations of a dimer.
  - `computeRelaxedEnergies` computes the energies of the dimers and monomers of molecules for training.
  - `compute` and `computeOnce` are convenience functions for a single structure, with and without relaxing the shells.
- **forces_simd.h / forces_simd.cpp** — The implementation of all potentials, for one conformation at a time or with the conformations in the inner loop, which is vectorized.
- **forces_explicit_simd.h / forces_explicit_simd.cpp** — Explicit SIMD versions of the most expensive potentials, see `simd_avx2.h`, `simd_avx512.h`, `simd_neon.h` and `simd_math.h`. `simd_information.h` reports which instruction set is used.
- **conformationdata.h** — `ConformationData`: the coordinates, forces and energies of the conformations of one topology, stored such that they can be vectorized, and the energies that are reused if nothing changed.
- **forcecomputerstatistics.h** — Statistics of the number of evaluations and of shell relaxation.
- **forcecomputerutils.h / forcecomputerutils.cpp** — Shared mathematical helpers (distance and angle calculations, potential-energy derivatives, etc.).
- **vsitehandler.h / vsitehandler.cpp** — `VsiteHandler`: constructs virtual-site positions from the positions of the real atoms that define them, and distributes virtual-site forces back onto those real atoms.
