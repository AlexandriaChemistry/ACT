# properties

This directory calculates thermodynamic and structural molecular properties that go beyond single-point energies, including second virial coefficients, normal modes, and dimer interactions.

## Contents

- **secondvirial.h / secondvirial.cpp** — `SecondVirial`: computes the second virial coefficient *B*(T) by numerical integration over dimer interaction energies at a range of temperatures.
- **b2data.h / b2data.cpp** — Storage and I/O for second virial coefficient data (computed values and experimental references).
- **dimergenerator.h / dimergenerator.cpp** — `DimerGenerator`: systematically generates dimer configurations (orientations and separations) for two molecules, used in *B*(T) and interaction-energy calculations.
- **ringpolymer.h / ringpolymer.cpp** — `RingPolymerSampler`: path integral molecular dynamics of one molecule as a ring polymer of beads, with shells relaxed and virtual sites constructed for each bead. It samples the structures of flexible monomers for *B*(T) with nuclear quantum effects (`b2 -nbeads`) and is the engine of simulations at constant temperature (`simulate -nbeads`).
- **pathintegralmd.h / pathintegralmd.cpp** — `pathIntegralMD`: simulation of a compound, dimer or cluster at constant temperature with ring polymers, with the energy file, the trajectory of centroids, averages with standard errors from blocks (`BlockAverager`), the histogram of the distance between two compounds and the binding energy (`simulate -nbeads`).
- **fragmentwall.h / fragmentwall.cpp** — `FragmentWall`: the distance between the centers of mass of the compounds of a dimer or cluster, and a wall on it that keeps them together in vacuum (`simulate -rwall`).
- **normalmodes.h / normalmodes.cpp** — Normal-mode analysis: diagonalises the Hessian to obtain harmonic vibrational frequencies and zero-point energies, and computes thermochemical corrections (ZPE, *H*, *S*, *G*).
- **rotator.h / rotator.cpp** — Rotates a molecule to a standard orientation (e.g. principal axes) before property calculations.
- **velocityhandler.h / velocityhandler.cpp** — Assigns or scales atomic velocities for MD initialisation or kinetic-energy analysis.

### CLI commands
- **simulate.cpp** — Implementation of the `simulate` sub-command: runs a short MD or energy-minimisation trajectory, at constant energy or, with beads, at constant temperature.
- **min_complex.cpp** — Implementation of the `min_complex` sub-command: energy-minimises a molecular complex and reports interaction energies.
