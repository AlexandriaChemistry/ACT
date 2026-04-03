# properties

This directory calculates thermodynamic and structural molecular properties that go beyond single-point energies, including second virial coefficients, normal modes, and dimer interactions.

## Contents

- **secondvirial.h / secondvirial.cpp** — `SecondVirial`: computes the second virial coefficient *B*(T) by numerical integration over dimer interaction energies at a range of temperatures.
- **b2data.h / b2data.cpp** — Storage and I/O for second virial coefficient data (computed values and experimental references).
- **dimergenerator.h / dimergenerator.cpp** — `DimerGenerator`: systematically generates dimer configurations (orientations and separations) for two molecules, used in *B*(T) and interaction-energy calculations.
- **normalmodes.h / normalmodes.cpp** — Normal-mode analysis: diagonalises the Hessian to obtain harmonic vibrational frequencies and zero-point energies, and computes thermochemical corrections (ZPE, *H*, *S*, *G*).
- **rotator.h / rotator.cpp** — Rotates a molecule to a standard orientation (e.g. principal axes) before property calculations.
- **velocityhandler.h / velocityhandler.cpp** — Assigns or scales atomic velocities for MD initialisation or kinetic-energy analysis.

### CLI commands
- **simulate.cpp** — Implementation of the `simulate` sub-command: runs a short MD or energy-minimisation trajectory.
- **min_complex.cpp** — Implementation of the `min_complex` sub-command: energy-minimises a molecular complex and reports interaction energies.
