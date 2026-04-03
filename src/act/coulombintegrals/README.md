# coulombintegrals

This directory provides analytical Coulomb integrals over Slater-type and Gaussian-type orbitals, used in the charge-model calculations of ACT.

## Contents

- **coulomb_gaussian.h** — Interface for Coulomb integrals computed with Gaussian basis functions.
- **gaussian_integrals.h** — Declarations for the individual Gaussian-orbital integral routines.
- **slater_integrals.h / slater_integrals.cpp** — High-level interface and dispatcher for Slater-type orbital (STO) Coulomb integrals; selects the correct pre-generated integral function based on orbital quantum numbers.
- **slater_low.h** — Low-level declarations shared by the generated Slater integral files.
- **Slater_nS_mS.cpp** — Pre-generated C++ source files (one per orbital-pair combination, e.g. `Slater_1S_1S.cpp`, `Slater_2S_3S.cpp`, …) containing closed-form expressions for the Coulomb interaction between two Slater *s*-type orbitals of principal quantum numbers *n* and *m*.
- **DSlater_nS_mS.cpp** — Corresponding files with the analytical derivatives of the Slater integrals with respect to interatomic distance.

These integrals are used by the charge-generation modules (`qgen`) to evaluate the electrostatic interaction energy between charge distributions represented as Slater or Gaussian functions, enabling parameter fitting without numerical quadrature.
