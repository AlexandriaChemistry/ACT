# qgen

This directory implements partial atomic charge generation using the Alexandria Charge Model (ACM) and the Restrained Electrostatic Potential (RESP) method.

## Contents

- **qgen_acm.h / qgen_acm.cpp** — `QgenAcm`: generates charges via the ACM, which solves a linear system that balances electronegativity-equalisation constraints with Slater/Gaussian Coulomb interactions between atom-centred charge distributions. Bond indices within fragments are 0-based and exclude shell particles; `nonFixed_[b.aI()]` maps pre-build core positions to post-build atom indices when polarisation shells are interleaved.
- **qgen_resp.h / qgen_resp.cpp** — `QgenResp`: fits charges by minimising the least-squares deviation between the electrostatic potential computed from point charges and a reference potential (from QM), optionally with Bayesian restraints.
- **qgen_cube.cpp** — Reads Gaussian cube files containing electrostatic potential grids; supplies the reference ESP data for RESP fitting.
- **qtype.h / qtype.cpp** — `QType` enum and utilities that classify atom charge types (core, shell, virtual site) and control which atoms participate in each charge-generation step.
