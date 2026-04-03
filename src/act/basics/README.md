# basics

This directory contains fundamental data types and infrastructure used throughout the ACT toolkit.

## Contents

- **act_particle.h / act_particle.cpp** — Defines the `ActParticle` enum, which classifies simulation particles into types such as `Atom`, `Vsite`, `Shell`, and `SigmaHole`.
- **allmols.h / allmols.cpp** — Container for storing and retrieving sets of molecules by name.
- **atomization_energy.h / atomization_energy.cpp** — Reads and provides atomisation energies from a CSV data file.
- **atomprops.h / atomprops.cpp** — Atomic properties (element symbols, masses, van-der-Waals radii, etc.) read from `atomprops.csv`.
- **basecontainer.h** — Generic base class used as a common interface for various container types.
- **chargemodel.h / chargemodel.cpp** — Enumeration and utilities for the charge model (e.g. ACM, RESP) in use.
- **dataset.h / dataset.cpp** — Lightweight wrapper that associates a name and file path with a collection of molecular data.
- **identifier.h / identifier.cpp** — `Identifier` class representing a named, typed entity (atom type, interaction type, etc.) used as a key throughout the codebase.
- **interactiontype.h / interactiontype.cpp** — Central `InteractionType` enum listing all bonded and non-bonded interaction types (`BONDS`, `ANGLES`, `VDW`, `ELECTROSTATICS`, …) together with utility functions.
- **libraryfile.h / libraryfile.cpp** — Locates ACT data files (force-field XMLs, CSV tables) at runtime via the install path or `ACTDATA` environment variable.
- **msg_handler.h / msg_handler.cpp** — `MsgHandler` class for unified logging, warnings, and fatal-error reporting with configurable verbosity levels.
- **mutability.h / mutability.cpp** — `Mutability` enum (`Fixed`, `Bounded`, `Free`) that controls whether a force-field parameter may be modified during optimisation.
- **version.h / version.cpp** — Compile-time version string for the ACT library.
