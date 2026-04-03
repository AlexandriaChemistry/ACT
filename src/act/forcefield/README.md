# forcefield

This directory implements the force-field data model: parameter storage, combination rules, potential functions, atom/particle types, and XML-based I/O.

## Contents

### Core data model
- **forcefield.h / forcefield.cpp** — `ForceField` class: top-level container that maps `InteractionType` keys to `ForceFieldParameterList` objects; loaded from and saved to XML.
- **forcefield_parameter.h / forcefield_parameter.cpp** — `ForceFieldParameter`: a single typed parameter with value, uncertainty, minimum/maximum bounds, `Mutability`, and physical unit.
- **forcefield_parameterlist.h / forcefield_parameterlist.cpp** — `ForceFieldParameterList`: ordered collection of `ForceFieldParameter` objects, keyed by `Identifier`.
- **forcefield_parametername.h** — Canonical string names for standard force-field parameters (e.g. `"sigma"`, `"epsilon"`, `"charge"`).
- **particletype.h / particletype.cpp** — `ParticleType`: describes a named particle (atom type) together with its mass, element, and associated `ActParticle` classification.

### Potential functions and combination rules
- **potential.h / potential.cpp** — Enumeration of supported potential functions (Lennard-Jones, Buckingham, Morse, …) and helpers to evaluate them.
- **combinationrules.h / combinationrules.cpp** — Combination-rule implementations (Lorentz–Berthelot, geometric mean, etc.) for deriving pair interaction parameters.
- **combruleutil.h / combruleutil.cpp** — Higher-level utilities that apply combination rules across all particle-type pairs.
- **generate_dependent.h / generate_dependent.cpp** — Generates derived (dependent) parameters (e.g. pair parameters) from the base set using combination rules.
- **forcefield_utils.h / forcefield_utils.cpp** — Miscellaneous helpers: unit conversion, parameter lookup, validity checks.
- **forcefield_tables.h / forcefield_tables.cpp** — Builds tabulated parameter tables (e.g. for lookup during simulation).
- **symcharges.h / symcharges.cpp** — Reads and applies symmetry constraints that force chemically equivalent atoms to carry equal charges.

### Serialisation and file I/O
- **forcefield_xml.h / forcefield_xml.cpp** — Reads and writes force-field data in the ACT XML format.
- **act_checksum.h / act_checksum.cpp** — Computes and verifies a checksum of force-field XML files to detect inadvertent modifications.

### CLI commands
- **edit_ff.cpp** — Implementation of the `edit_ff` sub-command: manipulates and compares force-field XML files.
- **gen_ff.cpp** — Implementation of the `gen_ff` sub-command: generates a new force-field XML file from a specification.
- **merge_ff.cpp** — Implementation of the `merge_ff` sub-command: merges two or more force-field XML files.
