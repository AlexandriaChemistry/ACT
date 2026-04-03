# import

This directory handles the import of molecular structures and computed properties from external chemistry software and file formats into the ACT data model.

## Contents

- **compound_reader.h / compound_reader.cpp** — `CompoundReader`: reads molecular structures and associated data using RDKit, supporting SDF, MOL, and related formats; converts the result into ACT's internal `MolProp` representation.
- **fetch_charges.h / fetch_charges.cpp** — Extracts partial atomic charges from external sources (e.g. Gaussian output files) and attaches them to a molecule.
- **import.h** — Public umbrella header for the import module.
- **import_utils.h / import_utils.cpp** — Shared helper functions for parsing and unit conversion used by the importer classes.
- **rdkit_io.cpp** — Low-level RDKit integration: calls RDKit APIs to parse molecular file formats and retrieves atom connectivity, coordinates, and properties.

The `import` module is the primary entry point for bringing new molecules into ACT workflows (e.g. via the `gentop` CLI command or the `gauss2molprop` script).
