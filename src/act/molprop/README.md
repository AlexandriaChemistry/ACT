# molprop

This directory defines the molecular property data model: the structures used to store experimental and computed molecular data, together with serialisation to XML and SQLite, and related CLI utilities.

## Contents

### Core data model
- **molprop.h / molprop.cpp** — `MolProp`: central class that stores all properties of a single molecule (name, formula, charge, multiplicity, fragments, experiments, and computed observables).
- **experiment.h / experiment.cpp** — `Experiment`: one set of measurements or QM calculations for a molecule, each containing a collection of `MolPropObservable` entries.
- **fragment.h / fragment.cpp** — `Fragment`: a chemically meaningful sub-unit of a molecule used during charge generation and topology building.
- **molpropobservable.h / molpropobservable.cpp** — `MolPropObservable`: a named, typed property value (energy, dipole, quadrupole, polarizability, …) with optional uncertainty.
- **topologyentry.h / topologyentry.cpp** — `TopologyEntry`: a bonded-interaction entry (bond, angle, dihedral) stored in a `MolProp` before full topology construction.
- **composition.h / composition.cpp** — Elemental composition of a molecule; used for formula generation and atomisation-energy corrections.
- **phase.h / phase.cpp** — `Phase` enum (`Gas`, `Liquid`, `Solid`) and associated utilities.
- **categories.h / categories.cpp** — String-based category labels that classify molecules for training-set selection.
- **multipole_names.h / multipole_names.cpp** — Canonical string names for multipole components (Qxx, Qxy, …).

### Serialisation and I/O
- **molprop_xml.h / molprop_xml.cpp** — Reads and writes `MolProp` collections in the ACT XML format.
- **molprop_sqlite3.h / molprop_sqlite3.cpp** — Reads and writes `MolProp` collections from/to an SQLite3 database.
- **molprop_util.h / molprop_util.cpp** — Higher-level helpers: merging property files, filtering by category or phase, computing statistics across a dataset.
- **molprop_tables.h / molprop_tables.cpp** — Generates formatted (LaTeX/plain-text) tables of molecular properties.

### CLI command
- **edit_mp.cpp** — Implementation of the `edit_mp` sub-command: merges, filters, and edits molecular property XML files.
