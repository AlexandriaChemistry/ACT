# utility

This directory provides general-purpose support infrastructure used across the entire ACT toolkit: string and unit handling, XML/JSON serialisation, MPI communication, regression analysis, and LaTeX output.

## Contents

- **communicationrecord.h / communicationrecord.cpp** — `CommunicationRecord`: thin wrapper around MPI communicators; provides collective operations (broadcast, reduce, gather) used by parallel training and simulation runs.
- **jsontree.h / jsontree.cpp** — `JsonTree`: lightweight tree structure for reading and writing JSON configuration data. All values are stored and serialised as quoted strings; do not expect typed numeric values on read-back.
- **xml_util.h / xml_util.cpp** — Helper functions for reading and writing XML using libxml2 (node traversal, attribute access, error handling).
- **stringutil.h / stringutil.cpp** — String manipulation utilities: trimming, splitting, case conversion, number formatting, and tokenisation.
- **units.h / units.cpp** — Physical unit definitions and conversion factors; used throughout ACT to convert between internal (kJ mol⁻¹, nm, e) and external units.
- **regression.h / regression.cpp** — Ordinary least-squares linear regression; used for fitting property data and computing correlation metrics.
- **latex_util.h / latex_util.cpp** — Helpers for generating LaTeX table markup from numerical data, used by the `analyze` command.
- **memory_check.h / memory_check.cpp** — Lightweight memory-usage reporting utilities for profiling and debugging.
