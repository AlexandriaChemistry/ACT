# GitHub Copilot Instructions for the Alexandria Chemistry Toolkit (ACT)

## What is ACT?

ACT (Alexandria Chemistry Toolkit) is a C++20 scientific software package for developing and training molecular mechanics force fields. It is derived from a fork of [GROMACS](https://www.gromacs.org) and retains large portions of the GROMACS infrastructure (build system, utilities, test framework). The primary program is the `alexandria` executable, which exposes sub-commands for force-field training, molecular property analysis, topology generation, MD simulation, and more.

The canonical paper: *Digital Discovery* **4** (2025) 1925 — <https://doi.org/10.1039/D5DD00178A>

---

## Repository Layout

```
ACT/
├── src/
│   ├── act/                  # All ACT-specific C++ source code
│   │   ├── basics/           # Core types: MsgHandler, Identifier, ActParticle, InteractionType, …
│   │   ├── forcefield/       # ForceField, ForceFieldParameter, ParticleType, XML I/O
│   │   ├── molprop/          # MolProp (molecular properties), Experiment, Fragment, SQLite/XML I/O
│   │   ├── qgen/             # Charge generation: ACM (QgenAcm) and RESP (QgenResp)
│   │   ├── forces/           # ForceComputer, ForceComputerImpl, VsiteHandler
│   │   ├── alexandria/       # High-level classes: ACTMol, Topology, FragmentHandler,
│   │   │                     #   train_ff, analyze, gentop, simulate, … CLI modules
│   │   ├── ga/               # Generic genetic-algorithm framework
│   │   ├── import/           # RDKit-based compound reader, SDF/PDB/Gaussian log import
│   │   ├── properties/       # Second virial, normal modes, rotation, dimer generation
│   │   ├── statistics/       # Regression, statistical utilities
│   │   ├── utility/          # CommunicationRecord (MPI), JsonTree, XML utilities, units
│   │   └── coulombintegrals/ # Coulomb / Slater integrals
│   ├── gromacs/              # GROMACS infrastructure (math, utilities, commandline, mdtypes, …)
│   └── testutils/            # Shared Google Test helpers (refdata, TestFileManager, …)
├── share/                    # Data files: atomprops.csv, alexandria.csv, force field XMLs, …
├── cmake/                    # CMake modules inherited from GROMACS
├── docs/                     # Sphinx documentation source
├── scripts/                  # Shell/Python helper scripts (ACTRC env setup, …)
├── tests/                    # Top-level integration tests
├── examples/                 # Example input files
└── .github/workflows/        # CI: act-test.yml (Linux), act-test-macos.yml (macOS)
```

---

## Build System

ACT uses **CMake** (≥ 3.13). All ACT-specific code is under `src/act/`; GROMACS utilities are under `src/gromacs/`.

### Typical build (Linux, with MPI + double precision)

```bash
mkdir build && cd build
export MPICXX=/usr/bin/g++
export MPICC=/usr/bin/gcc
export CXXFLAGS='-std=c++20'
cmake -DGMX_DOUBLE=ON \
      -DCMAKE_PREFIX_PATH=/path/to/conda/env \
      -DBUILD_SHARED_LIBS=OFF \
      -DGMX_OPENMP=OFF \
      -DGMX_MPI=ON \
      -DCMAKE_INSTALL_PREFIX=../tools \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=${MPICXX} \
      -DCMAKE_C_COMPILER=${MPICC} \
      ..
make -j4 install tests
```

### Key CMake options

| Option | Typical value | Purpose |
|---|---|---|
| `GMX_DOUBLE` | `ON` | Use double precision (required for accuracy) |
| `GMX_MPI` | `ON` | Enable MPI parallelism |
| `GMX_OPENMP` | `OFF` | OpenMP (usually off) |
| `BUILD_SHARED_LIBS` | `OFF` | Static linking |
| `CMAKE_PREFIX_PATH` | conda env path | Finds RDKit, Boost, Eigen from conda |

### Required dependencies (via conda)

```
conda-forge::librdkit-dev=2025.09.4
conda-forge::libboost-devel
eigen
```

System packages: `gcc`, `cmake`, `libopenmpi-dev`, `openmpi-bin`, `libxml2`.

---

## Running Tests

After building, source the environment and run the per-module test binaries:

```bash
source ../tools/bin/ACTRC        # sets PATH and other env vars
cd build
./bin/basics-test
./bin/act-forces-test
./bin/alexandria-test
./bin/act-utility-test
./bin/coulombIntegrals-test
./bin/fileio-test
./bin/forcefield-test
./bin/import-test
./bin/molprop-test
./bin/properties-test
./bin/qgen-test
./bin/sobol-test
./bin/statistics-test
```

Tests use **Google Test** with the GROMACS reference-data framework (`testutils/refdata.h`). Reference XML files live in `tests/refdata/` subdirectories (e.g., `src/act/alexandria/tests/refdata/`). When adding a new test that generates output, run it once with `--gtest_filter=...` to create the reference XML, then commit both test and reference file.

---

## Code Conventions

### Language and standard
- **C++20** (`-std=c++20`).
- All ACT code lives in `namespace alexandria`.
- GROMACS utility code lives in `namespace gmx`.

### File headers
Every source and header file must start with the standard ACT copyright/license block (GPL v2+) followed by a Doxygen `/*! \internal \brief ... \author ... */` comment. See any existing `.cpp`/`.h` for the exact template.

### Include order (within a `.cpp`)
1. `"actpre.h"` (first, always)
2. The matching header (for `.cpp` files)
3. C/C++ standard headers
4. GROMACS headers (`gromacs/...`)
5. ACT headers (`act/...`)

### Naming conventions
- Classes: `PascalCase` (e.g., `ForceField`, `MsgHandler`, `ActAtom`)
- Private members: `trailingUnderscore_` (e.g., `value_`, `filename_`)
- Methods: `camelCase` (e.g., `chargeGenerationAlgorithm()`)
- Enums: `enum class` with `PascalCase` values
- Test classes: `<Topic>Test` (e.g., `ForceComputerTest`)

### Error handling
- Use `MsgHandler` (from `act/basics/msg_handler.h`) for all user-visible messages and errors.
  - `msghandler->fatal(ACTMessage::..., "description")` — fatal error (throws)
  - `msghandler->msg(ACTStatus::Warning, ACTMessage::..., "description")` — non-fatal warning
  - `msghandler->ok()` — check if no errors have been encountered
- For programmer errors, use GROMACS exception macros: `GMX_THROW(gmx::InvalidInputError(...))` or `GMX_THROW(gmx::InternalError(...))`.
- Guard debug-only output: `if (msghandler->debug()) { ... }` or use `msghandler->writeDebug(s)`.

### Verbosity levels (`ACTStatus` enum)
`Fatal < Error < Warning < Info < Verbose < Debug`

### Doxygen documentation
Use `/*! \brief ... */` for class/method docs and `//!` for brief inline member docs. `\param[in]`, `\param[out]`, `\return` for parameter documentation.

### clang-tidy
The project uses a `.clang-tidy` config at `src/.clang-tidy` (subset of bugprone, readability, modernize, cppcoreguidelines, google checks). Respect these when writing new code.

---

## Key Abstractions

### `ForceField` (`src/act/forcefield/forcefield.h`)
Holds all force field parameters. Loaded from XML with `readForceField()`. Contains `ForceFieldParameterList` maps keyed by `InteractionType`.

### `ForceFieldParameter` (`src/act/forcefield/forcefield_parameter.h`)
A single typed parameter with value, uncertainty, min/max, `Mutability` (Fixed/Bounded/Free), and unit string. Use `GMX_THROW(gmx::InvalidInputError(...))` when a value is out of bounds.

### `InteractionType` enum (`src/act/basics/interactiontype.h`)
Central enum for all interaction types: `BONDS`, `ANGLES`, `PROPER_DIHEDRALS`, `VDW`, `ELECTROSTATICS`, `POLARIZATION`, `VSITE1`…`VSITE4`, etc.

### `ActParticle` enum (`src/act/basics/act_particle.h`)
Particle types: `Atom`, `Vsite`, `Shell`, `SigmaHole`.

### `MolProp` (`src/act/molprop/molprop.h`)
Stores molecular properties from experiment or QM. Serialised to/from XML (`molprop_xml.h`) and SQLite (`molprop_sqlite3.h`).

### `ACTMol` (`src/act/alexandria/actmol.h`)
High-level molecule class. Combines `MolProp` with a `Topology`, `FragmentHandler` (for charge generation), and force computation capabilities.

### `Topology` (`src/act/alexandria/topology.h`)
Holds `ActAtom` list, bonded/non-bonded parameter lists (as `ForceFieldParameterList`), shells, virtual sites. Built from `MolProp` + `ForceField`.

### `FragmentHandler` (`src/act/alexandria/fragmenthandler.cpp`)
Manages charge generation for molecular fragments. Bond indices are 0-based within each fragment (excluding shell particles). Use `nonFixed_[b.aI()]` to map pre-build core positions to post-build atom indices when shells are interleaved.

### `ForceComputer` (`src/act/forces/forcecomputer.h`)
Computes energies and forces. `compute()` / `computeOnce()` expect the `forces` vector to be pre-sized to match `topology->atoms().size()`.

### `QgenAcm` / `QgenResp` (`src/act/qgen/`)
Charge generation via the ACM (Alexandria Charge Model) or RESP (Restrained Electrostatic Potential) methods.

### `MsgHandler` (`src/act/basics/msg_handler.h`)
Central logging/error-reporting object. Always pass as `MsgHandler *msghandler` pointer. Check `msghandler->ok()` after operations that may set warnings/errors.

### `CommunicationRecord` (`src/act/utility/communicationrecord.h`)
Thin wrapper around MPI communicators. Used throughout for parallel runs.

### Genetic Algorithm (`src/act/ga/`)
Generic GA framework. Specialised for force-field training via `AcmGa`, `AcmFitnessComputer`, `AcmIndividual`, `AcmInitializer` in `src/act/alexandria/`.

### `JsonTree` (`src/act/utility/jsontree.h`)
Simple JSON tree. All values are serialised as JSON strings (quoted). Do not expect typed values when parsing.

---

## CLI Modules (the `alexandria` executable)

Sub-commands registered in `src/act/alexandria/alex_modules.cpp`:

| Command | Purpose |
|---|---|
| `gentop` | Generate topology from structure / Gaussian output |
| `simulate` | Run MD simulation |
| `min_complex` | Energy scan inputs |
| `b2` | Second virial coefficient vs. temperature |
| `nma` | Normal mode analysis / thermochemistry |
| `train_ff` | Train force field parameters (GA / Bayesian) |
| `geometry_ff` | Derive bond/angle/dihedral distributions |
| `analyze` | Analyze molecular / FF properties, produce LaTeX tables |
| `edit_ff` | Manipulate / compare force field XML files |
| `gen_ff` | Generate a force field file from specification |
| `edit_mp` | Merge / edit molecular property files |
| `merge_ff` | Merge force field files |

---

## Data Files (`share/`)

- `alexandria.csv` — atom types
- `atomprops.csv` — element properties
- `atom_bond.xml` — bond parameters
- `atomization-energies.csv` — atomisation energies
- `symmetric_charges.csv` — symmetry constraints for charges
- `vsite*.csv` — virtual site configurations
- Force field XML files are loaded at runtime via `readForceField(filename, &msghandler)`.

---

## CI / Workflows

Two GitHub Actions workflows (`.github/workflows/`):
- `act-test.yml` — Ubuntu latest (gcc + OpenMPI)
- `act-test-macos.yml` — macOS latest (Homebrew gcc + Open MPI)

Both: checkout → conda (RDKit, Boost, Eigen) → cmake build → `make install tests` → run all test binaries.

---

## Common Pitfalls

1. **`actpre.h` must be first** in every `.cpp`. Omitting it causes cryptic include-order errors.
2. **Double precision** (`GMX_DOUBLE=ON`) is required. Single-precision builds are not supported.
3. **MPI is required** (`GMX_MPI=ON`). Serial builds are not tested.
4. **Forces vector sizing**: `ForceComputer::compute` does not resize the forces vector; the caller must size it to `topology->atoms().size()` before calling.
5. **Fragment bond indices** are 0-based within fragments and exclude shell particles. Use `nonFixed_[b.aI()]` for post-build atom positions.
6. **JsonTree** always writes values as strings; do not expect numeric JSON types on read-back.
7. **Reference data tests**: after adding a new Google Test that uses `TestReferenceChecker`, run the test once in update mode to generate the `.xml` reference file, then commit it.
8. **ACTRC sourcing**: tests that depend on shared data files require the `ACTRC` environment script to be sourced first.
