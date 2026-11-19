# Plan: a Python API and an HTML5 user interface for ACT

Status: in progress. The decisions of 2026-11-19 on the package name, the
repository, the default force field, Slurm and the login are in section 1.
Implemented so far, on the branch `python-api`: step A1, the first three
items of step A2 and the first version of step A4 (section 5).

This plan answers [issue 59](https://github.com/AlexandriaChemistry/ACT/issues/59)
("Python binder for ACT") and the specification for a new user interface. It
describes

- a Python package `act` that opens up ACT to users in their own scripts, and
- a user interface that runs in a browser, on a desktop and on a server, and
  that uses only that Python package to talk to ACT.

A clickable mockup of the interface is in `docs/act-gui-mockup.html` (open it
in a browser). It shows real ACT results for acetone and example numbers for
the molprop screens.

## 1. Decisions taken

| Question | Decision | Why |
|---|---|---|
| Binding tool | pybind11 | Portable, open source, header-only, on conda-forge. Tested here with pybind11 3.1 and Python 3.14 (section 3). |
| Depth of the API | Two layers. Layer 1 runs every `alexandria` subcommand from Python and returns its JSON results. Layer 2 exposes C++ classes (force field, molecule, energies, minimization, normal modes, simulation). | Layer 1 gives a complete API quickly and is safe: a fatal error in ACT becomes a Python exception. Layer 2 is what issue 59 asks for and makes the interface fast for small molecules. |
| Front end | Plain HTML, CSS and JavaScript modules, no build step. | No Node.js toolchain in the repository, the CI or the conda recipe. People who know Python and some JavaScript can maintain it. |
| Back end | Python with FastAPI, on top of the `act` package. The same code serves the desktop and the server. | One code path to test. The interface cannot do anything a Python user cannot do. |
| Molecular viewer | 3Dmol.js (BSD-3-Clause), shipped inside the package | An existing tool, as the specification asks. It rotates, zooms, colours by a property such as charge, shows volumetric data, animates vibrations and trajectories, and reports clicks on atoms. All of this works in the mockup. |
| Plots | uPlot (MIT) in the browser for live plots; plotXVG on the server for PDF files | The plot on screen updates while a simulation runs. The PDF is what users already get from plotXVG, so screen and paper match. |
| Server | In two stages: first for known users with a login, later a public service without accounts | The public stage needs much stricter isolation. The design of stage 1 prepares for it (section 8). |
| First deliverable | This plan and the design mockup | To agree on scope and design before code is written. |
| Name of the Python package | `act`. The class `ACT` of the module `act.py` in ACT-utilities stays available for one release as `act.legacy.ACT`. | `import act` is what issue 59 uses. Scripts that use the old class keep working for one release. |
| Repository of the interface | A new repository `ACT-GUI` with its own noarch conda package `act-gui` | ACT-utilities does not get the dependencies of a web server. |
| Default force field | `ACS-pg1.xml` from `src/act/testdata/forcefields`, until the new default force field is ready | A new default is under development. The name is one setting (`default_forcefield`), so that replacing it is a change of one line. The interface labels it as a test force field. |
| Long jobs | The server can send them to Slurm on csb for the local users (step B11) | `train_ff -nooptimize` on large molprop files and long simulations do not belong on the web server. |
| Login for stage 1 | The university single sign-on, through OpenID Connect | The same identity can later own the jobs on csb. Single sign-on for csb itself is a separate issue (section 11). |
| Account for jobs on csb | The service account `act-user`, for now | Jobs under the user's own csb account need single sign-on on csb. The name of the account is a setting of the server. |
| Where the server runs | On a machine connected to the csb cluster: it can call `sbatch` and shares a file system with the compute nodes | The Slurm back end needs no SSH connection and can read the progress of a job from its directory. |

## 2. What exists today

These findings decide the order of the work. File names are relative to the
ACT repository.

**The subcommands.** The twelve subcommands (`gentop`, `simulate`, `minimize`,
`min_complex`, `b2`, `nma`, `train_ff`, `geometry_ff`, `edit_ff`, `gen_ff`,
`edit_mp`, `merge_ff`) are registered in `src/act/alexandria/alex_modules.cpp`.
Each is one function `int tool(int argc, char *argv[])` that parses the command
line and does the work. None has a class that can be called without `argv`.

**JSON output.** `simulate`, `minimize`, `nma`, `b2` and `train_ff` can write
JSON through the class `JsonTree` (`src/act/utility/jsontree.h`), with
nlohmann/json underneath. The other tools write none. Details that matter for
the interface:

| Tool | What is there | What is missing |
|---|---|---|
| `minimize` | Energies before and after, per term | The RMSD is only in the log file. `-json` is a flag that writes the fixed name `minimize.json`. |
| `nma` | Frequencies, intensities, thermochemistry at 0 K and 298.15 K; the IR spectrum as an xvg file | The eigenvectors are not written anywhere, so vibrations cannot be shown. The root key of the JSON is `simulate` instead of `nma`. Fixed file name `nma.json`. |
| `simulate` | Trajectory (PDB, one model per frame), energies (xvg), JSON summary | The files are not flushed during the run, so a plot cannot follow the simulation. |
| `train_ff -nooptimize` | Statistics per data set and property (N, RMSD, MSE, MAE, slope, R), per-molecule data, correlation data as xvg files | Nothing; this is enough for the correlation plots. |
| `gentop` | Topology, OpenMM XML, PDB, cube files of the electrostatic potential | No JSON with atoms, charges, shells and virtual sites. |
| `gen_ff`, `geometry_ff`, `edit_ff`, `edit_mp`, `merge_ff` | Log files only; `gen_ff -vs` reads JSON definitions of virtual sites | No machine-readable results, for instance of `edit_ff -validate` and `edit_mp -validate`. |

All leaf values in the JSON are strings (`"value": "-3596.84"`), also the
numbers.

**Reading molecules.** `CompoundReader` (`src/act/import/compound_reader.cpp`)
reads sdf, pdb and xyz files with RDKit. There is no SMILES input anywhere in
ACT, although the RDKit SMILES parser is already linked. Its file name and
other settings can only be set from the command line.

**Classes that already work without a command line.** `ForceField` with
`readForceField` and `writeForceField`, `ACTMol`, `ForceComputer`,
`MolHandler` (`minimizeCoordinates`, `nma`, `computeHessian`, `simulate`,
`coordinateRmsd`), `ThermoChemistry`, `validateForceField`,
`validateMolProps`, `writeOpenMM`, `MolPropRead` and `MolPropWrite`.

**Errors.** Three mechanisms are mixed: C++ exceptions (about 250 `GMX_THROW`),
a status flag in `MsgHandler` that the caller must check, and termination of
the process (`gmx_fatal`, 42 calls, and `GMX_RELEASE_ASSERT`, 126 calls, in
`src/act`). The last one ends in `MPI_Abort` or `_Exit`.

**MPI.** `CommunicationRecord` needs an initialised MPI library, and every
tool makes one. The classes listed above do not need it.

**Other state.** About 170 `printf` and `fprintf(stderr)` calls in `src/act`
write directly to the terminal. The data directory is found through the
program context or the environment variable `ACTDATA`. The choice of GPU is
global to the process.

**The XML readers** call `xmlParseFile` with the default options
(`forcefield_xml.cpp`, `molprop_xml.cpp`, `import_utils.cpp`). On a server,
where the files come from strangers, that is not strict enough (section 8).

**Python today.** ACT-utilities has `act.py`, a class that builds command
lines and runs them with `os.system`, and `actgui`, a Tk program of 1600
lines that makes the input for `gen_ff`. Both are replaced by this plan.
`act_openmm.py` (class `ActOpenMMSim`) and plotXVG (`plotxvg.plot(files,
save="plot.pdf")`) are used as they are.

**Build.** C++20. `libactLib` and `libactgromacs` are shared libraries by
default. No headers are installed, which is fine: the bindings are built
inside the ACT build. Doxygen comments use `/*! \brief ... */`. `make tidy`
runs `admin/clang_tidy_fix.py` over `src/act`.

## 3. Feasibility test (done)

Before writing this plan I tested the riskiest assumption: that a pybind11
module can use the ACT library as it is today. A module of 150 lines, linked
against `libactLib.a` of `build_Release_SINGLE`, was called from Python 3.14.
The source is not in the repository; step A5 starts from it.

| Test | Result |
|---|---|
| Read a force field and return a summary to Python | Works: `ACS-pg1.xml`, 77 particle types, 10 interaction types |
| Read `acetone.sdf`, make the topology, generate charges, compute the energy, minimize, compute the RMSD | Works. Same numbers as `alexandria minimize`: RMSD 0.0192 nm, energy from -3599.17 to -3764.84 kJ/mol |
| Return the Hessian and get the normal modes with numpy | Works: 24 modes, which match the 24 frequencies of `alexandria nma` one to one |
| Time per call | 0.10 s inside Python, against 0.37 s for the `alexandria minimize` program |
| Start-up without `MPI_Init` | Works. Layer 2 does not need MPI, so `import act` is safe on login nodes where MPI cannot start. |
| A force field file that does not exist | The whole Python process is killed (`gmx_fatal`, then `MPI_Abort`). This confirms that the error handling must change before layer 2 is usable (step A4). |

The mockup uses the output of this test and of `alexandria minimize`, `nma`
and `simulate` for acetone.

## 4. Architecture

```
 Browser (desktop or remote)                     User's own Python script
   HTML5 front end, plain JS modules                        |
   3Dmol.js viewer, uPlot plots                             |
        |  HTTPS: JSON, files, event stream                 |
        v                                                   |
 act-gui server (Python, FastAPI)                           |
   routes, login, workspaces, job queue                     |
        |  starts one worker process per job                |
        v                                                   v
 +----------------------------------------------------------------+
 | Python package "act"                                           |
 |   Layer 1  act.commands   one function per subcommand,         |
 |            runs "alexandria <tool>" as a child process,        |
 |            returns the JSON results and the output files       |
 |   Layer 2  act.ForceField, act.Molecule, ...  (pybind11)       |
 +----------------------------------------------------------------+
        |                                   |
        v                                   v
   alexandria executable             libactLib (C++), in process
```

Three rules follow from the specification and from section 2.

1. **The interface only uses the Python package.** Every button in the
   interface corresponds to a call a user can make in a script. The interface
   shows that call (the mockup shows the command line for `gen_ff`).
2. **The web server never runs ACT code in its own process.** A job runs in a
   worker process. If ACT stops with a fatal error, or a file reader crashes
   on a broken file, the job fails and the server continues.
3. **Results travel as JSON and files.** The subcommands write JSON; the
   Python layer turns it into objects with numbers and units; the server
   passes it to the browser.

### Where the code lives

| Part | Place | Package |
|---|---|---|
| pybind11 bindings (C++) | ACT, new directory `src/act/python` | conda package `alexandria`, option `ACT_PYTHON=ON` |
| Python code of the API (layer 1, result classes) | ACT, `python/act` | same package, in `site-packages/act` |
| Server and front end | New repository `ACT-GUI` | new noarch conda package `act-gui`, depends on `alexandria`, `fastapi`, `uvicorn`, `plotxvg` |
| `actgui` in ACT-utilities | removed when `act-gui` is released (step B12) | |
| `act.py` in ACT-utilities | moved to `act.legacy` in step A4, removed one release later | |

The API has to live in ACT because it must match the C++ code of the same
commit and is tested in the same CI. The interface gets its own repository,
like plotXVG, because its tools and tests are different.

The new package takes the name `act`. Step A4 moves the class `ACT` of
ACT-utilities to `act.legacy` unchanged, with a deprecation warning, and
ACT-utilities stops installing its own `act.py` in the same release, so that
the two packages never install the same file.

## 5. Part A: the Python API

Every step starts with its tests. C++ features get a gtest with reference
data in the directory of the tool; Python code gets pytest tests that run in
the CI of ACT. Every step ends with `make tidy` and Doxygen comments.

### Layer 1: the subcommands from Python

**A1. A description of the options of every subcommand, as JSON.** *(done)*
New in ACT: `alexandria help -export json` writes the file
`alexandria-options.json` with, for every subcommand, its description, its
help text and its options (name, type, default, allowed values, help text,
whether it is hidden or required and, for files, whether the file is read or
written and which extensions it may have). It is generated from the same
option objects that print the help and parse the command line, so it cannot
differ from what the program accepts, and nothing is listed twice. Python
uses it to check keyword arguments before starting a program, and the
interface uses it to build forms and help texts.
The export sits next to the existing exports for reStructuredText and shell
completion (`src/gromacs/commandline/optionsjsonwriter.cpp`), and gives all
subcommands in one call instead of the `-options json` per tool that the
first version of this plan proposed.
*Test first:* `optionsJsonTest.cpp` checks that the output is valid JSON,
contains every subcommand and, for each, exactly the options of
`-h -hidden`, with the hidden ones marked.

**A2. Complete the JSON output of the tools.** Each item is a small change in
ACT with its test written first.

| Change | Needed for | Status |
|---|---|---|
| `JsonTree` writes numbers as JSON numbers; also booleans and arrays | Everything; removes conversions in Python and JavaScript | done |
| `minimize`: RMSD and convergence status in the JSON; `-json` takes a file name; the molecule (particles, charges, bonds) and the coordinates before and after, superimposed | Minimization panel, viewer and charge display | done |
| `nma`: root key `nma`; for each mode the displacement of each particle and the reduced mass; the molecule and the minimized coordinates; `-json` takes a file name | Animation of vibrations | done |
| `simulate`: flush the energy and trajectory files at every output step, and write a progress line (step, time, energies) to a JSON-lines file | Energy plot that follows the simulation | to do |
| `gentop`: JSON with particles (atoms, shells, virtual sites), types, charges, bonds and coordinates, with `moleculeToJson` as in `minimize` | Viewer and charge display | to do |
| `edit_ff`: parameters as JSON (now `-dump` writes a .dat file) and the result of `-validate` as JSON | Force field editor | to do |
| `edit_mp`: the result of `-validate` as JSON, and a summary of the content (compounds, dimers, number of conformations, properties) | Molprop screens | to do |
| `gen_ff`, `geometry_ff`: a JSON summary of the force field that was made | Force field generator | to do |

The existing tests (`simulateTest`, `nmaTest`, `b2MpiTest`, `trainffTestUtil`)
already read JSON, so each change extends a test that is there.

**A3. SMILES input.** New in ACT: `-smiles` next to `-f` in `CompoundReader`.
RDKit parses the string, adds hydrogens, makes 3D coordinates (ETKDG) and
optimizes them roughly; then the molecule follows the same path as a file.
This needs two more RDKit components in CMake (`DistGeomHelpers`,
`ForceFieldHelpers`). It also gives `CompoundReader` setters, so that it can
be used without a command line.
*Test first:* SMILES of ten molecules of the test set give the same formula,
bonds and atom types as their sdf files; invalid SMILES give an error
message, not a crash.

**A4. Layer 1 in Python.** Package `act` in the directory `python/act` of
ACT, module `act.commands`. *The first version is done; what remains is
listed below.* What works now:

```python
import act

ff  = "ACS-pg1.xml"
res = act.commands.minimize(ff=ff, f="acetone.sdf")
print(res.data["rmsd"], res.data["Energies after"]["EPOT"])  # value and unit
print(res.file("c"))                                         # minimized structure

modes = act.commands.nma(ff=ff, f="acetone.sdf", ir=True).data["modes"]
print(modes[0]["frequency"], modes[0]["displacements"])

stats = act.commands.train_ff(ff=ff, mp="dimers.xml", optimize=False,
                              processes=8)                   # runs under mpirun
help(act.commands.train_ff)                                  # the options
```

- One function per subcommand, generated from the option description of A1,
  with docstrings from the help texts.
- The program runs as a child process with an argument list, never through a
  shell, in a directory of its own. The function returns a result object
  with the parsed JSON, the paths of the output files and the log.
- A non-zero exit status becomes an exception `act.ActError` that carries the
  error text from the log.
- A callback or an iterator gives progress during long runs (`simulate`,
  `train_ff`), read from the progress file of A2.
- `processes=N` starts the tool with the MPI launcher of the installation.

- Values cannot pose as options: a string that starts with a dash is
  refused, unless it is a negative number. Input files are checked and made
  absolute before the program starts.

*Tests:* `python/tests`, run by ctest as `PythonInterfaceTests`. They use the
standard `unittest` module, so that they need no extra package; pytest runs
them as well. Most of them use a stand-in for the program and need no build
of ACT; the others run `minimize` and `nma` with the real program.

*Still to do in A4:* result classes with named fields (`res.rmsd`,
`modes.frequencies`) once the JSON of all tools has its final form; progress
during long runs; tests for the remaining subcommands; the installation of
the package (section 11, question 1).

After A1 to A4 the interface can be built. Layer 2 can proceed in parallel.

### Layer 2: classes

**A5. Errors that Python can catch.** This is the condition for everything
else in layer 2.

1. One switch in `src/gromacs/utility/fatalerror.cpp`: in "library mode"
   `gmx_fatal` and the assert handler throw a C++ exception instead of
   ending the process. The module turns the switch on when it is imported;
   the `alexandria` program keeps its behaviour. This one change covers all
   170 places.
2. Then, file by file along the paths the API uses (force field XML, molprop
   XML, import, topology, charge generation), replace `gmx_fatal` by
   exceptions of a specific type, so that Python gets `FileNotFoundError`,
   `act.ForceFieldError`, `act.TopologyError` and so on.
3. Functions that only set the status of `MsgHandler` get a check in the
   binding that raises if the status is not ok.

*Test first:* gtest and pytest that a missing force field, a broken XML file
and an unknown atom type raise the expected exception and leave the process
alive.

**A6. Build and start-up.** CMake option `ACT_PYTHON` (default OFF), finds
pybind11, builds the module `act._core` against `libactLib`, and sets
position-independent code for static builds. At import the module sets the
data directory from the location of the package, does not start MPI (see
section 3), and connects `MsgHandler` to the Python `logging` module, so that
messages appear where Python users expect them. The remaining `printf`
calls on the paths of the API move to `MsgHandler` as they are met.

**A7. The classes.** In this order, each with pytest tests and docstrings:

| Python | On top of | Main methods |
|---|---|---|
| `act.ForceField` | `ForceField`, `readForceField`, `writeForceField`, `validateForceField` | `read`, `write`, `particle_types`, `interactions`, `parameters(interaction)`, `parameter(id, name)` with value, bounds, unit and mutability, `set_value`, `validate()` |
| `act.Molecule` | `ACTMol`, `importFile`, `Topology` | `from_file`, `from_smiles`, `atoms`, `bonds`, `coordinates` (numpy), `charges`, `generate_charges(algorithm)`, `fragments`, `to_sdf`, `to_pdb` |
| energies | `ForceComputer` | `mol.energy()` (dictionary per term), `mol.forces()`, `mol.interaction_energy()` for dimers |
| `mol.minimize()` | `MolHandler::minimizeCoordinates`, `coordinateRmsd` | returns coordinates, energies before and after, RMSD, status |
| `mol.normal_modes()` | `MolHandler::nma`, `ThermoChemistry` | frequencies, intensities, eigenvectors, thermochemistry, `spectrum(linewidth)` |
| `mol.simulate()` | `MolHandler::simulate` | settings as keyword arguments, a callback per output step, returns trajectory and energies as numpy arrays |
| `mol.write_openmm()` | `writeOpenMM` | the force field file for OpenMM |
| `act.MolProp`, `act.read_molprop` | `MolPropRead`, `MolPropWrite`, `validateMolProps` | iterate over compounds, experiments and energies; `validate()` |

This is the API of issue 59, with other names:

```python
import act

ff  = act.ForceField.read("ACS-pg1.xml")
mol = act.Molecule.from_smiles("CC(=O)C", ff)      # topology and charges
print(mol.charges, mol.energy()["EPOT"])
res = mol.minimize()
print(res.rmsd)
mol.write_openmm("acetone-openmm.xml")
```

Work in C++ that this needs, each with a test first:

- `SimulationConfigHandler` and `CompoundReader` get setters for what is now
  only reachable through command-line options.
- `MolHandler::nma` returns the eigenvectors (now it discards them).
- `MolHandler::simulate` gets a callback per output step and does not need
  the output environment of the command line.
- Long calculations release the Python lock and check for Ctrl-C between
  steps.
- Arrays go to numpy without a copy where the layout allows it.

**A8. Split the tools into "parse" and "run".** For `minimize`, `nma`,
`simulate`, `gentop`, `gen_ff`, `edit_ff` and `edit_mp`: the function of the
tool keeps the command line and fills an options structure; a second
function does the work from that structure. The existing regression tests
protect each split. After this, layer 1 and layer 2 share one
implementation, and the subcommands are easier to test.

**A9. Documentation and packaging.** A chapter "Using ACT from Python" in the
manual with the examples above, an API reference generated from the
docstrings, and the option `ACT_PYTHON=ON` in the conda recipe with an import
test and one small calculation as package test.

## 6. Part B: the user interface

### How it is built

- **Server** (`act-gui`, Python): FastAPI application. Routes for molecules,
  force fields, molprop files, jobs and result files. A job is one call of
  the `act` package in a worker process with its own directory. Progress
  goes to the browser as a server-sent event stream.
- **Front end**: one page, JavaScript modules loaded directly by the browser,
  one module per workspace, small shared modules for the viewer, the plots,
  the tables and the forms. 3Dmol.js and uPlot are copied into the package,
  so the interface works without internet access and a strict content
  security policy is possible.
- **Desktop**: the command `act-gui` starts the server on `127.0.0.1` with a
  random port and a one-time token, and opens the browser, as Jupyter does.
  Files are in a directory the user chooses. No login.
- **Server**: `act-gui serve --config server.toml` behind a reverse proxy
  with TLS. Login, a workspace per user, quotas and a job queue (section 8).

### Steps

Each step ends with tests of the routes (pytest), of the JavaScript modules
(the test runner built into Node.js, which needs no build step) and one
end-to-end test in a real browser (Playwright) that does what a user would
do.

**B0. Skeleton and design system.** The server, the page with the three
workspaces, the tabs, the buttons, the tables and the plot and viewer
modules, taken from the mockup. Light and dark theme. The end-to-end test
set-up in CI.

**B1. Job engine.** Start, follow, cancel and list jobs; a directory per job;
results as JSON; download of output files; time and size limits. The engine
talks to a back end with three operations (submit, status, cancel). The
first back end starts child processes on the same machine; step B11 adds
Slurm without changes to the routes or the pages. Tested with `alexandria`
on the ACT test data and with a tool that fails on purpose.

**B2. Molecule input, viewer and charges.** Upload of sdf, pdb and xyz files
or a SMILES string (A3); choice of the default force field or one the user
uploads; choice of the charge method. The viewer shows the molecule with
rotation, zoom and translation by mouse. Charge distributions are shown in
four ways: atoms coloured by net charge, labels with the charges, shells and
virtual sites as separate small particles, and the electrostatic potential
on the molecular surface from the cube file that `gentop -pot` writes.

**B3. Minimization.** Runs `minimize`. Shows the minimized structure on top
of the original one, the RMSD and the energy per term before and after.

**B4. Normal modes.** Runs `nma`. Shows the spectrum, a table of modes, the
vibration of the selected mode in the viewer (from the eigenvectors of A2)
and the table of thermochemistry. The spectrum can be saved as PDF through
plotXVG.

**B5. Simulation.** Runs `simulate` with a limit on the number of steps.
The energy plot grows while the job runs; the trajectory plays in the viewer
with a slider; energy plots can be saved as PDF through plotXVG.

**B6. OpenMM package.** Runs `gentop -openmm` and collects the coordinates,
the force field file, a `simulation.dat` for `act_openmm.py`, a run script
and a README with the commands into one zip file. Tested by running ten
steps with OpenMM in CI.

**B7. Force field generation.** The functions of `actgui`: potentials,
charge distribution, polarizability, combination rules, selection of atom
types, the predefined virtual sites; runs `gen_ff` and then, with a molprop
file, `geometry_ff`. Shows the summary and offers the XML file.

**B8. Force field editing and validation.** A table of the parameters with
value, bounds and mutability, filters per interaction and atom type, changes
through `edit_ff`; comparison of two force fields; `edit_ff -validate` with
the result as a list of findings.

**B9. Virtual sites by point and click.** The user loads a molecule, clicks
two or three atoms, chooses the kind of site and moves it with a slider. The
interface writes the definition in the JSON format of `gen_ff -vs` and adds
it to the force field. Layer 2 gives the position of the site at once.

**B10. Molprop files.** Upload and `edit_mp -validate` with the findings;
evaluation of a force field with `train_ff -nooptimize`, correlation plots
per property and a table of statistics; selection of a compound or dimer and
a histogram of its energies; all plots can be printed through plotXVG.

**B11. Server, stage 1.** Login with the university single sign-on,
workspaces, quotas, queue, isolation of the workers (section 8); an
administration page; documentation for installing it.

Long jobs go to Slurm on csb through a second back end of the job engine:

- *Which jobs.* Interactive calculations (charges, minimize, nma, short
  simulations) stay on the server. A job goes to Slurm when its estimate is
  above a limit in the configuration: `train_ff -nooptimize` on a large
  molprop file, a long simulation. The page shows where the job runs and
  its place in the queue.
- *The batch script.* The server writes it from a fixed template. Partition,
  account and the maximum wall time come from the configuration of the
  server. The user chooses the number of cores and the time from a short
  list. No text of the user enters the script; `sbatch` is started with an
  argument list (section 8, measure 2).
- *Following a job.* The job writes its progress (the JSON lines of step A2)
  and its results in the job directory, which the server and the compute
  nodes share. The server reads these files. It asks Slurm for the state of
  all its jobs in one `squeue` call, at most once a minute, and never per
  job or per page view.
- *Whose job it is.* The aim is that a job runs under the csb account of the
  user who started it, so that accounting and file permissions are right.
  That needs single sign-on on csb, which is a separate issue. Until then
  the service account `act-user` submits the jobs, with the user in the
  `--comment` field and a quota per user on the server. The name of the
  account is a setting, and the back end keeps the submitting account in one
  place, so that the change to accounts per user stays small.
- *Where the server runs.* On a machine connected to the csb cluster, which
  can call `sbatch` and shares a file system with the compute nodes. The job
  directories are on that file system.
- *Tests.* The back end is tested against a small script that imitates
  `sbatch`, `squeue` and `scancel`, so that CI needs no cluster, and by hand
  on csb before the release.

The public server of stage 2 does not send jobs to Slurm.

**B12. Release.** Conda package `act-gui`, a manual chapter with screenshots,
removal of `actgui` from ACT-utilities.

**B13. Server, stage 2.** The public service (section 8).

## 7. Design

The mockup is the reference for the look. Its rules:

- **Three workspaces** in a top bar: Molecule, Force field, Molprop data.
  Inside a workspace, tabs for the tools. Every tool has its own link, so a
  page can be bookmarked or sent to a colleague.
- **Results get the room.** Inputs sit in one compact bar or a narrow
  column. The viewer and the plots take the rest. The button that starts a
  calculation is next to its results.
- **One accent colour**, the blue of the ACT logo, for the action that
  starts something and for the current selection. Colours for good, warning
  and error are separate and only used for status. Charges use the red to
  blue scale that chemists know.
- **Typography**: IBM Plex Sans for text, IBM Plex Mono for numbers, units,
  identifiers and file names, so that columns of numbers line up.
- **Motion** only when the interface changes: a panel fades in over 0.2 s,
  the underline slides to the selected tab, a progress bar grows. Nothing
  moves on its own, and motion is off for users who ask their system to
  reduce it.
- **Honest about what runs.** The interface shows the command or Python
  call behind a result, so that a result can be reproduced in a script.
- It works on a narrow screen, in light and dark theme, and with the
  keyboard.

The mockup was checked in a browser at desktop and phone width, in both
themes, for all tabs.

## 8. Security of the server

The danger is specific: the server runs a C++ program on files that
strangers upload. The measures, in order of importance:

1. **Isolation of every job.** Stage 1: the worker runs as a separate user
   without privileges, in the directory of the job, with limits on CPU
   time, memory, file size and number of processes, and a wall-time limit.
   Stage 2: in addition each job runs in a container or namespace
   (bubblewrap, nsjail or rootless podman) without network, with a
   read-only system and a temporary directory.
2. **No path from user text to a command.** Programs are started with an
   argument list, never a shell. Option names and values are checked against
   the description of step A1. There is no field for "extra options".
   File names on disk are made by the server; names from users are only
   labels.
3. **Strict reading of uploads.** Limits on size and number of atoms; the
   type is checked from the content. The XML readers move from
   `xmlParseFile` to `xmlReadFile` with network access and external
   entities switched off, and with a limit on the size after decompression.
   SMILES strings have a length limit and are parsed in the worker.
4. **Fuzz tests of the readers** (sdf, pdb, xyz, molprop XML, force field
   XML, virtual site JSON) with libFuzzer and the existing AddressSanitizer
   build, in CI. They find the crashes before a stranger does.
5. **The web application.** A content security policy without inline
   scripts and without other hosts; text from users and from files enters
   the page only as text; protection against request forgery; session
   cookies that scripts cannot read; downloads only through job numbers the
   user owns; rate limits.
6. **Python.** Only JSON for data from users: no pickle, no `eval`, no
   unsafe YAML. Dependencies are pinned and checked for known problems in CI.
7. **Login (stage 1).** The university single sign-on through OpenID
   Connect. A workspace and quotas per user. Only users on a list in the
   configuration (or in a group of the identity provider) may send jobs to
   Slurm.
8. **Jobs on the cluster (stage 1).** A job on a compute node reads the
   same uploaded files as a job on the server, so it gets the same limits:
   the batch script starts `alexandria` through the same wrapper that sets
   the limits on memory, file size and processes, and Slurm enforces the
   wall time and the cores. The service account `act-user` owns nothing but
   the job directories and may only submit to the partition in the
   configuration.
9. **Stage 2, public.** Sessions without account, strict quotas, no long
   jobs (`train_ff` only on small files, nothing on the cluster), removal
   of data after a fixed time, a privacy notice, and monitoring for abuse.

On the desktop the server only listens on the local machine and needs the
token, so other users and web pages cannot reach it.

## 9. Conventions for all new code

- **Tests first.** A feature in ACT starts with its failing test. Remember
  `make tests` before `ctest`.
- **Doxygen-style comments.** C++: `/*! \brief ... \param[in] ... \return */`
  as in the rest of ACT. Bindings: the same text as docstring, so that
  `help(act.Molecule)` works. Python: docstrings in Doxygen form
  (`"""! @brief ... @param ... @return ..."""`). JavaScript: JSDoc comments
  (`/** ... @param ... @returns ... */`), the nearest equivalent, since
  Doxygen does not read JavaScript.
- **`make tidy`** passes before every push; `src/act/python` is added to the
  paths it checks.
- **Semantic commit messages**, as in the rest of ACT.

## 10. Order of work and rough size

The sizes are estimates in working days for one developer, before review.
They are rough; the refactoring steps (A5, A8) are the least certain.

| Step | Work | Days |
|---|---|---|
| A1 | Option description as JSON | 2 |
| A2 | JSON output of the tools completed | 5 |
| A3 | SMILES input | 2 |
| A4 | Layer 1 in Python, tests | 4 |
| B0, B1 | Skeleton, design system, job engine | 6 |
| B2 to B6 | Molecule workspace | 10 |
| A5, A6 | Errors as exceptions, build of the module | 5 |
| A7 | Classes of layer 2 | 10 |
| B7 to B9 | Force field workspace | 9 |
| B10 | Molprop workspace | 5 |
| A8 | Tools split into parse and run | 6 |
| A9, B12 | Documentation, conda packages, release | 4 |
| B11 | Server, stage 1, with the Slurm back end | 12 |
| B13 | Server, stage 2 | 8 |

Milestones:

1. **API 0.1** (A1 to A4): every subcommand from Python with JSON results.
2. **Desktop alpha** (B0 to B6): the molecule workspace on the desktop.
3. **API 1.0** (A5 to A7): the classes of issue 59.
4. **Desktop 1.0** (B7 to B10, A8, A9, B12): all functions of the
   specification, conda package, `actgui` retired.
5. **Server for known users** (B11), with long jobs on csb.
6. **Public server** (B13).

## 11. Open questions

The questions on the name of the package, the repository, the default force
field, Slurm, the login, the account for jobs on csb and the place of the
server are answered; the answers are in section 1. What remains:

1. **Where the Python package is installed.** ACT-utilities installs its
   modules as loose files in `site-packages/act` and adds that directory to
   the path of Python with `act.pth`. The package `act` belongs in the same
   directory. `import act` then gives the new package, but the modules of
   the package would also be importable under their bare names (`commands`,
   `options`, `errors`), which can hide modules of other packages. Proposal:
   ACT-utilities moves its files to `site-packages/act_utilities`, with its
   `.pth` file, the check in the CMake files of ACT and ACTRC following. That
   is a change in ACT-utilities and in both conda recipes, so the package is
   not installed yet: it is used from the source tree and tested by ctest.
   The old class is reachable as `act.legacy.ACT` when ACT-utilities is
   installed, and `from act import ACT` still works with a warning.
2. **Single sign-on on csb (separate issue).** For jobs under the user's own
   csb account, the identity of the university login must map to a csb
   account, and the server must be allowed to submit for that user. Until
   that exists the jobs run as `act-user` (step B11).
3. **Precision.** The conda packages are single precision. Proposal: the
   Python module follows the build it belongs to, and the first release is
   single precision only.
4. **GPU.** Small molecules do not gain from the GPU. Proposal: the
   interactive calculations run on the CPU; `train_ff -nooptimize` on large
   molprop files may use the GPU as the program does now, also on csb.
5. **Which wavenumber direction for spectra?** The mockup plots increasing
   wavenumbers, as the xvg file of `nma` and plotXVG do. Spectroscopists
   often plot them decreasing.

None of these blocks the work on the steps A1 to A8 and B0 to B10. The
first has to be answered before a release (steps A9 and B12).
