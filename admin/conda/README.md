# Conda packages of ACT

`recipe.yaml` builds two conda packages with
[rattler-build](https://rattler.build): `alexandria` (programs, libraries and
data) and `alexandria-python` (the Python package `act`, see below). The plan
and the decisions behind it are in `docs/conda-packaging-plan.md`.

| Platform    | Explicit SIMD | GPU                       | Variants      |
|-------------|---------------|---------------------------|---------------|
| `osx-arm64` | ARM NEON      | none                      | mpich, openmpi |
| `linux-64`  | AVX2 (whole build for x86-64-v3) | SYCL with AdaptiveCpp, CUDA backend | mpich, openmpi |
| `linux-aarch64` | ARM NEON | SYCL with AdaptiveCpp, CUDA backend | mpich, openmpi |

## Building locally

Install rattler-build in an environment of its own, then build from the last
commit (untracked files, such as build directories, are not used):

    conda create -n rattler -c conda-forge rattler-build
    conda activate rattler
    admin/conda/build-local.sh output-conda --variant mpi=mpich

The build runs the unit tests (`ctest`) and then the package tests in
`recipe.yaml`. Rebuilding the same version into an output directory that
already has the packages can fail with a hash mismatch; remove the output
directory first.

To try the packages:

    conda create -n alex-test -c file://$PWD/output-conda -c conda-forge alexandria-python
    conda activate alex-test
    alexandria -version
    python -c "import act; print(act.classes_available(), act.data_directory())"

## The Python package

Both packages are built in single precision only, with the GPU support
where there is one (the build `sycl_cuda` on Linux, which also runs on the
processor). Double precision, which the development and the tests of ACT
use, needs a build from the source code.

`alexandria-python` has the Python package `act`: the alexandria program
from Python (`act.run`, `act.commands`) and the classes of ACT
(`act.ForceField`, `act.Molecule` and so on) in the compiled module
`act._core`, see the chapter on Python in the manual. It requires the
`alexandria` package it was built with, for the program and the data files.

- The module is compiled for each Python version in the key `act_python` of
  `variants.yaml`, for each MPI library, since it uses MPI like the
  libraries of ACT. The build string says which: `py313_mpi_mpich_...`.
- The code of ACT is in the module itself (static libraries, built by
  `build-python.sh` with `-DACT_GPU=NONE`), so that importing it loads
  nothing of a GPU. The classes compute on the processor, also where the
  programs of `alexandria` use a GPU.
- Each Python version compiles ACT once more, a few minutes each. To build
  for one version only:

      admin/conda/build-local.sh output-conda --variant mpi=mpich --variant act_python=3.13

- The build runs the tests of the classes (`python/tests/test_core.py`) with
  the installed module and the program of `alexandria`; the package test
  (`check-python-package.py`) checks that the installed package finds the
  program and the data files of its environment.

## Uploading

Test builds go to the channel `alexandriachemistry` on anaconda.org with the
label `dev`, releases with the label `main`:

    anaconda -t <token> upload -u alexandriachemistry -l dev output-conda/*/alexandria-*.conda

The workflow `.github/workflows/act-conda.yml` builds `act-utilities`,
`alexandria` and `alexandria-python` on macOS (Apple silicon) and Linux
(x86-64 and 64-bit ARM), and uploads them
with the secret `ANACONDA_API_TOKEN`:

- a pull request that changes the recipe builds and tests the packages,
  `alexandria-python` for one Python version only, without uploading them;
- a tag `v*` uploads them with the label `main`;
- a run by hand (Actions, Conda packages, Run workflow) uploads them if
  *upload* is ticked, with the label given, `dev` by default.

The packages of each run are also kept as artifacts of the run.
Releases on conda-forge are built by the feedstock `alexandria-feedstock`
from the archive of the release tag.
