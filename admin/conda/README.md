# Conda packages of ACT

`recipe.yaml` builds the conda package `alexandria` (programs, libraries and
data) with [rattler-build](https://rattler.build). The plan and the decisions
behind it are in `docs/conda-packaging-plan.md`.

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

    conda create -n alex-test -c file://$PWD/output-conda -c conda-forge alexandria
    conda activate alex-test
    alexandria -version

## Uploading

Test builds go to the channel `alexandriachemistry` on anaconda.org with the
label `dev`, releases with the label `main`:

    anaconda -t <token> upload -u alexandriachemistry -l dev output-conda/*/alexandria-*.conda

The workflow `.github/workflows/act-conda.yml` builds `act-utilities` and
`alexandria` on macOS (Apple silicon) and Linux (x86-64 and 64-bit ARM),
and uploads them
with the secret `ANACONDA_API_TOKEN`:

- a pull request that changes the recipe builds and tests the packages,
  without uploading them;
- a tag `v*` uploads them with the label `main`;
- a run by hand (Actions, Conda packages, Run workflow) uploads them if
  *upload* is ticked, with the label given, `dev` by default.

The packages of each run are also kept as artifacts of the run.
Releases on conda-forge are built by the feedstock `alexandria-feedstock`
from the archive of the release tag.
