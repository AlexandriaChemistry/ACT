# Plan: conda packages for ACT

Status: draft, 2026-11-17

Goal: users can install ACT with

    conda install -c conda-forge alexandria

or, for test builds and pre-releases,

    conda install -c alexandriachemistry -c conda-forge alexandria

and get working `alexandria`, `train_ff` and the other programs, the data
files in `share/act`, and the Python scripts. Nothing is compiled on the
user's machine.

Targets of the first release:

| Target | conda platform | Explicit SIMD | GPU | Runs on |
|---|---|---|---|---|
| A | `osx-arm64` | ARM NEON | none | Macs with Apple silicon (M1 and later) |
| B | `linux-64` | AVX2 with FMA | SYCL (AdaptiveCpp), CUDA backend, plus the CPU | x86-64 Linux with AVX2 (Intel Haswell, AMD Zen and later), with or without an NVIDIA GPU |
| C | `linux-aarch64` | ARM NEON | SYCL (AdaptiveCpp), CUDA backend, plus the CPU | 64-bit ARM Linux, such as NVIDIA Grace (GH200 nodes of arrhenius), with or without an NVIDIA GPU |

All are built from one recipe, as variants; `linux-aarch64` is built in the
workflow on a native ARM runner. SVE is not used, since ACT has no SVE
kernels. Other targets (section 10) only add variants.

## 1. Decisions

| Question | Proposal | Why |
|---|---|---|
| Build tool | `rattler-build` with a v1 `recipe.yaml` | It is fast, conda-forge accepts it, and it runs the same way on a Mac, on Linux and in GitHub Actions. |
| Package names | `alexandria` (programs, libraries, data) and `alexandria-python` (scripts) | `act` is taken on conda-forge. `alexandria` and `alexandria-python` were free on 2026-11-17. A package `alexandria` exists in the channel `astrorama` (a Euclid library), but that does not conflict with conda-forge. |
| Channels | Releases on conda-forge, through a pull request to `conda-forge/staged-recipes`. Test builds, pre-releases and urgent fixes on our own channel `alexandriachemistry` on anaconda.org. | Neither channel has license terms (only Anaconda's own `defaults` channel does). conda-forge gives users one channel and rebuilds automatically when RDKit, Boost or AdaptiveCpp change ABI. Our channel can be updated the same day, without review. The name `alexandria` on anaconda.org belongs to someone else. `alexandriachemistry` was free. |
| Maintainers | David van der Spoel (`dspoel`, main developer), Julian Köngeter (`JulianKoengeter`) and Måns Rosenbaum (`MansRosenbaum`) | conda-forge lists the maintainers in the recipe; they review the bot's pull requests for new versions and ABI rebuilds. |
| Dependencies | conda-forge only, no Homebrew or system packages | The package must only link to libraries that conda installs. The current CI mixes `brew`/`apt` (compilers, MPI, libxml2) with conda, so it does not test this. |
| Compiler | conda-forge `{{ compiler('cxx') }}`: clang with libc++ on macOS, GCC 13 or newer on Linux | This is needed for binary compatibility with `librdkit` and Boost. On Linux, AdaptiveCpp needs the C++ library of GCC 13 or newer. |
| Precision | Single precision (the default) in `alexandria`. Double precision optional, with the `_d` suffix. | The GPU only works in single precision. The suffixes let both be installed in one environment. |
| MPI | Variants `mpich` and `openmpi`. On Linux, document the `external` builds of conda-forge, which let users use the MPI library of their cluster. | `train_ff` needs MPI. On clusters, the MPI of conda does not use the fast network. |
| OpenMP | Off (`ACT_OPENMP=OFF`, the default) | It is not used for performance. |
| Python scripts | Separate output `alexandria-python` (`noarch: python`), required by `alexandria` | The scripts are currently installed into `lib/pythonX.Y/site-packages/act`, which would tie the binary package to one Python version. |

## 2. Target A: macOS on Apple silicon with NEON

Facts from `cmake/actManageSimd.cmake`:

- NEON is part of every aarch64 processor, so `ACT_SIMD=ARM_NEON` needs no
  compiler flags. The package runs on all Apple silicon processors.
- `ACT_SIMD=AUTO` runs a `try_run` program. When cross-compiling, this is
  impossible and AUTO falls back silently to **NONE**. conda-forge builds
  `osx-arm64` packages by cross-compiling.

Settings: `-DACT_SIMD=ARM_NEON -DACT_GPU=NONE`. AdaptiveCpp has no Metal
backend, so there is no GPU on Macs.

Clang on arm64 contracts `a*b+c` into FMA by default. The unit tests already
pass this way on Apple silicon, both locally and in the macOS CI, which runs
on arm64 runners.

`${CMAKE_ARGS}` sets `MACOSX_DEPLOYMENT_TARGET=11.0`. Check that ACT needs no
libc++ feature of a newer macOS. If it does, raise `c_stdlib_version` in
`variants.yaml`.

## 3. Target B: Linux x86-64 with AVX2 and SYCL

### 3.1 SIMD

- Set `-DACT_SIMD=AVX2` explicitly. With AUTO, the build machine decides, and
  the GitHub and Azure runners have AVX-512.
- Also build the whole code with `-march=x86-64-v3` (AVX2, FMA, BMI). This
  lets the compiler vectorize all code with AVX2, not only
  `forces_explicit_simd.cpp`. conda-forge does this with the
  `x86_64-microarch-level` package (level 3). That package also refuses to
  install on processors without AVX2, which matches the name of the target.
  Measure the gain with the BFdb benchmark. If there is no gain, drop the
  microarch level and rely on the run-time check of the explicit SIMD.
- Measure AVX-512 against AVX2 on a Zen 4 and on an Intel processor. If
  AVX-512 pays off, add a level 4 variant later (section 10).

### 3.2 SYCL with AdaptiveCpp

- AdaptiveCpp is on conda-forge, with variants: `cpu` (OpenMP backend only)
  and `cuda` (also the CUDA backend). `cmake/actManageGpu.cmake` already
  lists the backends and warns if CUDA is missing.
- ACT uses the generic target of AdaptiveCpp. The kernels are compiled for
  the GPU at run time, so **no CUDA compiler is needed at build time**,
  only `acpp` and its LLVM. One binary then runs on any NVIDIA GPU, and
  without a GPU on the processor. This is why one package is enough.
- Host dependency `adaptivecpp=*=cuda*`. Its `run_exports` pin the run-time
  libraries and the `cuda-version`. The NVIDIA driver (`libcuda.so`) comes
  from the system. Add the virtual package `__cuda` only as a constraint
  (`run_constrained`), so the package also installs on machines without a
  GPU.
- Variant `gpu: [none, sycl_cuda]`. The users choose with
  `conda install alexandria=*=*sycl_cuda*`. Make `sycl_cuda` the default on
  `linux-64` (highest priority in the build string), since it also runs
  without a GPU.
- Does conda-forge have a HIP (ROCm) variant of AdaptiveCpp? If yes, add
  `gpu: sycl_hip` the same way, for AMD GPUs such as the MI250X on dardel.
  If not, AMD GPU users build from source, as now (manual, section "GPU support with SYCL").
- At run time, the programs already set `ACPP_DEBUG_LEVEL=1`, bind each
  process to one GPU (`CUDA_VISIBLE_DEVICES`), and say in the log file why
  they found no GPU. The package needs no further code for this.
- The manual (section "GPU support with SYCL") says that the multi-process service of NVIDIA is
  needed when several ranks share a GPU. This also goes into the conda
  section of the manual.

### 3.3 glibc and the system

- conda-forge builds against an old glibc sysroot (2.17 or 2.28). AdaptiveCpp
  on conda-forge decides which one. Use the same with `c_stdlib_version` in
  `variants.yaml`.
- The CUDA backend needs an NVIDIA driver that supports the `cuda-version`
  that AdaptiveCpp was built with. State the minimum driver in the manual.

### 3.4 MPI on clusters

The `mpich` package of conda-forge has `external` builds
(`mpich=*=external_*`), which contain no MPI library but use the one of the
system. This is the usual way to run conda packages with the MPI of a
cluster (Cray MPICH on dardel, MPICH on arrhenius, Open MPI on csb). Test it on
csb and dardel. The cluster workflow in the manual (`srun --mpi=pmi2`,
`--gpu-bind=single:N`) stays the same.

## 4. Changes in the ACT source tree (branch `conda-package`)

1. **`ACT_SIMD` without `try_run`.** Done: on `arm64`/`aarch64`, AUTO
   chooses `ARM_NEON` without running anything. On x86 when cross-compiling,
   AUTO warns and names `ACT_SIMD` instead of falling back silently to NONE.
2. **SIMD and GPU in the version output.** Done: `alexandria -version`
   printed only the welcome text, since the program runs its command manager
   quietly. It now prints the build information, including the new line
   `Explicit SIMD:` and `GPU support:`, so the package tests can check
   them. Note that all programs are subcommands of `alexandria`
   (`alexandria train_ff`).
3. **Relocatable install.** The programs find `share/act` relative to the
   executable through `getProgramContext().installationPrefix()` in
   `src/gromacs/utility/datafilefinder.cpp`. Check that this still works after
   conda moves the prefix. Check that the libraries use a relative rpath
   (`@loader_path/../lib` on macOS, `$ORIGIN/../lib` on Linux), including
   `libactGpu`, which acpp links. Check that `otool -L` or `ldd` show no
   build paths. Checked on this Mac: after moving the install, `alexandria`
   finds `share/act` and its own libraries through `@loader_path/../lib`.
4. **Files configured with the prefix.** `ACTRC`, `mol_csv_api.py`,
   `gauss2molprop`, `generate_mp`, `molselect` and `reshuffle_selection` come
   from `.cmakein` files. conda replaces the prefix in text files. Check that
   each one works in a different prefix. `ACTRC` is not needed in a conda
   environment.
5. **Python install directory.** Done: `-DACT_INSTALL_PYTHON=OFF` leaves
   out the modules and scripts. `-DACT_PYTHON_ONLY=ON` installs only them,
   without looking for any library, so the `alexandria-python` output needs
   only `cmake` and `python` at build time. A file `act.pth` in
   site-packages puts the modules (`import actutils` and so on) on the path
   of the Python of the environment, also for another Python version than
   the one at build time, so `ACTRC` is not needed in conda.
6. **Version.** The version is set by hand in `cmake/actVersionInfo.cmake`
   (now 1.9b). A GitHub tarball of a tag contains `admin/.isreposource`, so
   CMake adds `-dev`. The recipe passes `-DACT_BUILD_TARBALL=ON`, which
   removes it. No code change is needed.
6a. **Warnings are errors.** The top-level `CMakeLists.txt` sets
   `CMAKE_COMPILE_WARNING_AS_ERROR`. With the MPI wrapper of conda, the
   headers of the environment are included with `-I`, not as system headers,
   and the deprecation warnings in `nlohmann/json.hpp` stop the build. The
   recipe configures with `cmake --compile-no-warning-as-error`, so that new
   warnings of new compilers or libraries do not break the package.
7. **Optional dependencies.** SQLite3 is found automatically if present.
   Either depend on `libsqlite` or turn it off, so the build does not depend
   on what happens to be installed. CLN stays off.
8. **Tests without a GPU.** The GPU tests must skip, not fail, when no GPU is
   found (GitHub runners have none). With the `cuda` variant of AdaptiveCpp
   on a machine without a driver, the programs must fall back to the CPU
   quietly.

## 5. Recipe layout

New directory `admin/conda/` (later a feedstock repository):

    admin/conda/
      recipe.yaml      # outputs alexandria and alexandria-python
      variants.yaml    # per platform: mpi, gpu, precision, c_stdlib_version
      build.sh         # cmake configure, build, test, install
      README.md        # how to build locally and how to upload

`variants.yaml`:

```yaml
mpi: [mpich, openmpi]
precision: [single]          # add double later
gpu:
  - if: osx
    then: [none]
  - if: linux and x86_64
    then: [sycl_cuda, none]
```

Sketch of `recipe.yaml`:

```yaml
context:
  version: "1.9b"

recipe:
  name: alexandria-split
  version: ${{ version }}

source:
  url: https://github.com/AlexandriaChemistry/ACT/archive/refs/tags/v${{ version }}.tar.gz
  sha256: <fill in>

build:
  number: 0
  skip:
    - not ((osx and arm64) or (linux and x86_64))

outputs:
  - package:
      name: alexandria
    build:
      string: mpi_${{ mpi }}_${{ gpu }}_h${{ hash }}_${{ build_number }}
      script: build.sh
    requirements:
      build:
        - ${{ compiler('cxx') }}
        - ${{ compiler('c') }}
        - ${{ stdlib('c') }}
        - cmake >=3.24
        - ninja
        - if: build_platform != target_platform
          then: [ "${{ mpi }}" ]
      host:
        - libxml2-devel           # headers and libxml2.dylib; otherwise CMake takes the one of the macOS SDK
        - eigen >=5
        - nlohmann_json
        - libboost-headers
        - librdkit-dev
        - libsqlite
        - ${{ mpi }}
        - if: gpu == "sycl_cuda"
          then: [ "adaptivecpp * cuda*" ]
        - if: linux and x86_64
          then: [ "x86_64-microarch-level 3" ]
      run:
        - ${{ mpi }}
        - alexandria-python ==${{ version }}
      run_constrained:
        - if: gpu == "sycl_cuda"
          then: [ "__cuda" ]
    tests:
      - script:
          - alexandria -h
          - alexandria train_ff -h
          - test -f $PREFIX/share/act/alexandria.json
          - if: osx
            then: alexandria -version | grep -q "Explicit SIMD: *ARM_NEON"
          - if: linux
            then: alexandria -version | grep -q "Explicit SIMD: *AVX2"
          - if: gpu == "sycl_cuda"
            then: alexandria -version | grep -q "GPU support: *SYCL"
          - mpirun -n 2 alexandria train_ff <short run of the MPI test set>

  - package:
      name: alexandria-python
    build:
      noarch: python
      script: |
        cmake -S . -B build-python -DACT_PYTHON_ONLY=ON -DCMAKE_INSTALL_PREFIX=$PREFIX
        cmake --install build-python
    requirements:
      build: [ cmake >=3.24 ]
      host: [ python >=3.10 ]
      run: [ python >=3.10, numpy, matplotlib-base, rdkit, pyyaml ]
    tests:
      - python:
          imports: [ actutils, molprops ]

about:
  homepage: https://github.com/AlexandriaChemistry/ACT
  license: GPL-2.0-or-later AND LGPL-2.1-or-later   # ACT files GPL, GROMACS parts LGPL; COPYING has both texts
  license_file: COPYING
  summary: Alexandria Chemistry Toolkit, machine learning of force fields

extra:
  recipe-maintainers:
    - dspoel
    - JulianKoengeter
    - MansRosenbaum
```

OpenMM is not a dependency of `alexandria-python`. Users who want `act_openmm.py`
install `openmm` themselves.

## 6. Build script (`build.sh`)

```bash
set -euxo pipefail
mkdir build && cd build

case "${target_platform}" in
  osx-arm64) SIMD=ARM_NEON ;;
  linux-64)  SIMD=AVX2 ;;
esac
if [[ "${gpu}" == sycl_* ]]; then
  GPU="-DACT_GPU=SYCL -DAdaptiveCpp_DIR=${PREFIX}/lib/cmake/AdaptiveCpp"
else
  GPU="-DACT_GPU=NONE"
fi

cmake ${CMAKE_ARGS} --compile-no-warning-as-error -G Ninja .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
  -DBUILD_SHARED_LIBS=ON \
  -DACT_MPI=ON \
  -DACT_SIMD=${SIMD} \
  ${GPU} \
  -DACT_OPENMP=OFF \
  -DACT_DOUBLE=OFF \
  -DACT_INSTALL_PYTHON=OFF \
  -DACT_BUILD_TARBALL=ON

ninja
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-0}" != "1" ]]; then
  ninja tests                     # plain ninja leaves stale test binaries
  ctest --output-on-failure -j${CPU_COUNT}
fi
ninja install
```

## 7. Tests

1. **In the build** (native builds only): the 29 unit tests with `ctest`, as
   on every branch. The GPU tests skip without a GPU (section 4, item 8).
2. **In the package test** (fresh environment): the commands in
   `recipe.yaml`. These check the rpath, the data files, the SIMD, the GPU
   support, and a short `train_ff` run on two processes, compared with a
   reference chi2 with a tolerance.
3. **By hand, before the first upload:**
   - Target A: a Mac with M1 and one with M3 or M4, without Homebrew in
     `PATH`, running the tutorial of the manual.
   - Target B, CPU: a Linux machine with AVX2 and no GPU (a GitHub runner or
     a CPU node of arrhenius or dardel).
   - Target B, GPU: csb (x86-64, RTX A4500), with `mpich=*=external_*` and the
     MPI of the cluster. Compare the PC+GS benchmark (lockstep 4 and 8) with
     the build from source: the same chi2, and no slower.
   - Target B, AMD: dardel (MI250X), only if there is a `sycl_hip` variant.

## 8. Continuous integration and release

New workflow `.github/workflows/act-conda.yml`, triggered by a tag `v*` and
by hand (`workflow_dispatch`):

- Job `osx-arm64` on `macos-latest` (arm64). This is a native build that runs
  all tests.
- Job `linux-64` on `ubuntu-latest`, using the conda-forge sysroot. It builds
  both `gpu` variants. The GPU code is compiled, but only runs on the CPU
  backend there.
- Both jobs: install `rattler-build`, run
  `rattler-build build --recipe admin/conda/recipe.yaml -m admin/conda/variants.yaml`,
  and keep the packages as workflow artifacts.
- Upload to `alexandriachemistry` on anaconda.org (secret
  `ANACONDA_API_TOKEN`): on a tag `v*` with the label `main`, and by hand
  with the label `dev` for test builds
  (`conda install -c alexandriachemistry/label/dev ...`).
- Releases on conda-forge are built and uploaded by the feedstock, from the
  same tag. The recipe in `admin/conda/` and the one in the feedstock are
  kept the same; the feedstock version is the reference.
- Build time: about 10 to 20 minutes per variant, longer with acpp.
  Start with `mpich` only: two variants on Linux and one on macOS.
- The GPU check (section 7, item 3) stays a manual step on csb before each
  release, since GitHub has no GPU runners for us.

## 9. Documentation

- In `docs/manual/userman/installation.rst`, add a section "Installing with
  conda" before building from source. It covers the commands, choosing the
  MPI and GPU variants, `mpich=*=external_*` on clusters, the minimum NVIDIA
  driver, the multi-process service, and the `_d` programs if a double
  package is made.
- In `README.md`, add a short install line.
- Use minus signs, not em or en dashes, in the manual.

## 10. Later targets

- `linux-64` with AVX-512 (`x86_64-microarch-level 4`), if the measurements
  of section 3.1 show a gain.
- `linux-64` with `sycl_hip` for AMD GPUs, when AdaptiveCpp with HIP is on
  conda-forge.
- `osx-64` (Intel Macs, AVX2), only if users ask.
- Double precision packages for all targets, without GPU.

## 11. Order of work and estimate

| Step | Work | Estimate |
|---|---|---|
| 1 | Branch `conda-package`: changes of section 4, tests | 1–2 days |
| 2 | Relocation check on macOS and Linux: install, move the prefix, run the tests | 0.5 day |
| 3 | Recipe for target A, built locally on this Mac with `rattler-build` | 1 day |
| 4 | Target B without GPU, built in a Linux container or on csb | 1 day |
| 5 | Target B with `sycl_cuda`: linking with AdaptiveCpp from conda, rpath of `libactGpu`, test on csb with external MPI | 1–2 days |
| 6 | `alexandria-python` with `ACT_PYTHON_ONLY` (done in step 1) | - |
| 7 | Create the organization `alexandriachemistry` on anaconda.org; GitHub workflow that builds the recipe and uploads to it (label `dev`) | 1 day |
| 8 | Manual section, README | 0.5 day |
| 9 | Release tag, PR to `conda-forge/staged-recipes`, review by conda-forge (days to a few weeks) | 1 day plus waiting |
| 10 | Feedstock `alexandria-feedstock` created; later releases are PRs there (the bot opens them for new tags) | 0.5 day per release |
| 11 | (Later) further targets | open |

## 12. Open questions

1. Decided: package `alexandria` on conda-forge, plus the channel
   `alexandriachemistry` on anaconda.org. Maintainers: David van der Spoel
   (`dspoel`), Julian Köngeter (`JulianKoengeter`) and Måns Rosenbaum
   (`MansRosenbaum`), listed in the recipe under `extra`.
2. Decided: no double precision package in the first release. Later, ship
   both, the double one with the `_d` programs.
3. Decided: the `linux-64` package always has CUDA support (`sycl_cuda`);
   there is no variant without it. It also runs without a GPU.
4. Decided: the whole `linux-64` build is for `x86-64-v3` (AVX2), with
   `x86_64-microarch-level ==3`. AVX-512 may be added later if it pays off.
5. Decided: the Python modules and scripts move to the repository
   [ACT-utilities](https://github.com/AlexandriaChemistry/ACT-utilities),
   which gets its own conda package. ACT requires it: the package
   `alexandria` depends on it. This replaces the output `alexandria-python`
   and the options `ACT_INSTALL_PYTHON` and `ACT_PYTHON_ONLY`. See section 13.
6. Decided: the package version follows the ACT version, now 1.9b.
7. Done, the choices can still be changed: the Python package `act` of ACT
   (the alexandria program from Python and the classes of ACT, see
   `docs/python-api-gui-plan.md`) is the output `alexandria-python` of the
   recipe. The name was reserved for the scripts that went to
   `act-utilities`; the import name stays `act`. See section 14.

## 13. Moving the Python code to ACT-utilities

Decided and done: the package is `act-utilities`, the history of
`src/act/python` came along (`git filter-repo`), and CMake stops if the
modules cannot be imported. ACT-utilities must be merged before ACT, since
the CI of ACT installs it from its main branch. The steps:

1. Package name for ACT-utilities, for instance `act-utilities` (free on
   conda-forge), built `noarch: python` from a `pyproject.toml` there, with
   the modules as a package and the scripts as entry points or plain
   scripts.
2. Move `src/act/python` (modules, scripts and the tests
   `test_act_openmm.py` and `test_coords2molprop.py`) to ACT-utilities, with
   its history if wanted (`git filter-repo`).
3. In ACT: remove `src/act/python`, the two Python tests from CTest (or run
   them in ACT-utilities' CI with an installed `alexandria`), the Python part
   of `ACTRC`, and the options `ACT_INSTALL_PYTHON` and `ACT_PYTHON_ONLY`.
4. Source builds of ACT: CMake checks that ACT-utilities can be imported and
   stops with a message how to install it (`pip install` or `conda install`).
5. The manual: the chapters that use the scripts say where they come from.
6. The recipe: `alexandria` depends on the ACT-utilities package instead of
   `alexandria-python`.

## 14. The package `alexandria-python`

The Python package `act` has a compiled module, `act._core`, with the
classes of ACT (pybind11, CMake option `ACT_PYTHON`). A compiled module is
for one Python version, so it cannot be in `alexandria`, which does not
depend on Python, and it cannot be `noarch`.

| Topic | Choice | Reason |
|-------|--------|--------|
| Output | `alexandria-python` in the same recipe, built by `build-python.sh` | One version and one source for the program and its Python package. |
| Contents | The whole package `act`: the modules that run the program and the module `_core` | One package to install; `alexandria` stays free of Python files. |
| Python versions | The key `act_python` in `variants.yaml`, now 3.11 to 3.15 | The key `python` of conda-forge would also multiply `alexandria`, which has a Python in its build for a check of CMake only. The feedstock of conda-forge has to map its `python` key to this. |
| Libraries of ACT | Static, inside the module, without GPU support | Importing the module loads nothing of AdaptiveCpp or CUDA, and the module does not depend on how the libraries of `alexandria` were built. The price: ACT is compiled once for each Python version, and the classes do not use a GPU. Linking the module to the shared libraries of `alexandria` is the alternative when the classes should use the GPU. |
| MPI | The same variants as `alexandria` | The libraries of ACT call MPI, so the module links to the MPI library of the environment. |
| Dependency | The exact build of `alexandria` of the same recipe | The package runs its program and reads its data files in `share/act`. |
| Tests | `python/tests/test_core.py` in the build script, with the installed module and the program of `alexandria`; imports and `check-python-package.py` as package tests | The package test has only the two test files it needs, such that the package stays small. |
| Pull requests | One Python version (3.13); tags and runs by hand build all | Compile time. |

Not tested on a machine of ours: Linux, which only the workflow builds.
