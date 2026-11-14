************************
Verifying changes to ACT
************************
The test set described in the previous chapter checks that ACT works as intended. When changing the code, two scripts in the *admin* directory of the source code help to check more thoroughly that nothing broke: *verify.sh* runs the tests in different ways and compares a set of real calculations with a reference, and *coverage.sh* reports which parts of the code the tests do not run. In addition, *make tidy* cleans up the code with *clang-tidy*, see the last section.

==========================
Running the test set fully
==========================
From the directory with the ACT source code, the command::

  admin/verify.sh

builds ACT in the build directory *build_Release_DOUBLE*, including the test programs, and runs all tests twice: once with the explicit SIMD routines for the energies and forces and once without, by setting *ACT_NO_EXPLICIT_SIMD* (see :numref:`tab-cmake`). Since both ways should give the same results within the tolerances of the tests, this checks the explicit SIMD routines against the others. The options of the script are listed in :numref:`tab-verify`.

The script does not configure the build directories, except the one without vectorization of *-\-novec*: it uses the *cmake* and *ctest* programs with which a directory was configured, and the MPI launcher found there. Configure them once, from the directory with the source code, with the compilers and libraries you normally use, for instance::

  cmake -S . -B build_Release_DOUBLE -DCMAKE_BUILD_TYPE=Release -DACT_DOUBLE=ON
  cmake -S . -B build_Release_SINGLE -DCMAKE_BUILD_TYPE=Release -DACT_DOUBLE=OFF
  cmake -S . -B build_ASAN_DOUBLE -DCMAKE_BUILD_TYPE=ASAN -DACT_DOUBLE=ON

where the second one is needed for *-\-single* and the third one for *-\-asan*. Another directory than *build_Release_DOUBLE* can be given with *-b*; the directories of *-\-single* and *-\-asan* always have these names.

.. table:: Options of admin/verify.sh.
   :name: tab-verify

   +------------------+------------------------------------------------------------+
   | **Option**       | **Description**                                            |
   +==================+============================================================+
   | -b DIR           | The build directory, by default build_Release_DOUBLE.      |
   +------------------+------------------------------------------------------------+
   | -j N             | The number of parallel build and test jobs, default 8.     |
   +------------------+------------------------------------------------------------+
   | -\-novec         | Also build in DIR_NOVEC without any vectorization, that is |
   |                  | with ACT_SIMD=NONE and with vectorization by the compiler  |
   |                  | turned off, and run the tests there with and without       |
   |                  | explicit SIMD.                                             |
   +------------------+------------------------------------------------------------+
   | -\-asan          | Also build and test build_ASAN_DOUBLE, with the address    |
   |                  | sanitizer, which finds errors in the use of memory.        |
   +------------------+------------------------------------------------------------+
   | -\-single        | Also build and test build_Release_SINGLE, in single        |
   |                  | precision, with and without explicit SIMD.                 |
   +------------------+------------------------------------------------------------+
   | -\-gpu           | Also build and test the GPU build, on a computer with a    |
   |                  | GPU, where its tests may not be skipped, and run the       |
   |                  | regression set there as well, see below.                   |
   +------------------+------------------------------------------------------------+
   | -g DIR           | The GPU build directory, by default                        |
   |                  | build_Release_SINGLE_SYCL.                                 |
   +------------------+------------------------------------------------------------+
   | -\-save-reference| Run the regression set, see below, and keep the results as |
   |                  | the reference.                                             |
   +------------------+------------------------------------------------------------+
   | -\-regression    | Run the regression set and compare it with the reference.  |
   +------------------+------------------------------------------------------------+
   | -\-tolerance T   | The largest relative difference from the reference that is |
   |                  | accepted, default 0, that is identical results.            |
   +------------------+------------------------------------------------------------+
   | -\-no-tests      | Only build and run the regression set.                     |
   +------------------+------------------------------------------------------------+

The script ends with a list of what failed, if anything, and has a nonzero exit status in that case. The output of the build and of the tests is in the files *verify-build.log* and *verify-tests-simd.log* or *verify-tests-nosimd.log* in each build directory, and that of configuring the directory of *-\-novec* in its *verify-configure.log*. A test that fails is listed with its output in these files.

Which checks to run depends on the change. Before committing a small change, *admin/verify.sh* is usually enough. Before merging a branch, or after changes to the energy and force routines, also check the other builds::

  admin/verify.sh --novec --asan --single

The build without vectorization shows whether a difference comes from vectorization, the address sanitizer finds reading or writing outside arrays and memory leaks, and the single precision build checks the tolerances of the tests in single precision, which is used with a GPU.

After changes to the GPU code, run the script on a node with a GPU with *-\-gpu* as well. The GPU build directory has to be configured with *-DACT_GPU=SYCL -DACT_DOUBLE=OFF*, see :numref:`tab-cmake`. The tests of the GPU, whose names end in *Gpu*, run as separate *ctest* entries, for instance *ACTForcesTestsGpu*, which *ctest* reports as skipped without a GPU. With *-\-gpu*, a skipped GPU entry counts as a failure, since then no GPU was found. With *-\-save-reference* or *-\-regression*, the regression set is also run in the GPU build, with the polarizable models trained on the GPU, and compared with a reference made with that build (mode *gpu*).

==================
The regression set
==================
The tests check small parts of the code. The regression set checks whole calculations: training of the hydrogen fluoride dimer, of the halgren-PG force field for alcohols and water, of the polarizable alc-pol force field and of bonded parameters for monomers, as well as the second virial coefficient of the BK3, TIP5P and hydrogen fluoride models. These are small versions of the benchmarks of ACT and take about a minute together. They need the input files of the ACTdata repository, by default in *$HOME/Alexandria/ACTdata*, which can be changed with the environment variable *ACT_VERIFY_DATA*: the force fields, molecules and selections in *Performance/train_ff*, the dimer structures in *Performance/b2*, and the force fields *BK3.xml* and *TIP5P.xml* in *LegacyForceFields*. Each calculation is run with 8 MPI processes, with the MPI launcher found when the build directory was configured, so the build needs MPI (*ACT_MPI*, on by default).

The regression set is useful for changes that should not change the results, such as making the code faster or easier to maintain. Before such a change, keep the results as the reference::

  admin/verify.sh --save-reference

and after the change compare with it::

  admin/verify.sh --regression

Both run the regression set twice, with and without explicit SIMD, after building and testing as above unless *-\-no-tests* is given. For each calculation, the script prints whether all parameters of the trained force field (in *ff_inter.xml*), or all second virial coefficients (in *b2.xvg*), are identical to those of the reference, or else the largest relative difference, for instance::

  === Regression set, simd
      bk3-b2           b2.xvg         identical
      halgren-train    ff_inter.xml   identical
      hf-b2            b2.xvg         identical
      hf-train         ff_inter.xml   identical
      ljpc-intra       ff_inter.xml   identical
      pol-train        ff_inter.xml   identical
      tip5p-b2         b2.xvg         identical

and the same for *nosimd*, followed by *=== All passed*. A difference is printed as, for instance, *relative difference 3.2e-07*.
A calculation fails if the difference is larger than the tolerance (*-\-tolerance*, default 0), if its output is missing, or if it has a different number of values than the reference. The results are in *verify/simd* and *verify/nosimd* in the build directory, the reference in *verify-reference*. If a calculation fails while saving the reference, the old reference is kept. The output of each calculation is in *out.txt* in its directory, for instance *verify/simd/hf-train/out.txt*. The training uses a genetic algorithm with random numbers, which amplifies the smallest differences in the energies, so the trained force fields can only be compared within the same build and with explicit SIMD either on or off. A change that only reorders floating point operations, for instance, can change the trained force field considerably, while each energy changes in the last digits only. For such a change, the tests show whether it is correct, and a tolerance can be given, which applies to the trained force fields and the second virial coefficients alike; in practice it is only useful for the latter.

=============
Test coverage
=============
The command::

  admin/coverage.sh

builds ACT in *build_COVERAGE* with the source-based coverage of the clang compiler, runs all tests and reports which lines and functions of ACT they run, leaving out the GROMACS and external code and the tests themselves. It writes three reports in *build_COVERAGE/coverage*:

* *summary.txt*, the coverage of each directory of *src/act* and the files with the most lines that are not run by any test,
* *report.txt*, the coverage of each file, and
* *html/index.html*, the source code with the number of times each line was run, in which code that is not run by any test stands out in red.

This shows where new tests are most needed. The script needs the clang compiler and its tools *llvm-profdata* and *llvm-cov*, which come with Xcode on macOS.

====================================
Cleaning up the code with clang-tidy
====================================
The build target *tidy* cleans up the code of *src/act* with the *clang-tidy* program of the LLVM project and changes the source files in place. From a configured build directory, the command::

  make tidy

runs *clang-tidy* on all source files of *src/act* that are compiled in that build directory, collects the fixes it proposes and applies them. The checks are listed in *admin/clang-tidy-fix.yaml*; only checks that come with automatic fixes are used. They

* make local variables and references *const* where they are not changed, but only for variables that are not scalars, such as objects, strings, containers and iterators; variables of the built-in types, enumerations and pointers are left as they are (see *-\-const-scalars* below), as are the parameters of functions, so that no function signatures change,
* avoid needless copies, for instance by using *const auto &* for a variable that is copied from a constant reference, or *emplace_back* rather than *push_back*,
* use *empty()* rather than comparing *size()* with zero, and remove calls of *c_str()* and initializations of strings that are not needed,
* write *= default* for empty constructors and destructors and *override* for virtual functions that override another,
* simplify boolean expressions, use *nullptr*, and add braces around the statement of *if*, *for*, *while* and *else* when they are missing.

The fixes are applied in the layout of the ACT code, with 4 spaces of indentation: *const* is written in front of the type (*const auto it*), aligned declarations stay aligned, braces that are added get lines of their own, and *= default* follows the parameter list. Only the lines that are fixed change, the rest of each file is left as it is. Fixes that overlap with another fix are left for the next run, and fixes that cannot be written in this layout, such as braces in a function written on one line, are left out; the output lists how many fixes of each check were applied or left. Run *make tidy* until it reports no more fixes.

The clean-up can be limited to some files or directories, relative to the directory with the source code, and the options of :numref:`tab-tidy` can be given, with two environment variables::

  make tidy ACT_TIDY_PATHS=src/act/forces
  make tidy ACT_TIDY_PATHS="src/act/topology src/act/qgen" ACT_TIDY_ARGS=--dry-run

Header files are changed only when they are in the selected files or directories. The script that does the work, *admin/clang_tidy_fix.py*, can also be run directly, with the build directory given as *-\-build-dir*.

.. table:: Options of make tidy, given with ACT_TIDY_ARGS.
   :name: tab-tidy

   +------------------+------------------------------------------------------------+
   | **Option**       | **Description**                                            |
   +==================+============================================================+
   | -\-dry-run       | Only report how many fixes of each check would be applied, |
   |                  | without changing any file.                                 |
   +------------------+------------------------------------------------------------+
   | -\-const-scalars | Also make variables of scalar types (built-in types,       |
   |                  | enumerations and pointers) const.                          |
   +------------------+------------------------------------------------------------+
   | -j N             | The number of clang-tidy processes run at the same time,   |
   |                  | by default the number of cores.                            |
   +------------------+------------------------------------------------------------+

*clang-tidy* reads how each file is compiled from *compile_commands.json* in the build directory, which *cmake* writes when it configures the build. It therefore only sees the code that is compiled there: code for another configuration, for instance the SYCL code of a build with GPU support or code for single precision, is only cleaned up from a build directory with that configuration. The target needs *clang-tidy* and Python with the *PyYAML* module. Any recent version of *clang-tidy* will do; if it is not in the path, or has another name, give it when configuring::

  cmake -DACT_CLANG_TIDY=$HOME/miniforge3/envs/tidy/bin/clang-tidy .

where, for instance, the conda environment *tidy* was made with::

  conda create -n tidy -c conda-forge clang-tools clang

The package *clang* is needed for the header files of the compiler that *clang-tidy* uses. As with any change, check the result with *git diff*, then build ACT and run the tests, for instance with *admin/verify.sh*, before committing it.

======================
Code testing on GitHub
======================
When a branch is pushed to github five `CI`_ (continuous integration) tests are deployed, so-called *github-workflows*. Four of them download ACT to a test machine, build it and run the tests, all of which should pass. If a single test does not pass the branch cannot be included into the main branch.

.. _CI: https://en.wikipedia.org/wiki/Continuous_integration
 
The four different tests are

+ Double precision build on Linux
+ Double precision build on MacOS
+ Single precision build on Linux with explicit SIMD
+ Double precision build on Linux without explicit SIMD
+ *make tidy* on macOS, which fails if it would change the code of the directories of *src/act* that have been cleaned up with it, listed in *.github/workflows/act-tidy.yml*. It uses *clang-tidy* 23.1.5 in single precision, as other versions or configurations may propose other fixes. If it fails, run *make tidy* on the directories it names and commit the result.

If a pull-request is made to the main branch, an additional tests will be done, to verify that the manual will build. This test too, has to pass.

Finally, if a branch is merged into main, all tests are run again, and the both the on-line manual and the pdf manual are updated.
