************************
Verifying changes to ACT
************************
The test set described in the previous chapter checks that ACT works as intended. When changing the code, two scripts in the *admin* directory of the source code help to check more thoroughly that nothing broke: *verify.sh* runs the tests in different ways and compares a set of real calculations with a reference, and *coverage.sh* reports which parts of the code the tests do not run.

==========================
Running the test set fully
==========================
From the directory with the ACT source code, the command::

  admin/verify.sh

builds ACT in the directory *build_Release_DOUBLE*, including the test programs, and runs all tests twice: once with the explicit SIMD routines for the energies and forces and once without, by setting *ACT_NO_EXPLICIT_SIMD* (see :numref:`tab-cmake`). Since both ways should give the same results within the tolerances of the tests, this checks the explicit SIMD routines against the others. The options of the script are listed in :numref:`tab-verify`.

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

The script ends with a list of what failed, if anything, and has a nonzero exit status in that case. The output of the build and of the tests is in the files *verify-build.log* and *verify-tests-simd.log* or *verify-tests-nosimd.log* in the build directory.

==================
The regression set
==================
The tests check small parts of the code. The regression set checks whole calculations: training of the hydrogen fluoride dimer, of the halgren-PG force field for alcohols and water, of the polarizable alc-pol force field and of bonded parameters for monomers, as well as the second virial coefficient of the BK3, TIP5P and hydrogen fluoride models. These are small versions of the benchmarks of ACT and take about a minute together. They need the input files of the ACTdata repository, by default in *$HOME/Alexandria/ACTdata*, which can be changed with the environment variable *ACT_VERIFY_DATA*.

The regression set is useful for changes that should not change the results, such as making the code faster or easier to maintain. Before such a change, keep the results as the reference::

  admin/verify.sh --save-reference

and after the change compare with it::

  admin/verify.sh --regression

For each calculation, the script prints whether all parameters of the trained force field, or all second virial coefficients, are identical to those of the reference, or else the largest relative difference. The results are in *build_Release_DOUBLE/verify*, the reference in *build_Release_DOUBLE/verify-reference*, both with and without explicit SIMD. The training uses a genetic algorithm with random numbers, which amplifies the smallest differences in the energies, so the trained force fields can only be compared within the same build and with explicit SIMD either on or off. A change that only reorders floating point operations, for instance, can change the trained force field considerably, while each energy changes in the last digits only. For such a change, the tests show whether it is correct, and a tolerance can be given for the second virial coefficients.

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
