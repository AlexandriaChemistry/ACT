***********************
Installation of the ACT
***********************
The Alexandria Chemistry Toolkit (ACT) relies on a number of libraries. Even though we tried to keep it to a minimum, some more or less standard libraries are needed. ACT should compile fine on any UNIX (including MacOs) or Linux machine. 
%Most of the libraries can be installed using Anaconda or even Miniconda which has the advantage of running in user-space entirely, that is you do not super-user access to install it. 

=============
Prerequisites
=============
The following software packages are required for the ACT to work:

* C and C++ compilers supporting C++20 at least. On Linux a GNU c++ version newer than 12.0 is recommended.
* Some version of a library that supports the message passing interface (MPI) for parallel programming. A popular version is OpenMPI.
* The cmake tools (at least version 3.13.0) are needed for compiling the code.
* For linear algebra operations we use the Eigen library, version 5 or better.
* The RDKit library (at least version 2025.09.4)
* The Boost developer library of version 1.86.0 (matching RDKit)
* The LibXml2 is needed for processing the XML data files used by the ACT.
* Python, version 3.8 or better, and a number of Python libraries, namely, NumPy, Matplotlib, and PubChemPy.


The following libraries are optional only, but may be useful for developing the ACT:

* The Class Library for Numbers is used in an optional part of the code (Slater-distributed charges :cite:p:`Ghahremanpour2018a`) and can be omitted.
* The SQLite database engine is needed to process experimental data as well as quantum chemistry data from the Alexandria Library :cite:p:`Ghahremanpour2018a`, available from `Zenodo`_.

For developers, the following additional packages are needed:

* doxygen, it is used for generating software documentation,
* graphviz can be optionally added for generating graphs and tree-structures in the doxygen documentation,
* pygments, for source code listings,
* sphinx, which is used for building the manual, and
* sphinxcontrib-bibtex, for the references.

.. _Zenodo: https://doi.org/10.5281/zenodo.1004710

=================
Conda Environment
=================
There are multiple ways to fulfill the prerequisites.
The simplest way that should suffice on a single computer (i.e. not a cluster), is to first install miniconda on your computer, download the ACT conda environment file `Yml`_ and create and activate a new conda environment.::

  conda create -n ACT
  conda activate ACT
  conda config --add channels anaconda
  conda config --add channels conda-forge
  conda install librdkit-dev=2025.09.4 libboost-devel=1.86.0 cmake eigen=5 libxml2 numpy matplotlib pubchempy pillow xmltodict

This should install the libraries mentioned above (note: it will take some time!). If you are installing ACT on a high-performance computing cluster, there likely is support for compilers and a MPI library already. If not, then add the *openmpi* package to your conda install line.
Most Linux installations come bundled with the GNU compiler suite (`GCC`_) and for macOS the Xcode package can be downloaded free of charge from `Xcode`_. If you do not have a compiler, add *gcc* to the conda install line::

  conda install gcc=14 gxx=14 openmpi=5

For developers, please additionally install these packages::

  conda install doxygen graphviz pygments sphinx sphinxcontrib-bibtex

.. _Yml: https://github.com/dspoel/ACT/blob/main/ACT.yml
.. _GCC: https://gcc.gnu.org
.. _Xcode: https://developer.apple.com/xcode/ 

.. attention:: If you are installing ACT in a cluster we recommend to use the cluster-provided compilers and in particular the MPI library since it may be tuned to make optimal use of the communication hardware. 
High-performance computer centers typically provide compilers libraries using some kind of module system.

========================
Running the Installation
========================
Once you have download ACT either as a release version or by cloning the git repository, enter the directory where the ACT source code is installed, and issue the following commands::

  mkdir build-Release
  cd build-Release
  cmake -DCMAKE_INSTALL_PREFIX=$HOME/tools ..
  make -j 8 install

where the *CMAKE_INSTALL_PREFIX* flag points to the location where ACT will be installed and the *-j 8* flag instructs the make command to utilize 8 cores to speed up compilation.


In order to start using the software, run the following command::

  source $HOME/tools/bin/ACTRC

or add it to your .bash_profile (or equivalent, for remote machines) or .bashrc (or equivalent, for local machines), and restart the shell or log in again.
Then you can run the alexandria executable using::

  alexandria -h

To make sure you do have the correct commands in your path, please try the command::

  which alexandria

which should give you something like::

  % which alexandria
  ~/tools/bin/alexandria

There are some building options available that are mainly of use for developers (see :ref:`tab-cmake`). These options have to be specified using::

  -DOPTION=Value
  
where *Value* can be *ON* or *OFF* or something more option specific. 

.. table:: cmake flags available to build the alexandria program.
   :name: tab-cmake

   +----------------------+------------------------------------------------------------+
   | **Flag**             | **Description**                                            |
   +======================+============================================================+
   | CMAKE_BUILD_TYPE     | Build Type: either Release (default) Debug (for use with   |
   +----------------------+------------------------------------------------------------+
   |                      | a debugger) or ASAN (Adress Sanitizer, for uncovering      |
   +----------------------+------------------------------------------------------------+
   |                      | memory leaks and debugging crashes).                       |
   +----------------------+------------------------------------------------------------+
   | CMAKE_INSTALL_PREFIX | Path where to install the ACT, see above                   |
   +----------------------+------------------------------------------------------------+
   | CMAKE_PREFIX_PATH    | Path where cmake can look for libraries. Multiple paths    |
   +----------------------+------------------------------------------------------------+
   |                      | can be specified, separated by semicolons, for instance    |
   +----------------------+------------------------------------------------------------+
   |                      | -DCMAKE_PREFIX_PATH=${CONDA_PREFIX}/lib                    |
   +----------------------+------------------------------------------------------------+
   | ACT_CLN              | Install ACT using the Class Library for Numbers for high   |
   +----------------------+------------------------------------------------------------+
   |                      | precision calculations, used for Slater integrals.         |
   +----------------------+------------------------------------------------------------+
   |                      | Default OFF, activated when ON.                            |   
   +----------------------+------------------------------------------------------------+
   | ACT_BUILD_MANUAL     | Whether or not to provision for building the manual.       |
   +----------------------+------------------------------------------------------------+
   |                      | Requires installing developer tools, default OFF.          |
   +----------------------+------------------------------------------------------------+
   
   
============================
Troubleshooting Installation
============================
Sometimes, the *cmake* process or building using *make* does not work as described above.
For instance, there may be library mismatches like this::

   undefined reference to `__cxa_call_terminate@CXXABI_1.3.15'
   
which is caused by the fact that conda libraries expect specific versions of system libraries, combined with the make process supplying a wrong version of that library.
In that case you can instruct cmake to change the order of libraries by adding these flags to the cmake command line::

  


===============
Testing the ACT
===============

To start testing, you first want to familiarize yourself with the test set. If the ACT is in your home directory, you can::

  cd ACT/build-Release

Then you can build the test set using::

  make tests

and run it using::

  make test

which should give the following output::

  Running tests...
  Test project /Users/spoel/GG/ACT/build_Release_DOUBLE
        Start  1: TestUtilsUnitTests
   1/17 Test  #1: TestUtilsUnitTests ...............   Passed    3.10 sec
        Start  2: WangBuckinghamTests
   2/17 Test  #2: WangBuckinghamTests ..............   Passed    0.33 sec
        Start 16: AlexandriaTests

(more tests)::

  16/17 Test #16: AlexandriaTests ..................   Passed    8.80 sec
        Start 17: SobolTests
  17/17 Test #17: SobolTests .......................   Passed    0.16 sec
  
  100% tests passed, 0 tests failed out of 17

You can also run an individual test, like this bin/sobol-test which should give this output::

  [==========] Running 2 tests from 1 test case.
  [----------] Global test environment set-up.
  [----------] 2 tests from SobolTest
  [ RUN      ] SobolTest.Test08
  [       OK ] SobolTest.Test08 (0 ms)
  [ RUN      ] SobolTest.Test09
  [       OK ] SobolTest.Test09 (0 ms)
  [----------] 2 tests from SobolTest (0 ms total)
  
  [----------] Global test environment tear-down
  [==========] 2 tests from 1 test case ran. (0 ms total)
  [  PASSED  ] 2 tests.

Note that these tests are run every time a change in the ACT source code is uploaded to github, to prevent errors in the code from being introduced.

