***********************
Installation of the ACT
***********************
The Alexandria Chemistry Toolkit (ACT) relies on a number of libraries. Even though we tried to keep it to a minimum, some more or less standard libraries are needed. ACT should compile fine on any UNIX (including MacOs) or Linux machine. 
%Most of the libraries can be installed using Anaconda or even Miniconda which has the advantage of running in user-space entirely, that is you do not super-user access to install it. 

=============
Prerequisites
=============
The following software packages are required for the ACT to work:

* C and C++ compilers supporting C++20 at least. On Linux a GNU c++ version newer than 13.0 is recommended.
* Some version of a library that supports the message passing interface (MPI) for parallel programming. A popular version is OpenMPI.
* The cmake tools (at least version 3.13.0) are needed for compiling the code.
* For linear algebra operations we use the Eigen library, version 5 or better.
* The RDKit library (at least version 2025.09.4)
* The LibXml2 is needed for processing the XML data files used by the ACT.
* Python, version 3.8 or better, and a number of Python libraries, namely, NumPy, Matplotlib, and PubChemPy.


The following libraries are optional only, but may be useful for developing the ACT:

* The Class Library for Numbers is used in an optional part of the code (Slater-distributed charges :cite:p:`Ghahremanpour2018a`) and can be omitted.
* The SQLite database engine is needed to process experimental data as well as quantum chemistry data from the Alexandria Library :cite:p:`Ghahremanpour2018a`, available from `Zenodo`_.
* The doxygen package can be used for generating documentation. In that case you need the Graphviz package as well.

.. _Zenodo: https://doi.org/10.5281/zenodo.1004710

=================
Conda Environment
=================
There are multiple ways to fulfill the prerequisites.
The simplest way that should suffice on a single computer (i.e. not a cluster), is to first install miniconda on your computer, download the ACT conda environment file `Yml`_ and create and activate a new conda environment.::

  conda env create -f ACT.yml
  conda activate ACT

This should install the libraries mentioned above (note: it will take some time!). You still need a C++ compiler but most Linux installations come bundled with the GNU compiler suite (`GCC`_) and for macOS the Xcode package can be downloaded free of charge from `Xcode`_.

.. _Yml: https://github.com/dspoel/ACT/blob/main/ACT.yml
.. _GCC: https://gcc.gnu.org
.. _Xcode: https://developer.apple.com/xcode/ 

If you are installing ACT in a cluster we recommend to use the cluster-provided compilers and in particular the MPI library since it may be tuned to make optimal use of the communication hardware. High-performance computer centers typically provide compilers libraries using some kind of module system.

========================
Running the Installation
========================
Once the prerequisites  are met, the easiest way to get going is to fetch the 
`install_act`_ script to your working directory of choice.
To fetch the script, once you have clicked on the link, click on the download icon (Fig. :ref:`fig-download`).

.. _install_act: https://github.com/dspoel/ACT/blob/main/src/act/python/install_act

.. figure:: ../images/link.png
   :align: center
   :name: fig-download
   
   Image showing where to click to download the installation script.

Once you have downloaded the script, start by executing::

  ./install_act -h

and studying the options.
Let's say you have a four core machine and you want to install first time around, then this command should do the trick::

  ./install_act -ncores 4

In the best of worlds, the script will have created a directory under your home catalog, with the name tools. 

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

===============
Testing the ACT
===============

To start testing, you first want to familiarize yourself with the test set. If the ACT is in your home directory, you can::

  cd ACT/build_Release_DOUBLE

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

