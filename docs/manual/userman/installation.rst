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
* The nlohmann/json library is required for JSON I/O and can be installed via conda as ``nlohmann_json``.
* Python, version 3.8 or better, and a number of Python libraries, namely, NumPy, Matplotlib, and PubChemPy.


The following libraries are optional only, but may be useful for developing the ACT:

* The Class Library for Numbers is used in an optional part of the code (Slater-distributed charges :cite:p:`Ghahremanpour2018a`) and can be omitted.
* The SQLite database engine is needed to process experimental data as well as quantum chemistry data from the Alexandria Library :cite:p:`Ghahremanpour2018a`, available from `Zenodo`_.

For developers, the following additional packages are needed:

* doxygen, it is used for generating software documentation,
* graphviz can be optionally added for generating graphs and tree-structures in the doxygen documentation,
* pygments, for source code listings,
* sphinx, which is used for building the manual, with
* sphinxcontrib-bibtex, for the references and
* nbsphinx for including `Jupyter`_ Notebooks.

.. _Zenodo: https://doi.org/10.5281/zenodo.1004710
.. _Jupyter: https://jupyter.org/

=================
Conda Environment
=================
There are multiple ways to fulfill the prerequisites.
The simplest way that should suffice on a single computer (i.e. not a cluster), is to first install miniconda on your computer, download the ACT conda environment file `Yml`_ and create and activate a new conda environment.::

  conda create -n ACT
  conda activate ACT
  conda config --add channels anaconda
  conda config --add channels conda-forge
  conda install librdkit-dev>=2025.09.4 libboost-devel>=1.86.0 cmake eigen>=5 libxml2 numpy matplotlib pubchempy pillow nlohmann_json

This should install the libraries mentioned above (note: it will take some time!). If you are installing ACT on a high-performance computing cluster, there likely is support for compilers and a MPI library already. If not, then add the *openmpi* package to your conda install line.
Most Linux installations come bundled with the GNU compiler suite (`GCC`_) and for macOS the Xcode package can be downloaded free of charge from `Xcode`_. If you do not have a compiler, add *gcc* to the conda install line::

  conda install gcc=14 gxx=14 openmpi=5

For developers, please additionally install these packages::

  conda install doxygen graphviz pygments sphinx sphinxcontrib-bibtex nbsphinx

.. _Yml: https://github.com/dspoel/ACT/blob/main/ACT.yml
.. _GCC: https://gcc.gnu.org
.. _Xcode: https://developer.apple.com/xcode/ 

.. attention:: If you are installing ACT in a cluster we recommend to use the cluster-provided compilers and in particular the MPI library since it may be tuned to make optimal use of the communication hardware. High-performance computer centers typically provide compilers libraries using some kind of module system.

Building on a cluster
~~~~~~~~~~~~~~~~~~~~~
On a cluster, use the MPI library of the cluster, which its administrators built for its network and queue system, and conda only for the libraries that the cluster does not provide, such as RDKit. Leave *openmpi* out of the conda environment, or tell *cmake* which MPI to use with the compiler wrappers of the module, for instance::

  module load <compiler and MPI module of the cluster>
  conda activate ACT
  cmake .. -DMPI_C_COMPILER=$(which mpicc) -DMPI_CXX_COMPILER=$(which mpicxx) \
           -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} -DCMAKE_INSTALL_PREFIX=$HOME/tools

The C and C++ compilers may be those of the cluster or those of conda. The compilers of conda link against a system root of their own, in which the linker also looks for the libraries that other libraries depend on. The MPI library of a cluster depends on libraries outside it, for instance those of the network, and of CUDA if it is GPU aware. For an MPI library outside the conda environment, *cmake* therefore gives the linker the library directories of the system, of the MPI library, and of *LIBRARY_PATH* and *LD_LIBRARY_PATH*, which modules set. The output of *cmake* tells when it does so.

Start the programs with the launcher of the queue system, which places the processes on the nodes and cores of the job, for instance with Slurm::

  srun --mpi=pmi2 alexandria train_ff ...

Which *--mpi* option works depends on the MPI library, *srun --mpi=list* lists those that the cluster supports. The tests that use MPI are started with the *mpiexec* of the MPI library, unless *-DMPIEXEC_EXECUTABLE=srun* is given to *cmake*.

For example, on the GH200 nodes of arrhenius (NAISS), with the GPU toolchain of the cluster (GCC 14.3, MPICH 4.3.2, CUDA 13.0) and the libraries and AdaptiveCpp from conda, in a job on such a node since the login node has another processor::

  module load GPU/buildenv-gcccuda/2026.03-cu13.0
  export PATH=<conda environment>/bin:$PATH
  cmake .. -DACT_DOUBLE=OFF -DACT_GPU=SYCL -DCMAKE_PREFIX_PATH=<conda environment> \
           -DAdaptiveCpp_DIR=<AdaptiveCpp environment>/lib/cmake/AdaptiveCpp \
           -DMPI_C_COMPILER=$(which mpicc) -DMPI_CXX_COMPILER=$(which mpicxx)

and the programs are started with *srun --mpi=pmi2*.

On a Cray system, such as dardel (PDC), the compiler wrappers *cc* and *CC* of the Cray programming environment include the MPI library of the cluster already. Give them as the compilers and *cmake* finds MPI through them, nothing else is needed. For example, on the AMD MI250X nodes of dardel, with the Cray compilers, the libxml2 of the cluster and the other libraries and AdaptiveCpp from conda::

  cmake .. -DACT_DOUBLE=OFF -DACT_GPU=SYCL -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \
           -DCMAKE_PREFIX_PATH=<conda environment> \
           -DAdaptiveCpp_DIR=<conda environment>/lib/cmake/AdaptiveCpp \
           -DLIBXML2_LIBRARY=/pdc/software/eb/software/libxml2/2.12.3/lib/libxml2.so \
           -DLIBXML2_INCLUDE_DIR=/pdc/software/eb/software/libxml2/2.12.3/include/libxml2

The programs are started with plain *srun*. This libxml2 needs the ICU libraries of the cluster, so load them with *module load icu/74.1* in the job script as well.

========================
Running the Installation
========================
Once you have download ACT either as a release version or by cloning the git repository, enter the directory where the ACT source code is installed, and issue the following commands::

  mkdir build-Release
  cd build-Release
  cmake -DCMAKE_INSTALL_PREFIX=$HOME/tools ..
  make -j 8 install

where the *CMAKE_INSTALL_PREFIX* flag points to the location where ACT will be installed and the *-j 8* flag instructs the make command to utilize 8 cores to speed up compilation.


The programs and python scripts find the data files of ACT, in *share/act* of the installation, by themselves, also when started from the build directory, in which case they use those in *share* of the source code. To use other data files, set the environment variable *ACTDATA* to their directory; files that are not there are still taken from the installation. All that remains is to add the directory of the programs to your *PATH*, or to run the following command, which does that as well::

  source $HOME/tools/bin/ACTRC

Put either in your .bash_profile (or equivalent, for remote machines) or .bashrc (or equivalent, for local machines), and restart the shell or log in again.
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
   | ACT_DOUBLE           | Double precision, default ON. Single precision (OFF) is    |
   +----------------------+------------------------------------------------------------+
   |                      | faster and is needed to use a GPU. Before, this option was |
   +----------------------+------------------------------------------------------------+
   |                      | called GMX_DOUBLE; the old names of this and the next two  |
   +----------------------+------------------------------------------------------------+
   |                      | options still work, with a warning.                        |
   +----------------------+------------------------------------------------------------+
   | ACT_MPI              | Build with MPI, needed to train with several processes,    |
   +----------------------+------------------------------------------------------------+
   |                      | default ON.                                                |
   +----------------------+------------------------------------------------------------+
   | ACT_OPENMP           | OpenMP multithreading, default OFF. ACT parallelizes with  |
   +----------------------+------------------------------------------------------------+
   |                      | MPI processes, see ACT_OPENMP_SIMD for the vectorization   |
   +----------------------+------------------------------------------------------------+
   |                      | of loops.                                                  |
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
   | CMAKE_C_FLAGS        | By default nothing is needed, but sometimes this linker    |
   +----------------------+------------------------------------------------------------+
   |                      | flag helps to overcome problems:                           |
   +----------------------+------------------------------------------------------------+
   |                      | "-L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib"     |
   +----------------------+------------------------------------------------------------+
   | CMAKE_CXX_FLAGS      | See previous entry, you likely need to add both.           |
   +----------------------+------------------------------------------------------------+
   | ACT_SIMD             | Explicit SIMD instructions for the energy and force        |
   +----------------------+------------------------------------------------------------+
   |                      | routines: AUTO (default), NONE, ARM_NEON, AVX2 (with FMA,  |
   +----------------------+------------------------------------------------------------+
   |                      | for instance Intel Haswell and AMD Zen and later) or       |
   +----------------------+------------------------------------------------------------+
   |                      | AVX_512. AUTO selects the best SIMD for the processor that |
   +----------------------+------------------------------------------------------------+
   |                      | cmake runs on, or AVX2 if the compiler cannot produce      |
   +----------------------+------------------------------------------------------------+
   |                      | AVX-512 code. NONE turns explicit SIMD off. Only the       |
   +----------------------+------------------------------------------------------------+
   |                      | explicit SIMD routines are compiled for this SIMD, and     |
   +----------------------+------------------------------------------------------------+
   |                      | they are only used if the processor running the program    |
   +----------------------+------------------------------------------------------------+
   |                      | supports it. Otherwise the routines vectorized by the      |
   +----------------------+------------------------------------------------------------+
   |                      | compiler are used, such that the program also runs on      |
   +----------------------+------------------------------------------------------------+
   |                      | processors with less SIMD.                                 |
   +----------------------+------------------------------------------------------------+
   | ACT_PAIR_UNROLL      | Unroll factor, 1 to 16, of the loops over conformations    |
   +----------------------+------------------------------------------------------------+
   |                      | of the explicit SIMD routines for pair interactions, such  |
   +----------------------+------------------------------------------------------------+
   |                      | as Coulomb, van der Waals and bonds. Empty (default) uses  |
   +----------------------+------------------------------------------------------------+
   |                      | 4 for each SIMD. Changing it compiles only one file again  |
   +----------------------+------------------------------------------------------------+
   |                      | and is meant for finding the fastest value on a computer,  |
   +----------------------+------------------------------------------------------------+
   |                      | see below.                                                 |
   +----------------------+------------------------------------------------------------+
   | ACT_OPENMP_SIMD      | Let the compiler vectorize the loops over conformations    |
   +----------------------+------------------------------------------------------------+
   |                      | marked with OpenMP SIMD directives, using -fopenmp-simd.   |
   +----------------------+------------------------------------------------------------+
   |                      | This does not turn on OpenMP multithreading. Default ON.   |
   +----------------------+------------------------------------------------------------+
   | ACT_VECTOR_MATH      | Use a vector math library for mathematical functions such  |
   +----------------------+------------------------------------------------------------+
   |                      | as exp and erf in the loops vectorized by the compiler:    |
   +----------------------+------------------------------------------------------------+
   |                      | the system math library on macOS, libmvec or SVML          |
   +----------------------+------------------------------------------------------------+
   |                      | otherwise. Without such a library these loops are not      |
   +----------------------+------------------------------------------------------------+
   |                      | vectorized. Default ON.                                    |
   +----------------------+------------------------------------------------------------+
   | ACT_COUNT_FLOPS      | Count the floating point operations, such that the timings |
   +----------------------+------------------------------------------------------------+
   |                      | at the end of the log file report the GFlop/s of the parts |
   +----------------------+------------------------------------------------------------+
   |                      | of the program and the fraction of the peak performance.   |
   +----------------------+------------------------------------------------------------+
   |                      | Default OFF, since counting costs up to 10% of the time of |
   +----------------------+------------------------------------------------------------+
   |                      | molecular dynamics. Meant for developers.                  |
   +----------------------+------------------------------------------------------------+
   | ACT_GPU              | GPU support: NONE (default) or SYCL, which needs           |
   +----------------------+------------------------------------------------------------+
   |                      | AdaptiveCpp, see below.                                    |
   +----------------------+------------------------------------------------------------+

The SIMD (single instruction, multiple data) instructions of the processor are used to compute the energies and forces of several conformations of a molecule at once, which speeds up training considerably.
At the start of each program the log file reports the SIMD of the processor, the SIMD the compiler could use for the rest of the code and whether explicit SIMD is used, for instance::

  Explicit SIMD:      AVX2, 4 conformations at once, pair loops unrolled 4 times, in use

For most users the defaults are best. To build for a cluster with processors different from the one where you compile, set *ACT_SIMD* to the SIMD of the compute nodes. Since the explicit SIMD routines are only used if the processor supports them, a program built with AVX_512 still runs correctly on a processor with only AVX2, albeit more slowly. The rest of the code is vectorized by the compiler for the instructions given in *CMAKE_CXX_FLAGS*, for instance *-march=native* to use all instructions of the processor where you compile.
Setting the environment variable *ACT_NO_EXPLICIT_SIMD* when running a program turns off the explicit SIMD routines, which can be used to check the results or the speed of the different routines::

  ACT_NO_EXPLICIT_SIMD=1 alexandria train_ff ...

The log file also reports the peak performance of a core, estimated from the SIMD instructions of the processor (floating point operations per cycle) and its maximum clock frequency, which the timings at the end of the log file use to report the efficiency of the potentials, if ACT was configured with *-DACT_COUNT_FLOPS=ON*. For a processor that ACT does not recognize, the environment variable *ACT_PEAK_GFLOPS* sets the peak of a core in GFlop/s in double precision.

The loops over conformations of the explicit SIMD routines for pair interactions are unrolled, such that the processor has the instructions for several SIMD vectors of conformations at once. Which unroll factor is fastest depends on the processor, since a larger factor helps when the loop overhead matters, but helps less with wide SIMD such as AVX-512 and not when only few conformations are computed, as in the last iterations of relaxing shells. To find the best value on your computer, configure with *-DACT_COUNT_FLOPS=ON*, build with several values of *ACT_PAIR_UNROLL* in the same build directory and compare the GFlop/s of the potentials, such as *Potential COULOMB_POINT*, in the timings at the end of the log file of, for instance, a calculation of the second virial coefficient::

  cmake -DACT_PAIR_UNROLL=2 . && make -j 8 install
  cmake -DACT_PAIR_UNROLL=8 . && make -j 8 install
  cmake -DACT_PAIR_UNROLL= . && make -j 8 install

The last command returns to the default. The results do not depend on the unroll factor.

GPU support with SYCL
~~~~~~~~~~~~~~~~~~~~~
GPU support uses SYCL with `AdaptiveCpp <https://adaptivecpp.github.io>`_, which compiles the GPU code at run time for the GPU found, such that the same build runs on NVIDIA and AMD GPUs. AdaptiveCpp is available from conda-forge, in its own environment, for instance for NVIDIA GPUs::

  conda create -n act-sycl -c conda-forge "adaptivecpp=*=cuda*"
  cmake .. -DACT_GPU=SYCL -DAdaptiveCpp_DIR=$HOME/miniforge3/envs/act-sycl/lib/cmake/AdaptiveCpp

The variant matters: AdaptiveCpp uses a GPU through a backend, CUDA for NVIDIA GPUs and HIP for AMD GPUs, and the *cpu* variants of conda-forge, which other packages may pull into an environment, lack the CUDA backend. *cmake* reports the backends of the AdaptiveCpp it uses, which the line *GPU support* at the start of the log file repeats, and warns if the CUDA backend is missing. If a program finds no GPU, the log file says why, as far as it can tell: for instance that the node has NVIDIA GPUs but AdaptiveCpp has no CUDA backend, that the queue system gave the job no GPU, or that the node has no GPUs. The command *acpp-info -l* of AdaptiveCpp lists the devices it finds.

Only the GPU code is compiled with AdaptiveCpp's compiler acpp, the rest of ACT with the C++ compiler, which needs to be GCC 13 or newer, since the libraries of AdaptiveCpp need its C++ library. The GPU is used in single precision builds only (*-DACT_DOUBLE=OFF*). The section *Code used* at the start of the log file lists the GPUs found, with an estimate of their peak performance in single precision: the number of compute units times the single precision lanes of a compute unit, which depends on the type of GPU, times two operations per fused multiply-add, times the maximum clock frequency. The processes on a node are spread over its GPUs, by their number on the node, which Slurm or the MPI library sets.

In training with a polarizable force field, by default (*train_ff -gpu_resident*), each process keeps the structures of its molecules on the GPU, and computes there in each evaluation the charges, if they are generated with SQE, the virtual sites, all interactions and the relaxation of the shells. Only the parameters go to the GPU. In the evaluations of training the GPU also computes the deviations of the energies from the reference data, such that only those come back; the energies of the molecules themselves are computed for the output at the end. If the GPU does not support the potentials of the force field, the processor computes everything, as with *-nogpu_resident*. Otherwise the GPU computes the Coulomb interactions of several structures of all molecules at once and relaxes the shells of all molecules at once: the interactions with shells, which may be Coulomb, polarization, Lennard-Jones 12-6, buffered 14-7, Wang-Buckingham and MacDaniel-Schmidt, are computed on the GPU until the shells converged, the other interactions once on the processor. Since starting the GPU takes some tens of microseconds, this is only done for at least 100000 interactions times structures of all molecules, which the environment variable *ACT_GPU_MIN_WORK* changes.

Many processes per GPU keep it busy, since each process waits for the GPU part of the time. When several MPI ranks share a GPU, start the NVIDIA multi-process service in the job script, otherwise the GPU switches between the processes, which makes each operation on it about twenty times slower. The service takes at most 48 or 60 processes per GPU, depending on the version of the NVIDIA driver (48 on an RTX A4500 with driver 550, while NVIDIA documents 60 for current drivers), and the ranks beyond that compute without the GPU. *CUDA_DEVICE_MAX_CONNECTIONS* larger than its default of 2 lowers that number::

  export CUDA_MPS_PIPE_DIRECTORY=/tmp/mps-$SLURM_JOB_ID CUDA_MPS_LOG_DIRECTORY=/tmp/mpslog-$SLURM_JOB_ID
  mkdir -p $CUDA_MPS_PIPE_DIRECTORY $CUDA_MPS_LOG_DIRECTORY
  nvidia-cuda-mps-control -d
  srun ... alexandria train_ff ...
  echo quit | nvidia-cuda-mps-control

In *train_ff* with *-nogpu_master* the master process does not use the GPU, which leaves a place among the processes of the GPU and some hundreds of MB of its memory to the other processes. The master then computes the evaluations for the output, at the start and the end, on the processor, which took about 40 s more for the BFdb benchmark on one core.

Each process uses one GPU of the node, chosen by its number on the node, and ACT lets it see only that one: if *CUDA_VISIBLE_DEVICES*, *ROCR_VISIBLE_DEVICES* or *HIP_VISIBLE_DEVICES* lists several GPUs, as Slurm sets them, only that of the process is kept, and without such a list on a node with several NVIDIA GPUs *CUDA_VISIBLE_DEVICES* is set to it. Otherwise the GPU runtime would make a context on each GPU in each process, and with several GPUs per node the multi-process service would run out of processes. The log file tells which GPU each process uses.

With a population larger than the number of processes, each process can take several individuals at once with *train_ff -lockstep*. The GPU then runs the MCMC steps of these individuals together, in the same calculations, such that fewer processes use it better. The result of each individual is the same as when it is computed by itself, except for where the relaxation of the shells starts, see *-shell_keep*. With a population of 128, *-lockstep 8* with 17 processes was the fastest on a node with eight AMD MI250X GCDs and on one with four NVIDIA GH200, 1.5 to 4 times faster than one individual per process with 129 processes; see the chapter on performance for the measurements.

The GPU gives the same results as the processor, within single precision. The option *-nogpu*, which all programs have, turns the GPU off, for instance to compare the run time with and without the GPU. The programs set *ACPP_DEBUG_LEVEL* to 1, unless it is set already, such that AdaptiveCpp does not warn about the GPU backends it cannot load. Set it to 2 to see the messages of AdaptiveCpp.

============================
Troubleshooting Installation
============================
Sometimes, the *cmake* process or building using *make* does not work as described above.
For instance, there may be library mismatches like this::

   undefined reference to `__cxa_call_terminate@CXXABI_1.3.15'
   
which is caused by the fact that conda libraries expect specific versions of system libraries, combined with the make process supplying a wrong version of that library.
In that case you can instruct cmake to change the order of libraries by adding the flags to the cmake command line according to :numref:`tab-cmake`.


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


Testing the GPU code
~~~~~~~~~~~~~~~~~~~~
The tests of the GPU code need a GPU, in a single precision build with *-DACT_GPU=SYCL*. Their names end in *Gpu*, and *ctest* runs them separately from the other tests of a test program: *ACTForcesTestsGpu* (the pair potentials on the GPU), *AlexandriaTestsGpu* (the molecules on the GPU, several sets of parameters at once, the deviations on the GPU and *-shell_keep*) and *AlexandriaMpiTestsGpu* (training on the GPU with *-lockstep* and with molecules that are skipped). Without a GPU, or in double precision, *ctest* reports these three as skipped, for instance::

  20/29 Test #20: ACTForcesTestsGpu ................***Skipped   0.85 sec

and lists them at the end under *The following tests did not run*. A test run on a computer without a GPU, such as the head node of a cluster, therefore does not test the GPU code, even when all tests pass. *ACTGpuTests* checks the estimate of the peak performance and the explanation of why no GPU was found, which need no GPU. If a GPU is found but cannot compute the molecules of a test, that test fails.

On a cluster, compile everything on the head node, the tests included, and run the tests in a job on a node with a GPU. The head node needs no GPU to compile, since AdaptiveCpp compiles the GPU code when a program first uses the GPU. Set *ACT_SIMD* to the SIMD instructions of the GPU nodes if their processors differ from that of the head node, see above. If the processors of the GPU nodes are of another type than that of the head node, for instance ARM against x86, compile in a job on a GPU node instead. On csb, for instance, after::

  cd ACT/build_Release_SINGLE_SYCL
  make -j 20 alexandria tests

on the head node, a job script like this runs the tests on a GPU node, with the modules that were loaded for compiling::

  #!/bin/bash -l
  #SBATCH -p CLUSTER-AMD
  #SBATCH -N 1 -n 8
  #SBATCH --gres=gpu:1
  #SBATCH -t 00:40:00
  ml gcc/14 fftw openmpi/5.0.5
  cd $HOME/ACT/build_Release_SINGLE_SYCL
  $HOME/miniforge3/envs/ACT/bin/ctest -j 4 -V > ctest-gpu.log 2>&1
  grep -E "tests passed|Skipped|\*\*\*" ctest-gpu.log

The program *ctest* comes with *cmake*. On csb the *cmake* module does not work on the compute nodes, therefore the script uses the *ctest* of the conda environment of ACT, which has *cmake* as well. The last line should report that all tests passed, and no test should be skipped. Run *make tests* again after each change before submitting the job, since *make* by itself does not rebuild the test programs. To run a single test program on a GPU node interactively, for instance::

  srun -p CLUSTER-AMD --gres=gpu:1 -n 1 -c 8 bin/act-forces-test --gtest_filter='*Gpu'

The first run on a type of GPU takes longer, since AdaptiveCpp then compiles the GPU code, which it keeps in *~/.acpp* for later runs. If the MPI tests fail in the job because MPI programs are started differently there, run the others with *ctest -LE MpiTest*. Most tests have a time limit of 120 s, since some of them take 40 s on older processors.
