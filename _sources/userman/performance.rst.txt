.. _sec-performance:

*************************************
Getting the best performance from ACT
*************************************

Training a force field takes most of the computer time spent with ACT, so this chapter is mainly about *alexandria train_ff*. The run times quoted are examples from our own tests. The best settings depend on the computer and on the training set, so the most useful advice is to do a short test run with a few settings before a long one, see :ref:`sec-measuring`.

======================
Building ACT for speed
======================
* Build in single precision, with *-DACT_DOUBLE=OFF*. The explicit SIMD routines then compute twice as many interactions per instruction, and a GPU can only be used in single precision builds. Energies are still summed in double precision.
* On a cluster, set *ACT_SIMD* to the SIMD instructions of the compute nodes if you compile on a head node with another processor, and give the compiler the instructions of the compute nodes in *CMAKE_CXX_FLAGS*, for instance *-march=native* when compiling on a compute node.
* For a GPU, build with *-DACT_GPU=SYCL*, see the installation chapter.

The section *Code used* at the start of each log file tells which SIMD instructions, which precision and which GPU the program uses. Check it once on each new computer.

===============================
How train_ff uses the processes
===============================
One process, the master, runs the genetic algorithm. How the others work depends on the size of the population (*-pop_size*):

* **Population no larger than the number of processes.** The processes are divided into one group per individual, a middleman and its helpers, and the molecules of the training set are divided over the processes of a group. Each individual needs the same number of processes, so use a multiple of the population size; otherwise the processes left over are idle, which *train_ff* warns about. Splitting the molecules over helpers stays efficient down to 4-5 molecules per process :cite:p:`Spoel2025b`.
* **Population larger than the number of processes.** This is possible with the GA and HYBRID optimizers. Each process except the master is then a middleman that computes one individual at a time, and gets the next one from the master when it is done. The master itself computes no individuals.

In the second case, a generation ends when the last individual is done. If the number of middlemen does not divide the population, the last round of individuals keeps only some of the middlemen busy, and the others wait. Choose the number of processes such that every middleman computes the same number of individuals: one more than a divisor of the population. For example, a training with a population of 256 on a node with 288 cores (arrhenius, GH200) took:

.. table:: Run time of a training with a population of 256, by the number of processes.
   :name: tab-processes

   +-----------+------------+---------------------------+--------------+
   | Processes | Middlemen  | Individuals per middleman | Run time (s) |
   +===========+============+===========================+==============+
   | 64        | 63         | 4 or 5                    | 89.1         |
   +-----------+------------+---------------------------+--------------+
   | 65        | 64         | 4                         | 82.0         |
   +-----------+------------+---------------------------+--------------+
   | 129       | 128        | 2                         | 40.9         |
   +-----------+------------+---------------------------+--------------+

With 64 processes, 4 middlemen compute a fifth individual in each generation while the 59 others wait. All three runs gave the same result, which does not depend on the number of processes. In general, use as many processes as there are physical cores, one more than a divisor of the population. At the start, *train_ff* notes on the terminal and in the log file when the middlemen would share the individuals unevenly, how long they would wait on average, and which numbers of processes avoid it.

=================
Processor or GPU?
=================
The GPU computes the training only for **polarizable** force fields, in single precision builds, with *train_ff -gpu_resident* (the default). Each process then keeps its molecules on the GPU, which computes the interactions, the virtual sites, the relaxation of the shells and the deviations from the reference data, and the charges if they are generated with SQE. The relaxation of the shells repeats the force calculation many times in each evaluation, which is expensive on a processor, and this is where the GPU gains most.

For force fields **without shells**, the processor computes everything, also in a build with GPU support. Such a training takes the same time with and without *-nogpu*. Without shells, ACT only recomputes the molecules that a changed parameter affects, and of a dimer only the interactions between its two molecules, so one evaluation is short: about a hundred microseconds on a core for a training set of small dimers in our tests. That is shorter than the time needed just to start the work on a GPU and wait for the result, so for these force fields, use all cores of the node rather than the GPUs.

================
Training on GPUs
================
* **Several processes per GPU.** Each process waits for the GPU part of the time, so several processes per GPU keep it busy. On NVIDIA GPUs, start the multi-process service in the job script when several processes share a GPU, see the installation chapter. The service takes at most 48 or 60 processes per GPU, depending on the NVIDIA driver; the ranks beyond that compute without the GPU.
* **-lockstep.** With a population larger than the number of processes, *-lockstep k* gives each middleman *k* individuals at once. The GPU runs the MCMC steps of these individuals together, in the same calculations: each step takes as many operations on the GPU, and as much waiting for it, for *k* individuals as for one. Fewer processes are needed, *pop_size/k* + 1, which also shortens the start-up of the GPU in each process. With a population of 128 on a node with eight AMD MI250X GCDs (dardel) and on a node with four NVIDIA GH200 (arrhenius) we measured:

  .. table:: Wall time of the master in s of short trainings with a population of 128, by *-lockstep*, with *-shell_keep*.
     :name: tab-lockstep

     +---------------------------+------------+----------------+-----------+-----------+------------+
     | Computer                  | Training   | *-lockstep 1*  | *4*       | *8*       | *16*       |
     |                           |            | 129 processes  | 33        | 17        | 9          |
     +===========================+============+================+===========+===========+============+
     | dardel, 8 MI250X GCDs     | HHpol      | 63.3           | 32.4      | **21.0**  | 28.5       |
     +---------------------------+------------+----------------+-----------+-----------+------------+
     |                           | BFdb       | 243.1          | 167.6     | **156.4** | 172.8      |
     +---------------------------+------------+----------------+-----------+-----------+------------+
     | arrhenius, 4 GH200        | HHpol      | 101.0          | 34.8      | 23.5      | **22.0**   |
     +---------------------------+------------+----------------+-----------+-----------+------------+
     |                           | BFdb       | 156.5          | 84.6      | 76.2      | **75.5**   |
     +---------------------------+------------+----------------+-----------+-----------+------------+

  HHpol is a small training set (10 generations of 10 MCMC iterations), and most of what it gains from *-lockstep 8* is start-up time. BFdb is a large one (5 generations of 5 iterations); its GPUs are busy already with *-lockstep 4*, and computing more individuals at once does not make them faster. With 9 processes on 8 GCDs the GPUs have too little to do. *-lockstep 8* is a good choice on both computers. The result does not depend on *-lockstep*, except for where the relaxation of the shells starts, see *-shell_keep*, and through which middleman computes which individual: with the same seed, the BFdb runs above gave the same force field for all values, the HHpol runs did not.

  Without the GPU, or for a force field without shells, *train_ff* notes that *-lockstep* has no effect and hands out one individual at a time, which shares them more evenly over the middlemen: in a training with 256 individuals and 64 processes on dardel this took 1062 s, against 1601 s with groups of 4 individuals, since then one middleman computed a second group in each generation while the others waited. With *-lockstep*, the number of middlemen should divide the number of groups, see the note at the start. Each individual that a process computes at once also needs a copy of the molecules on the GPU, about 43 MB for BFdb; if the GPU runs out of memory, *train_ff* stops with a message that says so.
* **-shell_keep.** By default the relaxation of the shells starts from where the shells were the previous time, which needs fewer iterations when the parameters changed little. On an RTX A4500 this took a training from 179 s to 146 s. The relaxed shells differ within the tolerance on the forces, so the results differ slightly from those of *-noshell_keep*. In the trainings of :numref:`tab-lockstep`, which used *-shell_keep*, the shells needed about 25% fewer iterations than with *-noshell_keep*; BFdb took 7% less time on dardel (156 s against 167 s with *-lockstep 8*), and on arrhenius, where the relaxation is a smaller part of the time, the difference was 1-2%. Since where the relaxation starts depends on which individuals a process computed before, the result of a training with *-shell_keep* depends on the number of processes and on the computer. With *-noshell_keep*, the HHpol trainings reached the same :math:`\chi^2` on dardel and on arrhenius for each value of *-lockstep*, and with *-lockstep 4* the same force field; with *-shell_keep* they did not. Both are correct within the tolerance on the forces; use *-noshell_keep* when runs must be reproduced exactly.
* Molecules that a changed parameter does not affect are not computed again; this happens automatically, also on the GPU.

The best value of *-lockstep* depends on the GPU, on the number of GPUs per node and on the training set. Start with *-lockstep 8* and *pop_size/8* + 1 processes, and compare with 4 and 16 in a short run.

.. _sec-measuring:

====================================
Measuring and comparing the run time
====================================
At the end of the log file, *train_ff* prints how much time each part of the program took, summed over all processes and as an average per process for the master and for the middlemen, preceded by the wall time of the master. Useful lines are:

* *MPI waiting in Middleman loop*: middlemen waiting for work, for instance at the end of each generation. A large share means that the processes are not used evenly, see :numref:`tab-processes`.
* *MPI waiting in Genetic algorithm*: the master waiting for the middlemen, which is normal.
* *GPU initialization* and *Startup*: one-time costs. In short test runs they can be a large part of the time; leave them out when comparing settings.
* *GPU resident relaxation*, *GPU resident evaluation* and *Charge generation*: the computation itself.

When comparing settings, keep in mind that:

* The run time of the same training on a GPU can vary by about 10% between runs. Compare settings in one job, run them alternately, and repeat them.
* The training may stop before *-max_generations* when the fitness of the test set no longer improves (*-max_test_generations*). Runs that differ in the last digits, for instance on different computers, may then run a different number of generations. Compare the time per generation, from the lines *Generation n/N* in the log file.
* *Wall time of the master* starts after MPI has started. On some clusters starting MPI takes long, for instance about 0.7 s per process on the GH200 nodes of arrhenius, so the job takes longer than this time.
