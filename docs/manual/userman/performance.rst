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

* **Population no larger than the number of processes.** The processes are divided into one group per individual, a middleman and its helpers, and the molecules of the training set are divided over the processes of a group, see :numref:`fig-trainff-helpers`. The master is the middleman of the first group. Each individual needs the same number of processes, so use a multiple of the population size; otherwise the processes left over are idle, which *train_ff* warns about. Splitting the molecules over helpers stays efficient down to 4-5 molecules per process :cite:p:`Spoel2025b`.
* **Population larger than the number of processes.** This is possible with the GA and HYBRID optimizers. Each process except the master is then a middleman that computes one individual at a time, and gets the next one from a queue at the master when it is done, see :numref:`fig-trainff-queue`. The master itself computes no individuals. A middleman holds all molecules of the training set.

.. figure:: ../images/trainff-helpers.pdf
   :name: fig-trainff-helpers
   :width: 100%
   :align: center

   A population of 4 individuals on 12 processes. Each individual is computed by a middleman and two helpers, which each compute a third of the molecules and send their part of the :math:`\chi^2` to the middleman. The master is the middleman of individual A. All individuals are computed at the same time.

.. figure:: ../images/trainff-queue.pdf
   :name: fig-trainff-queue
   :width: 100%
   :align: center

   A population of 8 individuals on 5 processes. The master keeps the individuals of a generation in a queue and gives the next one to the first middleman that is done. Each of the 4 middlemen computes 2 individuals per generation, one per round.

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
* **-lockstep.** With a population larger than the number of processes, *-lockstep k* gives each middleman *k* individuals at once. The GPU runs the MCMC steps of these individuals together, in the same calculations: each step takes as many operations on the GPU, and as much waiting for it, for *k* individuals as for one. Fewer processes are needed, *pop_size/k* + 1, which also shortens the start-up of the GPU in each process. The master then hands out groups of *k* individuals, and a middleman keeps a copy of the molecules on the GPU for each individual of its group, a replica, see :numref:`fig-trainff-lockstep`. The replicas share the topology; each one adds its own parameters, coordinates and shells, about 43 MB for BFdb, on top of about 330 MB for the GPU context of the process and about 270 MB for the first copy of the molecules. With a population *P*, lockstep *k* and *M* middlemen (the processes minus the master), a generation has *P/k* groups, each middleman computes *P/(k M)* of them one after another, in rounds, and the GPU computes *k M* individuals at once. For example, *P* = 512 with 4 middlemen on one GPU and *-lockstep 16* gives 32 groups and 8 rounds per middleman, with 64 individuals on the GPU at once and about 1.3 GB of GPU memory per process for BFdb. *-lockstep 128* would do the same work in one round, but needs about 6 GB per process, which no longer fits 4 processes on a GPU with 21 GB. With a population of 128 on a node with eight AMD MI250X GCDs (dardel) and on a node with four NVIDIA GH200 (arrhenius) we measured:

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

  HHpol is a small training set (10 generations of 10 MCMC iterations), and most of what it gains from *-lockstep 8* is start-up time. BFdb is a large one (5 generations of 5 iterations); its GPUs are busy already with *-lockstep 4*, and computing more individuals at once does not make them faster. With 9 processes on 8 GCDs the GPUs have too little to do. *-lockstep 8* is a good choice on both computers. With *-noshell_keep*, each individual gets the same result as when it is computed by itself; with *-shell_keep* it may differ slightly, see there. With the same seed, the BFdb runs above gave the same force field for all values of *-lockstep*; the HHpol runs, made before each individual started from its own shells, did not.

  Without the GPU, or for a force field without shells, *train_ff* notes that *-lockstep* has no effect and hands out one individual at a time, which shares them more evenly over the middlemen: in a training with 256 individuals and 64 processes on dardel this took 1062 s, against 1601 s with groups of 4 individuals, since then one middleman computed a second group in each generation while the others waited. With *-lockstep*, the number of middlemen should divide the number of groups, see the note at the start. Each individual that a process computes at once also needs a copy of the molecules on the GPU, about 43 MB for BFdb; if the GPU runs out of memory, *train_ff* stops with a message that says so.
* **-shell_keep.** By default the relaxation of the shells starts from where the shells were the previous time, which needs fewer iterations when the parameters changed little. On an RTX A4500 this took a training from 179 s to 146 s. The relaxed shells differ within the tolerance on the forces, so the results differ slightly from those of *-noshell_keep*. In the trainings of :numref:`tab-lockstep`, which used *-shell_keep*, the shells needed about 25% fewer iterations than with *-noshell_keep*; BFdb took 7% less time on dardel (156 s against 167 s with *-lockstep 8*), and on arrhenius, where the relaxation is a smaller part of the time, the difference was 1-2%. With a population larger than the number of processes, the shells are kept only between the MCMC steps of an individual: each time a process gets an individual from the queue, its relaxation starts with the shells on their parents. The :math:`\chi^2` of an individual then does not depend on which individuals the process computed before, so the result does not depend on the number of processes: a PC+GS training with a population of 32 gave the same force field with 3 processes and *-lockstep 16* as with 5 processes and *-lockstep 8*, and took 160-180 s against 217 s with *-noshell_keep*. It can still depend on whether the GPU computes the individuals of a group together or one at a time, as with *-lockstep 1* or a group in parts (see below): the two leave out different molecules whose parameters did not change in an MCMC step, and a molecule that is computed again starts from its kept shells rather than reusing its energies. In the PC+GS training, *-lockstep 1* and a group in parts gave the same force field, which differed from that of *-lockstep 16*; with the environment variable ACT_NO_GPU_SKIP, which computes all molecules each time, all three were the same. It may still differ between computers. With *-noshell_keep*, the HHpol trainings reached the same :math:`\chi^2` on dardel and on arrhenius for each value of *-lockstep*, and with *-lockstep 4* the same force field. Both are correct within the tolerance on the forces; use *-noshell_keep* when runs must be reproduced exactly. At the end of a training, *train_ff* computes the :math:`\chi^2` of the best force field again. With *-shell_keep* this starts from other shells than the MCMC step that found the best force field, and with poor parameters, whose shells often do not converge, the result may differ much, since each conformation whose shells do not converge adds a penalty, see :ref:`sec-boundspenalty`; *train_ff* then warns and writes the force field. With *-noshell_keep* the two must be the same, and *train_ff* stops if they are not.
* Molecules that a changed parameter does not affect are not computed again; this happens automatically, also on the GPU.

The best value of *-lockstep* depends on the GPU, on the number of GPUs per node and on the training set. Start with *-lockstep 8* and *pop_size/8* + 1 processes, and compare with 4 and 16 in a short run.

With *-lockstep 0*, *train_ff* chooses the value itself. If the GPU trains the force field, which the master checks in the first force field before the processes get their tasks, it aims for at least two middlemen per GPU and groups of at most 16 individuals, with at most 40 middlemen per GPU, and uses the fewest middlemen that compute the groups in the same number of rounds. For a population of 128 this gives the settings that were the fastest in :numref:`tab-lockstep`: groups of 16 on 8 middlemen with four GPUs, and groups of 8 on 16 middlemen with eight. The processes left over are idle, with a warning that says how many to start; with too few processes, a warning says how many would compute the groups in one round. Without GPU training, each middleman computes one individual at a time with all processes. The GPUs counted are those that the job may use, as listed in *CUDA_VISIBLE_DEVICES*, *ROCR_VISIBLE_DEVICES* or *HIP_VISIBLE_DEVICES*, or else those of the node.

With any value of *-lockstep*, a middleman puts as many replicas on the GPU as fit into its share of the GPU memory: the memory of the GPU divided over the processes that use it, less about 350 MB for the context of each and 10% that is kept free. It measures the memory of one replica at the first evaluation. If fewer fit than the individuals of a group, *train_ff* notes it and the middleman computes the group in parts, each with its individuals together, with the same results with *-noshell_keep*. Fewer processes per GPU, or a smaller *-lockstep*, then let more individuals be computed at once.

.. figure:: ../images/trainff-lockstep.pdf
   :name: fig-trainff-lockstep
   :width: 100%
   :align: center

   A population of 16 individuals on 3 processes with *-lockstep 4*. The master hands out 4 groups of 4 individuals; each of the 2 middlemen computes 2 groups per generation, one per round, with the 4 individuals of a group as 4 replicas on the GPU, computed together. With 5 processes (4 middlemen), each middleman would compute one group, in one round. Below, the GPU memory of one middleman for BFdb with the 6P force field.

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
