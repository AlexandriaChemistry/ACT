#!/bin/sh

#SBATCH -n 16
#SBATCH -t 24:00:00

# Remove oversubscribe option if your mpirun does not support it
mpirun -n 16 -oversubscribe alexandria train_ff -ff Train-inter.xml -mp ../XML/alcohol-mp2.xml -sel ../SELECTIONS/alcohol-monomer.dat -optimizer HYBRID -max_generations 10 -max_iter 20 -temp 1 -anneal_globally -anneal 0  -printSP -pop_size 16 -fit ' De  bondlength  bondenergy  beta  kb  angle  kt  c0  c1  c2  c3  ' -fc_epot 1 -fc_force 0.1 -n_crossovers 2 -cp_gen_interval 5 -cp_pop_frac 0.2 -fc_bound 1000  -v 4 -o intra.xml -g train_intra.log

