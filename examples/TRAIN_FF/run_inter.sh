#!/bin/sh

#SBATCH -n 16
#SBATCH -t 24:00:00

# Remove oversubscribe option if your mpirun does not support it
mpirun -n 16 -oversubscribe alexandria train_ff -ff ../myff2.xml -mp ../XML/alcohol-sapt.xml -sel ../SELECTIONS/alcohol-dimer.dat -optimizer HYBRID -max_generations 10 -max_iter 20 -temp 1 -anneal_globally -anneal 0  -printSP -pop_size 16 -fit 'chi eta zeta sigma epsilon' -fc_inter 1 -n_crossovers 2 -cp_gen_interval 5 -cp_pop_frac 0.2 -fc_bound 1000  -fc_charge 20000 -v 4 -o inter.xml -g train_inter.log

