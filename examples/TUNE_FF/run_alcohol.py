#!/usr/bin/env python3

#SBATCH -t "24:00:00"
#SBATCH -n 16

import os
from act import *

# All these files should be in place or alexandria will crash.
sel    = "../SELECTIONS/alcohol.dat"
xml    = "../XML/alcohol.xml"
act = ACT(xml, sel, True)

# The force field file we start with
ForceFieldFileIn  = "../ACS-pg.xml"

# Generate Bonds, Angles etc.
act.bastat(ForceFieldFileIn, ForceFieldFileIn, "bastat.log", 
           { "-klin": 36000 })

# Now loop over the optimization targets.
# First, the EEM parameters will be optimized to reproduce the
# electrostatic potential.
# Second, the parameters for energy and forces will be tuned
# to get those properties correct.
# Third, another round of optimization where also the frequencies
# are targeted.
for target in Target:
    ForceFieldFileOut = ( "tune_ff_%s.xml" % ( target.name ) )
    LogFile           = ( "tune_ff_%s.log" % ( target.name ) )
    Chi2File          = ( "chi2_%s.xvg" % ( target.name ) )
    ConvFile          = ( "conv_%s.xvg" % ( target.name ) )
    # Make sure the pop_size is even and at most as large as the
    # number of core on the machine. 
    options           = { "-max_generations": 20, 
                          "-pop_size":        16,
                          "-chi2":            Chi2File,
                          "-conv":            ConvFile }
    # After the last step of the optimizations, where we tune the
    # frequencies, we will also print some files for analysis.
    if Target.Freq == target:
        options["-norandom_init"] = ""
        options["-printSP"]   = ""
        options["-alphacorr"] = "alpha_corr.xvg"
        options["-espcorr"]   = "esp_corr.xvg"
        options["-epotcorr"]  = "epot_corr.xvg"
        options["-freqcorr"]  = "freq_corr.xvg"
    else:
        options["-random_init"] = ""
        options["-nocalc_frequencies"] = ""
    act.tune_ff(ForceFieldFileIn, ForceFieldFileOut, 
                LogFile, target, False, options)
    ForceFieldFileIn  = "Train-" + ForceFieldFileOut
