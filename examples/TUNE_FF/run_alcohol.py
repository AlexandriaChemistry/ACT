#!/usr/bin/env python3

#SBATCH -t "24:00:00"
#SBATCH -n 32

import os
from act import *

# All these files should be in place or alexandria will crash.
sel    = "../SELECTIONS/alcohol.dat"
xml    = "../XML/alcohol.xml"
act = ACT(xml, sel, True)

# The force field file we start with
ForceFieldFileIn  = "../ACS-pg.xml"

# Generate Bonds, Angles etc.
act.bastat(ForceFieldFileIn, ForceFieldFileIn, "bastat.log", {})

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
    EpotFile          = ( "epot_%s.xvg" % ( target.name ) )
    ConvFile          = ( "conv_%s.xvg" % ( target.name ) )
    options           = { "-max_generations": 20, 
                          "-pop_size":        32,
                          "-epot":            EpotFile,
                          "-conv":            ConvFile }
    # After the last step of the optimizations, where we tune the
    # frequencies, we will also print some files for analysis.
    if Target.Freq == target:
        options["-printSP"]   = ""
        options["-alphacorr"] = "alpha_corr.xvg"
        options["-espcorr"]   = "esp_corr.xvg"
        options["-epotcorr"]  = "epot_corr.xvg"
        options["-freqcorr"]  = "freq_corr.xvg"
    else:
        options["-nocalc_frequencies"] = ""
    act.tune_ff(ForceFieldFileIn, ForceFieldFileOut, 
                LogFile, target, options)
    ForceFieldFileIn  = "Train-" + ForceFieldFileOut
