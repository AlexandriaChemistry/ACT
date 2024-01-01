#!/usr/bin/env python3

#SBATCH -t "24:00:00"
#SBATCH -n 16

import os
from act import *

# All these files should be in place or alexandria will crash.
selmonomer  = "../SELECTIONS/alcohol.dat"
seldimer    = "../SELECTIONS/alcoholdimer.dat"
xmlesp      = "../XML/alcohol-esp.xml"
xmlepot     = "../XML/alcohol-epot.xml"
xmldimer    = "../XML/alcohol-dimer.xml"
act_run     = ACT(xml, sel, True)

# The force field file we start with
ForceFieldFileIn  = "../ACS-pg.xml"

# Generate Bonds, Angles etc.
act_run.geometry_ff(ForceFieldFileIn, ForceFieldFileIn, "bastat.log", 
                { "-klin": 36000 })

# What bonded parameters to fit determines on what we used to generate the initial force field
fitmorse   = "'beta De D0 bondlength'"
fitcubic   = "'bondenergy rmax bondlength kb'"
fitbondeds = fitmorse

fitvdw = "'sigma epsilon gamma'"

# All these files should be in place or alexandria will crash.
run_flags = { Target.EEM: { "sel": selmonomer,  "xml": xmlesp },
              Target.Epot: { "sel": selmonomer, "xml": xmlepot,  "extra": { "-fit": fitbondeds } },
              Target.Inter: { "sel": seldimer,  "xml": xmldimer, "extra": { "-fit": fitvdw } } }

# Now loop over the optimization targets.
# First, the EEM parameters will be optimized to reproduce the
# electrostatic potential.
# Second, the parameters for energy and forces will be tuned
# to get those properties correct.
# Third, another round of optimization where also the frequencies
# are targeted.
ffACM = "Train-tune_ff_ACM.xml"
for target in [ Target.ACM, Target.Epot, Target.Inter ]:
    # Check whether we can use the charges determined earlier
    useACM = False
    if Target.ACM != target:
        # Set the charges in the xml file based on the Alexandria Charge Model
        xmlout = "temp_" + run_flags[target]["xml"]
        useACM = act.set_charges(ffACM, xmlesp, run_flags[target]["xml"], xmlout)
        if useACM:
            run_flags[target]["xml"] = xmlout
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
    # If we already have charges in the xml file, let's use them.
    if useACM:
        options["-qqm"] = "ACM"
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
<<<<<<< HEAD:examples/TRAIN_FF/run_alcohol.py
<<<<<<< HEAD:examples/TRAIN_FF/run_alcohol.py
    act_run.train_ff(ForceFieldFileIn, ForceFieldFileOut, 
                     LogFile, target, False, options)
=======
=======
>>>>>>> 568d6a250cc0251089c550564f33013aa8e8db69:examples/TUNE_FF/run_alcohol.py
    act_run.tune_ff(ForceFieldFileIn, ForceFieldFileOut, 
                LogFile, target, False, options)
>>>>>>> eafa6844c (Renamed variable in run_alcohol.py):examples/TUNE_FF/run_alcohol.py
    ForceFieldFileIn  = "Train-" + ForceFieldFileOut
