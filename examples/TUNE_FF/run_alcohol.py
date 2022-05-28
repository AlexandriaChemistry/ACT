#!/usr/bin/env python3

#SBATCH -n 32

import os
from act import *

# All these files should be in place or alexandria will crash.
sel    = "../SELECTIONS/alcohol.dat"
xml    = "../XML/alcohol.xml"
act = ACT(xml, sel, True)

ForceFieldFileIn  = "../ACS-pg.xml"
act.bastat(ForceFieldFileIn, ForceFieldFileIn, "bastat.log", {})    
for target in Target:
    ForceFieldFileOut = ( "tune_ff_%s.xml" % ( target.name ) )
    LogFile           = ( "tune_ff_%s.log" % ( target.name ) )
    EpotFile          = ( "epot_%s.xvg" % ( target.name ) )
    ConvFile          = ( "conv_%s.xvg" % ( target.name ) )
    act.tune_ff(ForceFieldFileIn, ForceFieldFileOut, 
                LogFile, target, 
                { "-max_generations": 20, 
                  "-pop_size":        32,
                  "-epot":            EpotFile,
                  "-conv":            ConvFile })
    ForceFieldFileIn  = "Train-" + ForceFieldFileOut
