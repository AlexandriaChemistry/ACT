#!/usr/bin/env python3

import os

# All these files should be in place or alexandria will crash.
sel    = "../SELECTIONS/alcohol.dat"
gentop = "../ACS-pg.dat"
xml    = "../XML/alcohol.dat"

os.system("mpirun -np 2 alexandria tune_eem -v -d %s -sel %s -f %s -fit 'alpha zeta' -maxiter 500 -temp 10" % ( gentop, sel, xml ))
