#!/usr/bin/env python3

import os

# All these files should be in place or alexandria will crash.
sel    = "../SELECTIONS/alcohol.dat"
gentop = "../ACS-pg.xml"
xml    = "../XML/alcohol.xml"

os.system("mpirun -np 2 alexandria tune_eem -v -fc_mu 1  -d %s -sel %s -f %s -fit 'alpha zeta' -maxiter 500 -temp 200" % ( gentop, sel, xml ))
