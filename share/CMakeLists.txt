#
# This file is part of the Alexandria Chemistry Toolkit
#
# Copyright (c) 2014-2026, by the Alexandria Chemistry team, led by
# David van der Spoel, Paul J. van Maaren and Mohamman M. Ghahremanpour. 
#
# See License file in top level directory for license info. 
#
# The other data files, which the Python tools use as well, are in
# ACT-utilities, which installs them in share/act too.
install(FILES
         residuetypes.dat
         symmetric_charges.csv
    DESTINATION ${ACT_INSTALL_ACTDATADIR}
    COMPONENT data)
