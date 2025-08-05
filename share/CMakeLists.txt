#
# This file is part of the Alexandria Chemistry Toolkit
#
# Copyright (c) 2014-2025, by the Alexandria Chemistry team, led by
# David van der Spoel, Paul J. van Maaren and Mohamman M. Ghahremanpour. 
#
# See License file in top level directory for license info. 
#
install(FILES 
         ACTlogo2.png
         alexandria.csv
         atomization-energies-dft.csv
         atomization-energies.csv
         atomprops.csv
         atomprops.csv
         atomtype_mapping.csv
         atomtypes_gui.csv
         bondvsites.csv
         excluded_molecules.csv
         residuetypes.dat 
         symmetric_charges.csv
    DESTINATION ${ACT_INSTALL_ACTDATADIR}
    COMPONENT data)
