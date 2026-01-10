#
# This file is part of the Alexandria Chemistry Toolkit
#
# Copyright (c) 2014-2026, by the Alexandria Chemistry team, led by
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
         atom_bond.xml
         atomtype_mapping.csv
         atomtypes_gui.csv
         excluded_molecules.csv
         residuetypes.dat 
         symmetric_charges.csv
         vsite2_halogen_sigma_holes.csv
         vsite2_water_sigma_holes.csv
         vsite3_water_bisector.csv
         vsite3_water_lone_pairs.csv
         vsite3_alcohol_bisector.csv
         vsite3_alcohol_lone_pairs.csv
         vsite4_amine.csv
    DESTINATION ${ACT_INSTALL_ACTDATADIR}
    COMPONENT data)
