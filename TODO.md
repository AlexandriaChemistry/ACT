## Open questions in the ACT code
+ file: ./basics/identifier.h line: 64 : return a boolean instead of crashing
+ file: ./basics/identifier.h line: 137 check implementation
+ file: ./basics/identifier.cpp line: 187 only insert ones we do not already have
+ file: ./molprop/molprop.h line: 96 implement this
+ file: ./molprop/molprop.h line: 111 check and double check
+ file: ./molprop/molprop_xml.h line: 50 implement using a serialized protocol rather than reading the
+ file: ./molprop/molprop_xml.h line: 68 implement using a serialized protocol rather than reading the
+ file: ./molprop/topologyentry.h line: 207 : this is a fake bond order.
+ file: ./molprop/molprop_util.cpp line: 59 : implement alexandria id here.
+ file: ./molprop/molprop_xml.cpp line: 640 make this more rigorous and less ugly.
+ file: ./molprop/edit_mp.cpp line: 80 check whether this is needed.
+ file: ./molprop/molprop.cpp line: 257 : this ignore the possibility that there could be
+ file: ./molprop/molprop.cpp line: 401 make this check more rigorous, check for overlaps etc.
+ file: ./molprop/molprop_tables.h line: 64 transform iqm to enum
+ file: ./molprop/molprop_tables.h line: 83 transform ims to enum
+ file: ./molprop/molprop_tables.h line: 117 transform ims to enum
+ file: ./molprop/molprop_tables.h line: 118 introduce enum to switch between absolute and relative tolerance
+ file: ./molprop/molprop_tables.h line: 141 more explanation text
+ file: ./ga/penalizer.h line: 140 : if needed, make this penalizer print the indices of the individuals it is
+ file: ./ga/tests/genetic_algorithmtest.cpp line: 126 check return value
+ file: ./forcefield/forcefield_parameter.h line: 71 : check unit
+ file: ./forcefield/forcefield.cpp line: 215 : add check for number of interactions?
+ file: ./forcefield/forcefield.cpp line: 723 do not ignore return value
+ file: ./python/act_openmm.py line: 805 : update comment
+ file: ./python/act_openmm.py line: 1071 fetch this number from system.context
+ file: ./python/act_openmm.py line: 1072 : this needs the erf function!
+ file: ./python/act_openmm.py line: 1450 check whether this if statement should be flipped.
+ file: ./python/act_gct.py line: 99 move this check out of this loop
+ file: ./python/mol_csv_api.py.cmakein line: 115 : fix this to be correct.
+ file: ./python/act.py line: 41 implement algorithm to suggest reasonable pop_size
+ file: ./alexandria/fragmenthandler.cpp line: 162 : check whether this works for polarizable models
+ file: ./alexandria/fragmenthandler.cpp line: 205 only copy the coordinates if there is more than one fragment.
+ file: ./alexandria/train_ff.cpp line: 159 rename function and make a coupling between targets from the command line
+ file: ./alexandria/train_ff.cpp line: 208 : what about the flags? here it is a bit more clear that they should be all false?
+ file: ./alexandria/train_ff.cpp line: 242 only open these files when we are optimizing in verbose mode.
+ file: ./alexandria/train_ff.cpp line: 590 : resetting the train parameters for the trainffprinter. we may have to work on that if we want to show the best test parameters too
+ file: ./alexandria/train_ff.cpp line: 608 : parallellize this. fixme: there is no need to do this i believe, it's done above, and parallel!
+ file: ./alexandria/train_ff.cpp line: 696 : pargs is a confighandler, maybe we could inherit the superclass?
+ file: ./alexandria/train_utility.cpp line: 83 add error checking for lsq
+ file: ./alexandria/train_utility.cpp line: 1291 add checks for existence
+ file: ./alexandria/devcomputer.cpp line: 422 : compute this only once if both dipole and quadrupole are used in fitting
+ file: ./alexandria/devcomputer.cpp line: 634 double check if the atomizationenergy is needed.
+ file: ./alexandria/devcomputer.cpp line: 669 fix beta (but how?)
+ file: ./alexandria/gen_table.cpp line: 279 fix this
+ file: ./alexandria/babel_io.cpp line: 816 : check when we need to read just one molecule
+ file: ./alexandria/staticindividualinfo.cpp line: 163 fix the uncertainty
+ file: ./alexandria/staticindividualinfo.cpp line: 167 : this will generate the whole matrix of parameters from scratch.
+ file: ./alexandria/dissociation_energy.cpp line: 454 free the gmx_stats_t
+ file: ./alexandria/normalmodes.cpp line: 243 this will crash
+ file: ./alexandria/devcomputer.h line: 204 move to cpp file
+ file: ./alexandria/molhandler.cpp line: 493 : make ftoler a parameter.
+ file: ./alexandria/molhandler.cpp line: 839 check which one it should be
+ file: ./alexandria/molhandler.cpp line: 989 : is this really needed?
+ file: ./alexandria/molhandler.cpp line: 1222 : rewrite without gromacs fluff.
+ file: ./alexandria/allbondeds.cpp line: 515 make this a parameter
+ file: ./alexandria/acm_ga.cpp line: 160 : have we already checked that the number of processors is the correct one?
+ file: ./alexandria/acm_ga.cpp line: 235 : check whether we need to update this at all here
+ file: ./alexandria/acm_ga.cpp line: 272 : is this necessary? it will always be train?
+ file: ./alexandria/acm_ga.cpp line: 316 : can we remove the > 0?
+ file: ./alexandria/acm_ga.cpp line: 396 : is this necessary? it will always be train?
+ file: ./alexandria/acm_ga.cpp line: 405 : can we just negate instead of comparing?
+ file: ./alexandria/acm_ga.cpp line: 429 : if we end up sending more stuff, it might be worth it to just send the entire genome
+ file: ./alexandria/mcmcmutator.h line: 54 : shouldn't we move sensitivity analysis somewhere else?
+ file: ./alexandria/actmol.cpp line: 352 store the interaction forces
+ file: ./alexandria/actmol.cpp line: 694 is this correct? should not each qp have it's own coordinates?
+ file: ./alexandria/actmol.cpp line: 713 check whether this needed
+ file: ./alexandria/actmol.cpp line: 770 , likely we should not change the coordinates here, just the charges
+ file: ./alexandria/actmol.cpp line: 822 do not use -1 here, but particle.core()
+ file: ./alexandria/actmol.cpp line: 873 what about shells?
+ file: ./alexandria/actmol.cpp line: 896 not sure whether this is needed but why not.
+ file: ./alexandria/actmol.cpp line: 1096 write a replacement for this function
+ file: ./alexandria/actmol.cpp line: 1164 check but we should likely not update coords
+ file: ./alexandria/actmol.cpp line: 1398 check whether this is needed. likely it is here, since it is the first time.
+ file: ./alexandria/actmol.cpp line: 1403 check whether this is needed. likely it is here, since it is the first time.
+ file: ./alexandria/acmfitnesscomputer.cpp line: 66 fix printing
+ file: ./alexandria/acmfitnesscomputer.cpp line: 111 : implement broadcast
+ file: ./alexandria/molhandler.h line: 70 : another method without this
+ file: ./alexandria/molhandler.h line: 73 : create another method without this argument, then
+ file: ./alexandria/confighandler.h line: 104 : should this be under the ga directory???
+ file: ./alexandria/confighandler.h line: 276 : this should be general
+ file: ./alexandria/openmm_xml.cpp line: 69 : comment each element
+ file: ./alexandria/openmm_xml.cpp line: 516 : optimize values
+ file: ./alexandria/openmm_xml.cpp line: 528 : optimize values
+ file: ./alexandria/openmm_xml.cpp line: 540 : optimize values
+ file: ./alexandria/openmm_xml.cpp line: 1030 get data from topology instead of making stuff up.
+ file: ./alexandria/openmm_xml.cpp line: 1050 look up this number!
+ file: ./alexandria/openmm_xml.cpp line: 1055 check that the parameters are correct.
+ file: ./alexandria/molgen.cpp line: 415 check the loop over multiple ids
+ file: ./alexandria/molgen.cpp line: 476 check multiple ids
+ file: ./alexandria/molgen.cpp line: 897 : make sure the master has a bit less work to do
+ file: ./alexandria/molgen.cpp line: 981 : free mycomms
+ file: ./alexandria/molgen.cpp line: 1066 checks for energy should be done only when energy is a target for fitting.
+ file: ./alexandria/molgen.cpp line: 1113 check that this is correct.
+ file: ./alexandria/topology.cpp line: 142 : allow adding multiple shells
+ file: ./alexandria/topology.cpp line: 155 multiple shell support
+ file: ./alexandria/train_ff.h line: 68 : the class that does all the optimization work.
+ file: ./alexandria/secondvirial.cpp line: 557 : write out the math.
+ file: ./alexandria/secondvirial.cpp line: 664 : fix the torque contribution
+ file: ./alexandria/actmiddleman.cpp line: 65 : what about those false flags?
+ file: ./alexandria/actmiddleman.cpp line: 76 : we need to make some logfiles for the middlemen, because they apparently cannot write to the global logfile
+ file: ./alexandria/actmiddleman.cpp line: 111 : is this really necessary?
+ file: ./alexandria/allbondeds.h line: 84 implement support for extracting the median value
+ file: ./alexandria/mcmcmutator.cpp line: 134 : make optional
+ file: ./alexandria/mcmcmutator.cpp line: 150 make optional
+ file: ./alexandria/b2data.cpp line: 242 : there is factor 0.5 here
+ file: ./qgen/qgen_resp.h line: 59 make clear what unit is used.
+ file: ./qgen/qgen_acm.cpp line: 168 if fixed charges are in the atoms struct we do not need to look them up
+ file: ./qgen/qgen_acm.cpp line: 497 check this stuff: delta_eta * q0
+ file: ./qgen/qgen_acm.cpp line: 738 check this! only to be done when there are shells!
+ file: ./qgen/qgen_resp.cpp line: 456 check with symmetric charges and virtual sites
+ file: ./forces/forcecomputer.cpp line: 130 optimize this protocol using overrelaxation
+ file: ./forces/forcecomputer.cpp line: 411 take a into account
+ file: ./forces/vsitehandler.cpp line: 1100 more checking.
