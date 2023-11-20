import os, random, sys
from act_openmm import ActOpenMMSim
from get_csv_rows import *

class GeneralCouplingTheory:
    '''
    Class to facilitate force field optimization in the liquid phase.
    '''
    def __init__(self, coupleTypes:str, logfile:str, actff:str):
        self.log = None
        if not os.path.exists(coupleTypes):
            sys.exit("File %s does not exist", coupleTypes)
        self.actff = actff
        self.log = open(logfile, "w")
        self.log.write("Welcome to the return of the son of General Coupling Theory!\n")
        self.log.write("Will update ACT force field file %s\n" % self.actff)
        self.get_couple_types(coupleTypes)

    def __del__(self):
        if None != self.log:
            # Close log file
            self.log.close()
            # Close parameter files
            for obs in self.couple_types.keys():
                for atom in self.couple_types[obs].keys():
                    for param in self.couple_types[obs][atom].keys():
                        self.conv[obs][atom][param].close()

    def get_couple_types(self, coupleTypes:str):
        self.couple_types = {}
        lineno = 1
        for row in get_csv_rows(coupleTypes, 7, delim=','):
            observable = row[0]
            if not observable in self.couple_types:
                self.couple_types[observable] = {}
            atom  = row[1]
            param = row[2]
            if not atom in self.couple_types[observable]:
                self.couple_types[observable][atom] = {}
            elif atom in self.couple_types[observable] and param in self.couple_types[observable][atom]:
                sys.exit("Inconsistency in %s, atom % with param %s occurs multiple times for observable %s" %
                         ( coupleTypes, atom, param, observable ))
            try:
                self.couple_types[observable][atom][param] = { "min": float(row[3]),
                                                               "max": float(row[4]),
                                                               "value": float(row[5]),
                                                               "slope": float(row[6]) }
            except ValueError:
                sys.exit("Could not interpret line %d in %s" % ( lineno, coupleTypes ) )
            lineno += 1

        self.conv = {}
        for obs in self.couple_types.keys():
            self.conv[obs] = {}
            for atom in self.couple_types[obs].keys():
                self.conv[obs][atom] = {}
                for param in self.couple_types[obs][atom].keys():
                    pp = self.couple_types[obs][atom][param]
                    self.log.write("Will adapt %s for %s to fit %s within %g and %g current %g\n" %
                                   ( param, atom, obs, pp["min"], pp["max"], pp["value"]))
                    self.conv[obs][atom][param] = open("%s-%s-%s.xvg" % ( obs, atom, param ), "w")

    def setTargets(self, targets:list):
        self.targets = targets
        for t in self.targets:
            self.log.write("Target %s value %s\n" % ( t["observable"], t["value"] ))

    def do_iter(self, monomer_pdb:str, bulk_pdb:str, monomer_dat:str, bulk_dat:str, iter:int, pfraction:float):
        # It is much faster to convert the ACT force field file for a monomer
        sim1 = ActOpenMMSim(pdbfile=monomer_pdb, datfile=monomer_dat, actfile=self.actff)
        sim1.setup()
        emonomer = sim1.minimize()
        # Now let's do the bulk
        sim = ActOpenMMSim(pdbfile=bulk_pdb, datfile=bulk_dat, xmlfile="act.xml")
        sim.set_monomer_energy(emonomer)
        sim.run()
        # Store coordinates in new file
        coords = ("bulk%d.pdb" % iter )
        sim.write_coordinates(coords)
        # Extract target properties
        tmap = { "rho": "Density (g/mL)", 
                 "dhvap": "Potential Energy (kJ/mole)"  }
        # Prepare to extract values from the energy file
        observations = sim.log_to_average(tmap)
        # Now loop over the fitting targets
        for t in self.targets:
            myobs = t["observable"]
            myval = observations[myobs]
            if "dhvap" == myobs:
                observations[myobs] = sim.dhvap(myval)
                deviation =  observations[myobs] - t["value"]
            elif "rho" == myobs:
                deviation = myval-t["value"]
            else:
                sys.exit("Unknown observable %s" % myobs)
            # Pick parameter to change for this observable
            myobs  = t["observable"]
            # TODO move this check out of this loop
            if not myobs in self.couple_types:
                sys.exit("Do not know how to optimize %s" % myobs)
            # Randomly pick which atom to update
            atypelist = list(self.couple_types[myobs].keys())
            pindex    = random.randint(0, len(atypelist)-1)
            pkey      = atypelist[pindex]
            # Update force field
            for param in self.couple_types[myobs][pkey].keys():
                myparam = self.couple_types[myobs][pkey][param]
                pstep   = random.random()*pfraction*(myparam["max"] - myparam["min"])*deviation*myparam["slope"]
                myparam["step"] = pstep
        return coords, observations

    def reset_step(self):
        for obs in self.couple_types.keys():
            for atom in self.couple_types[obs].keys():
                for param in self.couple_types[obs][atom].keys():
                    self.couple_types[obs][atom][param]["step"] = 0

    def update_ff(self, myiter:int):
        for obs in self.couple_types.keys():
            for atom in self.couple_types[obs].keys():
                for param in self.couple_types[obs][atom].keys():
                    myparam = self.couple_types[obs][atom][param]
                    oldval = myparam["value"]
                    myparam["value"] += myparam["step"]
                    myparam["value"]  = max(myparam["min"], min(myparam["value"], myparam["max"]))
                    if oldval != myparam["value"]:
                        self.log.write("Iter %d obs %s Changing %s %s from %g to %g\n" %
                                       ( myiter, obs, atom, param, oldval, myparam["value"]))
                        self.log.flush()
                        os.system("alexandria edit_ff -ff %s -o %s -p %s -val %g -a %s" % ( self.actff, self.actff, param,
                                                                                            myparam["value"], atom ))

    def print_convergence(self, myiter:int):
        for obs in self.couple_types.keys():
            for atom in self.couple_types[obs].keys():
                for param in self.couple_types[obs][atom].keys():
                    pp = self.couple_types[obs][atom][param]
                    self.conv[obs][atom][param].write("%5d  %10g\n" % ( myiter, pp["value"] ))

    def run(self, monomer_pdb:str, bulk_pdb:str, monomer_dat:str, bulk_dat:str,
            niter:int, pfraction:float):
        # Loop over the iterations
        coords = bulk_pdb
        outf   = {}
        for t in self.targets:
            myobs = t["observable"]
            outf[myobs] = open(("%s.xvg" % myobs), "w")
        for myiter in range(niter):
            # Set all the parameter change steps to 0
            self.reset_step()
            # Do a simulation and collect data afterwards
            coords, observations = self.do_iter(monomer_pdb, coords, monomer_dat, bulk_dat, myiter, pfraction)
            # Update the force field
            self.update_ff(myiter)
            # And print stuff!
            self.print_convergence(myiter)
            for t in self.targets:
                myobs = t["observable"]
                outf[myobs].write("%5d  %10g\n" % ( myiter, observations[myobs] ))
                outf[myobs].flush()
        for t in self.targets:
            myobs = t["observable"]
            outf[myobs].close()
