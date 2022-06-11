#!/usr/bin/env python3
        
import os, shutil, sys
from enum import Enum
from get_csv_rows import *

class Alg(Enum):
    GA     = 1
    MCMC   = 2
    HYBRID = 3
    
class Target(Enum):
    EEM   = 1
    Epot  = 2
    Force = 3
    Freq  = 4

class ACT:
    '''Simple class to run Alexandria Chemistry Toolkit programs.'''
    def __init__(self,
                 MolPropFile: str,
                 SelectionFile: str, 
                 verbose: bool=False):
        try:
            self.actdata = os.environ["ACTDATA"]
        except:
            sys.exit("You have not sourced the ACT environment")
        self.verbose = verbose
        if self.verbose:
            print("Welcome to the Alexandria Chemistry Toolkit Python interface.")
        for myfiletest in [ MolPropFile, SelectionFile ]:
            if not os.path.exists(myfiletest):
                sys.exit("File %s does not exist" % myfiletest)
        self.molpropfile   = MolPropFile
        self.selectionfile = SelectionFile
        self.analyze_nodes()
        self.analyze_selection()
        # Todo implement algorithm to suggest reasonable pop_size
        self.pop_size = self.node_count
        self.set_algorithm(Alg.HYBRID, self.pop_size)
    
    def analyze_nodes(self):
        scpt = "SLURM_CPUS_PER_TASK"
        if scpt in os.environ:
            self.node_count = int(os.environ[scpt])
        else:
            self.node_count = os.cpu_count()
        if self.verbose:
            print("Found %d available cores." % self.node_count)

    def analyze_selection(self):
        self.ntrain = 0
        self.ntest  = 0
        for row in get_csv_rows(self.selectionfile, 2):
            if row[1] == "Train":
                self.ntrain += 1
            elif row[1] == "Test":
                self.ntest += 1
        if self.verbose:
            print("There are %d Train and %d Test compounds." % 
                  ( self.ntrain, self.ntest ) )

    def set_algorithm(self, algorithm:Alg, popsize:int):
        self.algopts = { "-optimizer":       algorithm.name,
                         "-cp_gen_interval": "5",
                         "-cp_pop_frac":     "0.2",
                         "-max_generations": "20",
                         "-pop_size":        str(popsize),
                         "-pr_cross":        "0.5",
                         "-n_crossovers":    "5",
                         "-temp":            "20",
                         "-tweight":         "",
                         "-max_iter":        "20",
                         "-anneal":          "0.5",
                         "-noremovemol":     "",
                         "-mindata":         "1",
                         "-debug":           "0" }

    def runpar(self) -> str:
        runit = shutil.which("srun")
        if None != runit:
            return runit
        runit = shutil.which("mpirun")
        if None != runit:
            return ( "%s -n %d -oversubscribe " % ( runit, self.node_count ))
        return ""
        
    def bastat(self, ForceFieldFileIn: str, ForceFieldFileOut: str,
               LogFile:str, options: dict):
        if not os.path.exists(ForceFieldFileIn):
            sys.exit("No force field file %s" % ForceFieldFileIn)
        cmd = ( "alexandria bastat -ff %s -o %s -mp %s -sel %s -g %s" % 
                ( ForceFieldFileIn, ForceFieldFileOut,
                  self.molpropfile, self.selectionfile, LogFile ) )
        for opt in options:
            cmd += ( " %s %s " % ( opt, options[opt] ))
        if self.verbose:
            print(cmd)
        os.system(cmd)
        
    def tune_ff(self, ForceFieldFileIn: str, ForceFieldFileOut: str,
                LogFile:str, target: Target, OptimizeGeometry: bool, options: dict):
        if not os.path.exists(ForceFieldFileIn):
            sys.exit("No force field file %s" % ForceFieldFileIn)
        cmd = ( "%s alexandria tune_ff -ff %s -o %s -mp %s -sel %s -g %s" % 
                ( self.runpar(), ForceFieldFileIn, ForceFieldFileOut,
                  self.molpropfile, self.selectionfile, LogFile ) )
        for opt in options:
            cmd += ( " %s %s " % ( opt, options[opt] ))
        ener_params = [ "sigma", "epsilon", "gamma", "kt", "klin", "kimp", "De", "D0", "beta" ]
        if OptimizeGeometry:
            ener_params.append("bondlength")
            ener_params.append("angle")
        fit_params = "'"
        for e in ener_params:
            fit_params += " " + e
        fit_params += "'"
        if target == Target.EEM:
            myopts = { "-fc_esp":          "1", 
                       "-fc_charge":       "400", 
                       "-fc_unphysical":   "100", 
                       "-fc_polar":        "1",
                       "-fit":             "'alpha zeta charge eta chi delta_eta delta_chi'" }
        elif target == Target.Epot:
            myopts = { "-fc_epot":  "1",
                       "-fit":      fit_params }
        elif target == Target.Force:
            myopts = { "-fc_epot":  "1",
                       "-fc_force": "0.1",
                       "-fit":      fit_params }
        elif target == Target.Freq:
            myopts = { "-fc_epot":  "1",
                       "-fc_force": "0.1",
                       "-fc_freq":  "0.1",
                       "-fit":      fit_params }


        for opt in myopts:
            if not opt in options:
                cmd += ( " %s %s " % ( opt, myopts[opt] ))
        for opt in self.algopts:
            if not opt in options:
                cmd += ( " %s %s " % ( opt, self.algopts[opt] ))
      
        if self.verbose:
            print(cmd)
        os.system(cmd)
     
if __name__ == '__main__':
    MolPropFile       = "molprop.xml"
    SelectionFile     = "sel.dat"
    act = ACT(MolPropFile, SelectionFile, True)
    
    ForceFieldFileIn  = "ACS-pg.xml"
    ACT.bastat(ForceFieldFileIn, ForceFieldFileIn, "bastat.log", {})
    for target in Target:
        ForceFieldFileOut = ( "tune_ff_%s.xml" % ( target.name ) )
        LogFile           = ( "tune_ff_%s.log" % ( target.name ) )
        act.tune_ff(ForceFieldFileIn, ForceFieldFileOut,
                    LogFile, target, { "-max_generations": 5 })
        ForceFieldFileIn  = ForceFieldFileOut
    
