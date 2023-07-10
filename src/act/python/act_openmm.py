#!/usr/bin/env python3

# OpenMM python example script for running Alexandria force fields
# in OpenMM using a user-selected integrator.
# This script implements a modified Buckingham potential with
# Hogervorst combination rules and Gaussian distributed charges 
# for the nonbonded interactions. 
# A Morse or Cubic potential is used for the bonded interactions.
############################################################ 
#              PROCEED AT YOUR OWN RISK.                   #
############################################################
# Author: Marie-Madeleine Walz, David van der Spoel Group,
# Department of Cell and Molecular Biology, Uppsala University, Sweden. 
# marie-madeleine.walz@icm.uu.se
############################################################

from openmm import *
from openmm.app import *
from simtk.unit import *
from simtk import openmm, unit
import numpy as np
import argparse, math, sys, shutil
from enum import Enum

ONE_4PI_EPS0 = "138.93544561"

# To distinguish 
class VdW(Enum):
    LJ8_6  = 1
    LJ12_6 = 2
    LJ14_7 = 3
    WBH    = 4
    GBHAM  = 5

# Map strings to VdW entries.    
VdWdict = { 'LJ8_6':VdW.LJ8_6, 'LJ12_6': VdW.LJ12_6, 'LJ14_7': VdW.LJ14_7, 'WBH': VdW.WBH, 'GBHAM': VdW.GBHAM }
# Make reverse map as well.
dictVdW = {}
for key in VdWdict:
    dictVdW[VdWdict[key]] = key

class SimParams:
    def __init__(self, filename:str):
        self.filename = filename
        self.params    = {}
        with open(filename, 'r') as inFileStream:
            for line in inFileStream:
                try:
                    key, equal, value =line.split() 
                    self.params[key] = value
                except:
                    continue

    def getFloat(self, key:str, default=0) -> float:
        if key in self.params and len(self.params[key]) > 0:
            try:
                if self.params[key].find("*") > 0:
                    words = self.params[key].split("*")
                elif self.params[key].find("/") > 0:
                    words = self.params[key].split("/")
                else:
                    words = [ self.params[key] ]
                value = float(words[0])
            except ValueError:
                sys.exit("Incorrect float value '%s' for key '%s' in %s" % ( words[0], key, self.filename ))
            return value
        else:
            print("Unknown or empty key '%s' in %s, using default value = %g" % ( key, self.filename, default ))
            return default

    def getInt(self, key:str) -> int:
        if key in self.params and len(self.params[key]) > 0:
            try:
                words = self.params[key].split("*")
                value = int(words[0])
            except ValueError:
                sys.exit("Incorrect integer value '%s' for key '%s' in %s" % ( words[0], key, self.filename ))
            return value
        else:
            sys.exit("Unknown or empty key '%s' in %s" % ( key, self.filename ))
        
    def getStr(self, key:str) -> str:
        if key in self.params and len(self.params[key]) > 0:
            return self.params[key]
        else:
            sys.exit("Unknown or empty key '%s' in %s" % ( key, self.filename ))
        
    def getBool(self, key:str) -> bool:
        if key in self.params and len(self.params[key]) > 0:
            return self.params[key] in [ "True", "true" ]
        else:
            sys.exit("Unknown or empty key '%s' in %s" % ( key, self.filename ))
        
class CombinationRules:
    def __init__(self, qdist:str, comb:str, args):
        self.qdist = qdist
        self.comb  = comb
        self.args = args

    def geometricString(self, vara:str, varb:str)->str:
        return ("sqrt(%s*%s)" % ( vara, varb ))
        
    def geometric(self, vara:float, varb:float)->float:
        return eval(self.geometricString(str(vara), str(varb)))
        
    def arithmeticString(self, vara:str, varb:str)->str:
        return ("0.5*(%s+%s)" % ( vara, varb ))

    def arithmetic(self, vara:float, varb:float)->float:
        return eval(self.arithmeticString(str(vara), str(varb)))

    def combStrings(self, vdw:VdW):
        if "Hogervorst" == self.comb:
            cepsilon = "((2 * epsilon1 * epsilon2)/(epsilon1 + epsilon2))"
            cgamma   = "((gamma1 + gamma2)/2)"

            if vdw == VdW.WBH:
                csigma   = "(((sqrt(((epsilon1*gamma1*(sigma1^6))/(gamma1-6)) * ((epsilon2*gamma2*(sigma2^6))/(gamma2-6)))*(gamma-6))/(epsilon*gamma))^(1.0/6.0))"
                return csigma, cepsilon, cgamma 
            
            elif vdw == VdW.GBHAM:
                crmin   = "(((sqrt(((epsilon1*gamma1*rmin1^6)/(gamma1-6)) * ((epsilon2*gamma2*rmin2^6)/(gamma2-6)))*(gamma-6))/(epsilon*gamma))^(1.0/6.0))"
                #This combination rule is likely wrong for GWBH 4 parameter one
                print("WARNING, Hogervorst combination rule and 4-parameter generalized Wang-Buckingham potential has likely no basis for their usage together...")
                cdelta = "(sqrt(delta1*delta2))"
                return crmin, cepsilon, cgamma, cdelta
        elif "Geometric" == self.comb:
            cepsilon = "(sqrt(epsilon1*epsilon2))"
            cgamma   = "(sqrt(gamma1*gamma2))"
            if vdw == VdW.WBH:
                csigma   = "(sqrt(sigma1*sigma2))"
                return csigma, cepsilon, cgamma
            elif vdw == VdW.GBHAM:
               rmin   = "(sqrt(rmin11*rmin22))"
               cdelta = "(sqrt(delta1*delta2))"
               return crmin, cepsilon, cgamma, cdelta
            elif vdw == VdW.LJ14_7:
               csigma   = "(sqrt(sigma1*sigma2))"
               cdelta = "(sqrt(delta1*delta2))"
               return csigma, cepsilon, cgamma, cdelta
            else:
                sys.exit("No support for %s" % dictVdW[vdw])
        else:
            sys.exit("Unknown combination rule '%s'" % self.comb)
            
    def zetaString(self)->str:
        if self.qdist == "Gaussian":
            return ("(zeta1*zeta2/sqrt(zeta1^2+zeta2^2))")
        else:
            sys.exit("No support for charge distribution type %s" % self.qdist)

class ActForce:
    def __init__(self, fcname:str, fgnumber:int):
        self.fcname = fcname
        self.fgnumber = fgnumber

class ActOpenMMSim:
    def __init__(self, pdbfile:str=None, xmlfile:str=None, datfile:str=None, actfile:str=None,
                 txtfile=None, pdbtraj=None, dcdtraj=None, polarizable:bool=None):
        self.txt         = None
        self.args        = self.argParser(pdbfile, xmlfile, datfile, actfile, txtfile, pdbtraj, dcdtraj, polarizable)
        self.pdb         = PDBFile(self.args.pdb_file)
        self.sim_params  = SimParams(self.args.dat_file)
        self.comb        = CombinationRules(self.sim_params.getStr("charge_distribution"),
                                            self.sim_params.getStr("combination_rule"), self.args)
        self.force_group = None
        vdwopt           = 'vanderwaals'
        vdw              = self.sim_params.getStr(vdwopt)
        if not vdw in VdWdict:
            sys.exit("Unknown value for option %s in %s" % ( vdwopt, self.args.dat_file ))
        
        self.vdw         = VdWdict[vdw]
        self.txt_header()
        self.gen_ff()
    
    def __del__(self):
        if None != self.txt:
            print("Please check output in %s" % self.args.txt_file )
            print("Energies are in %s" % self.args.ene_file )
            if None != self.args.dcd_traj:
                print("DCD trajectory is in %s" % self.args.dcd_traj )
            if None != self.args.pdb_traj:
                print("PDB trajectory in %s" % self.args.pdb_traj )
            self.txt.close()
        
    def txt_header(self):
        self.txt = open(self.args.txt_file, "w")
        self.txt.write("Starting OpenMM calculation using the ActOpenMMSim interface.\n")
        self.txt.write("input pdbfile:         %s\n" % self.args.pdb_file)
        self.txt.write("simulation parameters: %s\n" % self.args.dat_file)
        self.txt.write("vanderwaals:           %s\n" % dictVdW[self.vdw])
        self.txt.write("charge distribution:   %s\n" % self.comb.qdist)

    def gen_ff(self):
        if None != self.args.xml_file:
            self.forcefield  = ForceField(self.args.xml_file)
            self.txt.write("force field            %s\n" % self.args.xml_file)
        elif None != self.args.act_file:
            # Run alexandria gentop, but first check whether we have alexandria
            if None == shutil.which("alexandria"):
                sys.exit("You provided and ACT force field file, but the alexandria program is not in your PATH")
            self.xmloutfile = "act.xml"
            if os.path.exists(self.xmloutfile):
                if self.args.verbose:
                    print("Removing existing OpenMM force field file %s" % self.xmloutfile)
                os.unlink(self.xmloutfile)
            mycmd = ("alexandria gentop -ff %s -f %s -openmm %s" % ( self.args.act_file, 
                                                                     self.args.pdb_file,
                                                                     self.xmloutfile ))
            os.system(mycmd)
            if not os.path.exists(self.xmloutfile):
                sys.exit("Failed running '%s'" % mycmd)
            self.txt.write("Succesfully generated an OpenMM force field file %s from ACT force field %s" % (self.xmloutfile, self.args.act_file))
            self.forcefield = ForceField(self.xmloutfile)

    def xmlOutFile(self):
        return self.xmloutfile
        
    def init_force_groups(self)->int:
        self.force_group = {}
        self.fgnumber    = {}
        # First, copy existing force groups, but with only one force per group
        max_fg = -1
        self.count_forces("Init force group 1")
        for force in self.system.getForces():
            fcname   = force.getName()
            fgnumber = force.getForceGroup()
            # Find largest group while we are at it
            max_fg   = max(max_fg, fgnumber)
            if not fgnumber in self.force_group:
                self.force_group[fgnumber] = fcname
                self.fgnumber[fcname]      = fgnumber
        # Now give remaining forces new force group numbers
        self.count_forces("Init force group 2")
        for force in self.system.getForces():
            fcname   = force.getName()
            if not fcname in self.fgnumber:
                max_fg += 1
                force.setForceGroup(max_fg)
                self.force_group[max_fg] = fcname
                self.fgnumber[fcname]    = max_fg
        self.count_forces("Init force group 3")
        if self.args.verbose:
            for force in self.system.getForces():
                print("System: %s group %d" % ( force.getName(), force.getForceGroup()))
            for group in self.force_group:
                print("Self: %s group %d" % ( self.force_group[group], group ))
        return len(self.force_group.keys())

    def add_force_group(self, force, nonbond:bool, newfg:bool):
        if None == self.force_group:
            new_fgnumber = self.init_force_groups()
        else:
            new_fgnumber = len(self.force_group.keys())
        fcname   = force.getName()
        fgnumber = force.getForceGroup()
        if nonbond and self.nonbondedMethod != NoCutoff:
            # Remove direct space nonbondeds and add them to the local book-keeping
            self.txt.write("Will add direct and reciprocal part for nonbonded forces.\n")
            directname = fcname + ' (direct space)'
            self.force_group[fgnumber] = directname
            self.fgnumber[directname]  = fgnumber
            force.setName(directname)
            # Now for the PME part
            new_fgnumber += 1
            recipname = fcname + ' (reciprocal space)'
            self.force_group[new_fgnumber] = recipname
            self.fgnumber[recipname]       = new_fgnumber
            force.setReciprocalSpaceForceGroup(new_fgnumber)
        elif newfg:
            new_fgnumber += 1
            self.force_group[new_fgnumber] = fcname
            self.fgnumber[fcname]          = new_fgnumber
            force.setForceGroup(new_fgnumber)
        else:
            if fgnumber in self.force_group:
                # Update the name just in case
                self.force_group[fgnumber] = fcname
        self.count_forces("Add force group %d %s" % ( new_fgnumber, fcname))

    def del_force(self, force, nonbond:bool=False):
        fcname   = force.getName()
        fgnumber = force.getForceGroup()
        if self.args.verbose:
            print('Will try to delete force %s group %d' % (fcname, fgnumber))
        if not nonbond:
            # Find the index belonging to the force to be deleted
            # rather than using the force group number.
            iforce = -1
            index  = 0
            for force in self.system.getForces():
                if force.getName() == fcname:
                    iforce = index
                    break
                index += 1
            if -1 != iforce:
                self.system.removeForce(iforce)
            del self.fgnumber[fcname]
            del self.force_group[fgnumber]
        self.count_forces("Deleted force group %d %s" % ( fgnumber, fcname))
        return
        fnumber = []
        if fcname in self.fgnumber:
            fnumber = [ self.fgnumber[fcname] ]
        elif nonbond:
            fnumber = [ 0, 1 ]
        else:
            print("Cannot find force '%s'" % fcname)
            print(self.fgnumber)
            
        for fnum in fnumber:
            fcname  = self.force_group[fnum]
        
    def count_forces(self, label:str):
        if self.args.verbose:
            print("%s: there are %d forces"  % (label, len(self.system.getForces())))
        if self.args.debug:
            for force in self.system.getForces():
                print("DBG: fcname %s fgnumber %d" % ( force.getName(), force.getForceGroup()))

    def argParser(self, pdbfile:str, xmlfile:str, datfile:str, actfile:str, txtfile:str, pdbtraj:str, dcdtraj:str, polarizable:bool):
        mydesc = '''
        Run an OpenMM simulation based on an Alexandria force field generated by the
        Alexandria Chemistry Toolkit. Alternatively, you can reuse an existing OpenMM
        force field file generated previously. You also need to provide a pdb file with
        atomic coordinates and CONECT records to indicate the bonds and a simulation parameter file.
        
        Note that the options selected for polarizability, Van der Waals potential etc. should match
        the force field or incorrect results will be produced. The following Van der Waals potentials
        are more or less supported:
        LJ12_6 for normalLennard-Jones,
        LJ8_6 for Lennard-Jones with softer repuslion,
        LJ14_7 for 14-7 potential of Qi 2016,
        WBH for modified Buckingham (Default), 
        GBHAM for generalized 4-parameter Buckingham.
        A couple other flags are available as shown below.
        '''
        parser = argparse.ArgumentParser(description=mydesc)
        parser.add_argument("-pdb", "--pdb_file", help="Input coordinates .pdb file", default=None)
        parser.add_argument("-act", "--act_file", help="Alexandria force field file, specify either this one or an OpenMM force field file", default=None)
        parser.add_argument("-xml", "--xml_file", help="OpenMM force field .xml file", default=None)
        parser.add_argument("-dat", "--dat_file", help="simulation parameter .dat file", default=None)
        defene = "energy.csv"
        parser.add_argument("-ene", "--ene_file", help="File for storing energies, default +"+defene, default=defene)
        deflog = "output.txt"
        parser.add_argument("-txt", "--txt_file", help="File for simulation information, default +"+deflog, default=deflog)
        pdbtraj = "trajectory.pdb"
        parser.add_argument("-pdbtraj", "--pdb_traj", help="Trajectory in pdb format, can become large, default "+pdbtraj, default=pdbtraj)
        dcdtraj = "trajectory.dcd"
        parser.add_argument("-dcdtraj", "--dcd_traj", help="Trajectory in dcd format, default "+dcdtraj, default=dcdtraj)
        parser.add_argument("-pol", "--polarizable", help="Turn on support for polarization", action="store_true")
        parser.add_argument("-emonomer", "--emonomer", help="Energy of monomer will be subtracted when printing energies", type=float, default=None)
        parser.add_argument("-v", "--verbose", help="Print more stuff", action="store_true")
        parser.add_argument("-dbg", "--debug", help="Print debug stuff", action="store_true")
        # Parse the arguments
        args = parser.parse_args()
        if args.debug:
            args.verbose = True

        if None != txtfile:
            args.txt_file = txtfile
        if None != pdbtraj:
            args.pdb_traj = pdbtraj
        if None != dcdtraj:
            args.dcd_traj = dcdtraj

        if None == args.pdb_file and None == pdbfile:
            sys.exit("Please pass a pdb file")
        else:
            if None == args.pdb_file:
                args.pdb_file = pdbfile
            if not os.path.exists(args.pdb_file):
                sys.exit("Error: pdb file %s does not exist" % args.pdb_file)

        if None == args.xml_file and None == xmlfile and None == args.act_file and None == actfile:
            sys.exit("Please pass an ACT or OpenMM force field file using the appropriate flag")
        else:
            if None == args.act_file:
                # ACT force field gets precedence
                if None != actfile:
                    args.act_file = actfile
                    if not os.path.exists(args.act_file):
                        sys.exit("Error: ACT force field file %s does not exist" % args.act_file)
                else:
                    if None == args.xml_file and None != xmlfile:
                        args.xml_file = xmlfile
                    if not os.path.exists(args.xml_file):
                        sys.exit("Error: OpenMM force field file %s does not exist" % args.xml_file)

        if None == args.dat_file and None == datfile:
            sys.exit("Please pass a simulation parameter file")
        else:
            if None == args.dat_file:
                args.dat_file = datfile
            if not os.path.exists(args.dat_file):
                sys.exit("Error: paramter file %s does not exist" % args.dat_file)
                
        if not args.polarizable and None != polarizable:
            args.polarizable = polarizable
            
        return args

    def nmol(self)->int:
        return self.topology.getNumResidues()
 
    def temperature(self)->float:
        return self.temperature_c
 
    def set_monomer_energy(self, emonomer:float):
        self.args.emonomer = emonomer

    def set_params(self):
        # SET SIMULATION PARAMETERS
        ################################################
        self.dt                 = self.sim_params.getFloat('dt')
        self.equilibrationSteps = self.sim_params.getInt('equilibrationSteps')
        self.steps              = self.sim_params.getInt('steps')
        
        nbmethod                        = { 'LJPME': LJPME, 'PME': PME, 'Ewald': Ewald,
                                            'CutoffPeriodic': CutoffPeriodic, 'NoCutoff': NoCutoff }
        self.nonbondedMethod           = nbmethod[self.sim_params.getStr('nonbondedMethod')]
        self.use_switching_function    = self.sim_params.getBool('use_switching_function')
        self.switch_width              = self.sim_params.getFloat('switch_width')
        self.nonbondedCutoff           = self.sim_params.getFloat('nonbondedCutoff')
        self.use_dispersion_correction = self.sim_params.getBool('use_dispersion_correction')
        self.col_freq                  = self.sim_params.getFloat('collision_frequency', 0.1) 
        self.maxDrudeDist              = self.sim_params.getFloat('maxDrudeDistance', 0.02)
        self.useAndersenThermostat     = self.sim_params.getBool('useAndersenThermostat')
        self.temperature_c             = self.sim_params.getFloat('temperature_c')
        self.useMonteCarloBarostat     = self.sim_params.getBool('useMonteCarloBarostat')
        self.useMonteCarloAnisotropicBarostat     = self.sim_params.getBool('useMonteCarloAnisotropicBarostat')
        if self.useMonteCarloAnisotropicBarostat:
            self.scaleX             = self.sim_params.getBool('scaleX')
            self.scaleY             = self.sim_params.getBool('scaleY')
            self.scaleZ             = self.sim_params.getBool('scaleZ')
            self.pressX             = self.sim_params.getFloat('pressX')
            self.pressY             = self.sim_params.getFloat('pressY')
            self.pressZ             = self.sim_params.getFloat('pressZ')
            self.pressvec           = [self.pressX,self.pressY,self.pressZ]

        constrmethod = { 'HBonds':HBonds, 'HAngles':HAngles, 'None':None }
        self.constraints            = constrmethod[self.sim_params.getStr('constraints')]
        self.rigidWater             = self.sim_params.getBool('rigidWater')
        self.constraintTolerance    = self.sim_params.getFloat('constraintTolerance')
        
        # COMPUTING PLATFORM
        ################################################
        plform = self.sim_params.getStr('usePlatform')
        self.platform = Platform.getPlatformByName(plform)
        if 'CUDA' == plform:
            properties = {'CudaPrecision': 'single'}
            self.usePrecisionCuda = self.sim_params.getStr('usePrecisionCuda')
        else:
            if self.platform.supportsDoublePrecision():
                self.txt.write("Setting precision to double\n")
#                self.platform.setPropertyValue("Precision", "double")
        self.txt.write("Using OpenMM version %s\n" % self.platform.getOpenMMVersion())
        self.txt.write("Integration time step %g ps\n" % self.dt)

    def start_output(self):
        # OUTPUT
        ################################################
        # Do not open files unnecessarily
        save = self.sim_params.getInt('saveDcd')
        self.dcdReporter = None
        if save > 0 and self.steps >= save:
            self.dcdReporter  = DCDReporter(self.args.dcd_traj, save)
        else:
            self.args.dcd_file = None
        self.dataReporter = StateDataReporter(self.args.ene_file, self.sim_params.getInt('saveEnergy'),
                                              totalSteps=self.steps,
                                              step=self.sim_params.getBool('outStep'),
                                              time=self.sim_params.getBool('outTime'),
                                              speed=self.sim_params.getBool('outSpeed'),
                                              progress=self.sim_params.getBool('outProgress'),
                                              potentialEnergy=self.sim_params.getBool('outPotentialEnergy'),
                                              kineticEnergy=self.sim_params.getBool('outKineticEnergy'),
                                              temperature=self.sim_params.getBool('outTemperature'),
                                              volume=self.sim_params.getBool('outVolume'),
                                              density=self.sim_params.getBool('outDensity'),
                                              separator=self.sim_params.getStr('outSeparator'))
        # Do not open files unnecessarily
        save = self.sim_params.getInt('checkPoint')
        self.chkReporter = None
        if save > 0 and self.steps > save:
            self.chkReporter = CheckpointReporter('checkpnt.chk', save)
        # Do not open files unnecessarily
        save = self.sim_params.getInt('savePdb')
        self.pdbReporter = None
        if save > 0 and self.steps >= save:
            self.pdbReporter = PDBReporter(self.args.pdb_traj, save)
        else:
            self.args.pdb_traj = None
        
    def make_system(self):
        # TOPOLOGY
        ################################################
        topology  = self.pdb.topology
        positions = self.pdb.positions
        self.modeller  = Modeller(topology, positions)
        self.modeller.addExtraParticles(self.forcefield)
        self.topology  = self.modeller.topology
        self.positions = self.modeller.positions
        myDrudeMass    = self.sim_params.getFloat('drudeMass')
        myEwaldErrorTolerance = self.sim_params.getFloat('ewaldErrorTolerance')
        self.rigidWater = False
        rmcom           = True
        if self.nonbondedMethod == NoCutoff:
            rmcom = False
        if self.args.verbose:
            print("Using flexible water (if present).")
        if self.args.polarizable:
            self.system = self.forcefield.createSystem(self.topology,
                                                       nonbondedMethod=self.nonbondedMethod,
                                                       nonbondedCutoff=self.nonbondedCutoff,
                                                       removeCMMotion=rmcom,
                                                       ewaldErrorTolerance=myEwaldErrorTolerance,
                                                       rigidWater=self.rigidWater,
                                                       drudeMass=myDrudeMass*unit.amu)
            if self.args.verbose:
                print("The force field is polarizable and the drude mass is %g.\nMake sure it is consistent with your force field file." % self.args.Drude_mass)
        else:
            self.system = self.forcefield.createSystem(self.topology,
                                                       nonbondedMethod=self.nonbondedMethod,
                                                       nonbondedCutoff=self.nonbondedCutoff,
                                                       removeCMMotion=rmcom,
                                                       ewaldErrorTolerance=myEwaldErrorTolerance,
                                                       rigidWater=self.rigidWater)
            if self.args.verbose:
                print("The force field is NOT polarizable.")

        # INITIAL SETTINGS FOR FORCES
        ################################################
        for force in self.system.getForces():
            if hasattr(force, 'setCutoffDistance'):
                force.setCutoffDistance(self.nonbondedCutoff)
            if hasattr(force, 'setUseSwitchingFunction'):
                force.setUseSwitchingFunction(self.use_switching_function)
            if hasattr(force, 'setSwitchingDistance'):
                switch_distance = self.nonbondedCutoff-self.switch_width
                force.setSwitchingDistance(switch_distance)
            if hasattr(force, 'setEwaldErrorTolerance'):
                force.setEwaldErrorTolerance(myEwaldErrorTolerance)
            if hasattr(force, 'setUseDispersionCorrection'):
                force.setUseDispersionCorrection(self.use_dispersion_correction)
        self.count_forces("Initial")
    
    def find_shells_cores(self, drudeforce):
        self.cores = []
        self.shells = []
        self.core_shell = []
        self.my_core  = {}
        self.my_shell = {}
        for index in range(drudeforce.getNumParticles()):
            if self.args.debug:
                print(f"Polforce {drudeforce.getParticleParameters(index)}")
            [particle, particle1, particle2, particle3, particle4, charge, pol, aniso12, aniso34] = drudeforce.getParticleParameters(index)
            self.shells.append(particle) # particle  = shell
            self.cores.append(particle1) # particle1 = core
            self.my_core[particle] = particle1
            self.my_shell[particle1] = particle
            self.core_shell.append((particle,particle1))
        if self.args.debug:
            # Checking correct atom/shell pairing
            print(f"cores      {self.cores}")
            print(f"shells     {self.shells}")
            print(f"core_shell {self.core_shell}")
            print("########################")
                
    # CODE FOR ALEXANDRIA NONBONDED FORCES
    ################################################
    def add_direct_space_force(self): 
        """
        Create a CustomNonbondedForce to calculate the direct-space force of the Alexandria
        Van der Waals potenti - Lennard-Jones and gaussian distributed charge Coulomb - point charge Coulomb,
        placing it in specified force group.
        The LJ and point charge is necessary for both the dispersion correction and for the LJPME, and for using PME
        Create a CustomBondForce to calculate the direct space force of WBH and gaussian Coulomb for interactions 
        that are excluded (besides core-shell interactions).
        """
        cnbname       = "CustomNonbondedForce"
        dforce        = "DrudeForce"
        forces        = {}
        self.customnb = None
        drudeforce    = None
        for force in self.system.getForces():
            fname = force.getName()
            if self.args.debug:
                print("Found force %s" % fname)
            forces[fname] = force
            if cnbname == fname:
                self.customnb = forces[cnbname]
            elif dforce == fname:
                drudeforce = forces[dforce]
        self.count_forces("Direct space 1")
        # There always is a regular NonbondedForce
        self.nonbondedforce  = forces['NonbondedForce']
        self.add_force_group(self.nonbondedforce, True, False)
        if drudeforce and not self.args.polarizable:
            sys.exit("There are drudes in the system but you forgot the -pol flag")
        if self.args.verbose:
            print("***************************")
            print(f"Number of particles (incl. drudes):  {self.system.getNumParticles()}")
        self.count_forces("Direct space 2")
        if self.args.polarizable:
            self.add_force_group(drudeforce, False, False)
            self.find_shells_cores(drudeforce)
        self.count_forces("Direct space 3")

        if not self.customnb:
            return
            
        #if self.nonbondedMethod != NoCutoff:
        cutoff_distance = self.nonbondedforce.getCutoffDistance()
        switch_distance = self.nonbondedforce.getSwitchingDistance()
        self.count_forces("Direct space 4")
    
        # Electrostatics is our screened Coulomb minus the point charge based potential
        expression = 'Coulomb_gauss - Coulomb_point;'
        self.qq_expression = ( "(%s*charge1*charge2*erf(zeta*r)/r)" % ONE_4PI_EPS0 )
        expression += ( 'Coulomb_gauss = %s;' % self.qq_expression )
        if self.nonbondedMethod == NoCutoff:
            expression += 'Coulomb_point = 0;'
        else:
            expression += ( 'Coulomb_point = (%s*charge1*charge2/r);' % ONE_4PI_EPS0 )
        expression += ( "zeta = %s;" % self.comb.zetaString())
            
        self.qq_correction = openmm.CustomNonbondedForce(expression)
        if self.nonbondedMethod == NoCutoff:
            self.qq_correction.setName("Coulomb"+self.comb.qdist)
        else:
            self.qq_correction.setName("CoulombCorrection"+self.comb.qdist)
        self.qq_correction.addPerParticleParameter("charge")
        self.qq_correction.addPerParticleParameter("zeta")
        self.qq_correction.setUseSwitchingFunction(self.use_switching_function)
        # We do not use self.nonbondedMethod because the exclusion correction
        # is in real space only. However, we do have to distinguish between
        # periodic and vacuum systems.
        if self.nonbondedMethod == NoCutoff:
            self.qq_correction.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
        else:
            self.qq_correction.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            self.qq_correction.setCutoffDistance(cutoff_distance)
        # TODO: Check whether we need this, likely not.
        # self.qq_correction.setSwitchingDistance(switch_distance)
        # Don't use dispersion correction for coulomb, it does not converge, see:
        # https://github.com/openmm/openmm/issues/3162
        self.qq_correction.setUseLongRangeCorrection(False)
        self.count_forces("Direct space 5")

        self.charges = []
        if self.args.verbose:
            print("There are %d particles in the nonbondedforce" % self.nonbondedforce.getNumParticles())
        vdw = self.sim_params.getStr("vanderwaals")
        for index in range(self.nonbondedforce.getNumParticles()):
            myparams = self.customnb.getParticleParameters(index)
            if vdw == "WBH":
                [_, _, _, _, charge, zeta] = myparams
                if self.args.debug:
                    print(f"nonbonded vdw sigma, epsilon, gamma, charge, zeta {myparams}")
            elif vdw == "GWBH":
                [_, _, _, _, _, charge, zeta] = myparams
                if self.args.debug:
                    print(f"nonbonded vdw rmin, epsilon, gamma, delta, charge, zeta {myparams}")
            elif vdw == "14_7":
                [_, _, _, _, _, charge, zeta] = myparams
                if self.args.debug:
                    print(f"nonbonded vdw sigma, epsilon, gamma, delta, charge, zeta {myparams}")
            else:
                sys.exit("Not implemented what to do")
            self.charges.append(charge)
            self.qq_correction.addParticle([charge, zeta])

        if self.args.debug and vdw == "WBH":
            np = self.nonbondedforce.getNumParticles()
            for i in range(np):
                [_, sigma1, epsilon1, gamma1, _, _] = self.customnb.getParticleParameters(i)
                for j in range(i,np):
                    [_, sigma2, epsilon2, gamma2, _, _] = self.customnb.getParticleParameters(j)
                    s12 = self.comb.geometric(sigma1, sigma2)
                    e12 = self.comb.geometric(epsilon1, epsilon2)
                    g12 = self.comb.geometric(gamma1, gamma2)
                    print("i %d j %d sigma %10g epsilon %10g gamma %10g" % ( i, j, s12, e12, g12 ))
                    
        for index in range(self.nonbondedforce.getNumExceptions()):
            [iatom, jatom, _, _, _] = self.nonbondedforce.getExceptionParameters(index)
            self.qq_correction.addExclusion(iatom, jatom)
            if self.args.debug:
                print("Coulomb excl %d iatom %d jatom %d" % ( index, iatom, jatom ))
        self.count_forces("Direct space 6")
        self.add_force_group(self.qq_correction, False, True)
        self.system.addForce(self.qq_correction)
        self.count_forces("Direct space 7")

        # Van der Waals, is our custom potential minus the default LJ.
        LJ_expression = 'U_LJ = 4*epsilon_LJ*((sigma_LJ/r)^12 -(sigma_LJ/r)^6);'
        LJ_expression += ('epsilon_LJ   = %s;' % self.comb.geometricString("epsilon_LJ1", "epsilon_LJ2"))
        LJ_expression += ('sigma_LJ     = %s;' % self.comb.arithmeticString("sigma_LJ1", "sigma_LJ2"))
        LJ_expression += ('sigma_LJ_rec = %s;' % self.comb.geometricString("sigma_LJ1", "sigma_LJ2"))
        if vdw == "WBH":
            expression = 'U_WBH-U_LJ;'
            if self.nonbondedMethod == NoCutoff:
                expression += 'U_LJ = 0;'
            else:
                expression += LJ_expression
            # Note that sigma really is 1/sigma, to prevent division by zero.
            self.vdw_expression =('vdW*(((2*epsilon)/(1-(3/(gamma+3)))) * (1.0/(1.0+(sigma*r)^6)) * ((3/(gamma+3))*exp(gamma*(1-(sigma*r)))-1));')
            #self.vdw_expression = ('vdW*(((2*epsilon)/(gamma3)) * (1.0/(1.0+(sigma*r)^6)) * ((3/(gamma+3))*(gamma*(1-(sigma*r)))-1));')
            #self.vdw_expression =('vdW(((2.0*epsilon)/(1.0-(3.0/(gamma+3.0)))) * ((sigma^6)/(sigma^6+r^6))* ((3.0/(gamma+3.0))*exp(gamma*(1.0-(r/sigma)))-1.0));')
            
            expression += ( 'U_WBH = %s;' % self.vdw_expression )
            csigma, cepsilon, cgamma = self.comb.combStrings(self.vdw)
            # The statements have to be in this order! They are evaluated in the reverse order apparently.
            expression += ( 'gamma3   = (gamma/(3+gamma));')
            expression += ( 'sigma    = %s;' % csigma )
            expression += ( 'epsilon  = %s;' % cepsilon )
            expression += ( 'gamma    = %s;' % cgamma )
            expression += 'vdW = (vdW1*vdW2);'
            self.vdw_correction = openmm.CustomNonbondedForce(expression)
            if self.nonbondedMethod == NoCutoff:
                self.vdw_correction.setName("VanderWaals"+vdw)
            else:
                self.vdw_correction.setName("VanderWaalsCorrection"+vdw)
            for pp in [ "sigma", "epsilon", "gamma", "vdW", "sigma_LJ", "epsilon_LJ" ]:
                self.vdw_correction.addPerParticleParameter(pp)

            self.vdw_correction.setUseSwitchingFunction(self.use_switching_function)
            # See comment above at qq_correction.
            if self.nonbondedMethod == NoCutoff:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
            else:    
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            self.vdw_correction.setCutoffDistance(cutoff_distance)

            self.vdw_correction.setSwitchingDistance(switch_distance)
            self.vdw_correction.setUseLongRangeCorrection(self.nonbondedforce.getUseDispersionCorrection())
            for index in range(self.nonbondedforce.getNumParticles()):
                [charge_LJ, sigma_LJ, epsilon_LJ] = self.nonbondedforce.getParticleParameters(index)
                if self.nonbondedMethod == NoCutoff:
                    self.nonbondedforce.setParticleParameters(index, sigma=sigma_LJ, epsilon=0, charge=0)
#                print(self.customnb.getParticleParameters(index))
                [vdW, sigma, epsilon, gamma, charge, zeta] = self.customnb.getParticleParameters(index)
                if sigma > 0:
                    sigma = 1.0/sigma
                self.vdw_correction.addParticle([sigma, epsilon, gamma, vdW, sigma_LJ, epsilon_LJ])
                if self.args.debug:
                    print("index %d sigma %g, epsilon %g, gamma %g, vdW %g, sigma_LJ %g, epsilon_LJ %g" %  (index, sigma, epsilon, gamma, vdW, sigma_LJ._value, epsilon_LJ._value ))
            for index in range(self.nonbondedforce.getNumExceptions()):
                [iatom, jatom, chargeprod, sigma, epsilon] = self.nonbondedforce.getExceptionParameters(index)
                self.vdw_correction.addExclusion(iatom, jatom)
                if self.args.debug:
                    print("VDW excl %d iatom %d jatom %d" % ( index, iatom, jatom ))
            self.add_force_group(self.vdw_correction, False, True)
            self.system.addForce(self.vdw_correction)
#################################################
        elif vdw == "GWBH":
            expression = 'U_GWBH-U_LJ;'
            if self.nonbondedMethod == NoCutoff:
                expression += 'U_LJ = 0;'
            else:
                expression += LJ_expression

            self.vdw_expression =('vdW*(        epsilon*((delta + 2*gamma + 6)/(2*gamma)) * (1/(1+((r/rmin)^6))) * (  ((6+delta)/(delta + 2*gamma + 6)) * exp(gamma*(1-(r/rmin))) -1 ) - (epsilon/(1+(r/rmin)^delta))           );')
            expression += ( 'U_GWBH = %s;' % self.vdw_expression )
            crmin, cepsilon, cgamma, cdelta = self.comb.combStrings(self.vdw)
            expression += ( 'rmin    = %s;' % crmin )
            expression += ( 'epsilon  = %s;' % cepsilon )
            expression += ( 'gamma    = %s;' % cgamma )
            expression += ( 'delta    = %s;' % cdelta )
            expression += 'vdW = vdW1*vdW2;'
            ############TODO put vdw correction in one block except for unique parameters?
            self.vdw_correction = openmm.CustomNonbondedForce(expression)
            self.vdw_correction.setName("VanderWaalsCorrection")
            for pp in [ "rmin", "epsilon", "gamma", "delta", "vdW", "sigma_LJ", "epsilon_LJ" ]:
                self.vdw_correction.addPerParticleParameter(pp)
            self.vdw_correction.setUseSwitchingFunction(self.use_switching_function)

            if self.nonbondedMethod == NoCutoff:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
            else:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            self.vdw_correction.setCutoffDistance(cutoff_distance)

            self.vdw_correction.setSwitchingDistance(switch_distance)
            self.vdw_correction.setUseLongRangeCorrection(self.nonbondedforce.getUseDispersionCorrection())

            for index in range(self.nonbondedforce.getNumParticles()):
                [charge_LJ, sigma_LJ, epsilon_LJ] = self.nonbondedforce.getParticleParameters(index)
                [vdW, rmin, epsilon, gamma, delta, charge, zeta] = self.customnb.getParticleParameters(index)
                self.vdw_correction.addParticle([rmin, epsilon, gamma, delta, vdW, sigma_LJ, epsilon_LJ])
                if self.args.debug:
                    print("index %d rmin %g, epsilon %g, gamma %g, delta %g, vdW %g, sigma_LJ %g, epsilon_LJ %g" %  (index, rmin, epsilon, gamma, delta, vdW, sigma_LJ._value, epsilon_LJ._value ))
            for index in range(self.nonbondedforce.getNumExceptions()):
                [iatom, jatom, chargeprod, sigma, epsilon] = self.nonbondedforce.getExceptionParameters(index)
                self.vdw_correction.addExclusion(iatom, jatom)
                if self.args.debug:
                    print("excl %d iatom %d jatom %d" % ( index, iatom, jatom ))
            self.add_force_group(self.vdw_correction, False, True)
            self.system.addForce(self.vdw_correction)
            self.count_forces("Direct space 8")
#################################################
        elif vdw == "14_7":
            expression = 'U_14_7-U_LJ;'
            expression += LJ_expression

            self.vdw_expression =( 'vdW*( epsilon*( ( (1+ delta)/((r/sigma)+ delta))^7 ) * ( ( (1+ gamma)/(((r/sigma)^7) +gamma )  ) -2       ) );')
            expression += ( 'U_14_7 = %s;' % self.vdw_expression )
            csigma, cepsilon, cgamma, cdelta = self.comb.combStrings(self.vdw)
            expression += ( 'sigma    = %s;' % csigma )
            expression += ( 'epsilon  = %s;' % cepsilon )
            expression += ( 'gamma    = %s;' % cgamma )
            expression += ( 'delta    = %s;' % cdelta )
            expression += 'vdW = vdW1*vdW2;'
            ############TODO put vdw correction in one block except for unique parameters?
            self.vdw_correction = openmm.CustomNonbondedForce(expression)
            self.vdw_correction.setName("VanderWaalsCorrection")
            for pp in [ "sigma", "epsilon", "gamma", "delta", "vdW", "sigma_LJ", "epsilon_LJ" ]:
                self.vdw_correction.addPerParticleParameter(pp)
            self.vdw_correction.setUseSwitchingFunction(self.use_switching_function)

            if self.nonbondedMethod == NoCutoff:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
            else:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            self.vdw_correction.setCutoffDistance(cutoff_distance)

            self.vdw_correction.setSwitchingDistance(switch_distance)
            self.vdw_correction.setUseLongRangeCorrection(self.nonbondedforce.getUseDispersionCorrection())

            for index in range(self.nonbondedforce.getNumParticles()):
                [charge_LJ, sigma_LJ, epsilon_LJ] = self.nonbondedforce.getParticleParameters(index)
                [vdW, sigma, epsilon, gamma, delta, charge, zeta] = self.customnb.getParticleParameters(index)
                self.vdw_correction.addParticle([sigma, epsilon, gamma, delta, vdW, sigma_LJ, epsilon_LJ])
                if self.args.debug:
                    print("index %d sigma %g, epsilon %g, gamma %g, delta %g, vdW %g, sigma_LJ %g, epsilon_LJ %g" %  (index, sigma, epsilon, gamma, delta, vdW, sigma_LJ._value, epsilon_LJ._value ))
            for index in range(self.nonbondedforce.getNumExceptions()):
                [iatom, jatom, chargeprod, sigma, epsilon] = self.nonbondedforce.getExceptionParameters(index)
                self.vdw_correction.addExclusion(iatom, jatom)
                if self.args.debug:
                    print("excl %d iatom %d jatom %d" % ( index, iatom, jatom ))
            self.add_force_group(self.vdw_correction, False, True)
            self.system.addForce(self.vdw_correction)

    #################################################
    def real_exclusion(self, nexcl:int, iatom:int, jatom:int)->bool:
        if self.system.isVirtualSite(iatom) or self.system.isVirtualSite(jatom):
            return True
        if nexcl == 0:
            return False
        elif nexcl == 1:
            # If we have an exclusion between two bonded atoms
            # we have to exclude the shells as well. Therefore
            # we first look up the cores for the atom numbers
            # that are passed to this routine.
            icore = iatom
            jcore = jatom
            if self.args.polarizable:
                if iatom in self.shells:
                    icore = self.my_core[iatom]
                if jatom in self.shells:
                    jcore = self.my_core[jatom]
            return ((icore,jcore) in self.bonds or (jcore,icore) in self.bonds)
        else:
            sys.exit("Cannot handle nexcl == %d" % nexcl)
        return False

    def add_excl_correction(self):
        # Add vdW and electrostactics that have been excluded.
        # This has to be done as the number of exclusions is 3 for 
        # nonbonded interactions in OpenMM and it likely less in ACT.
        # These interactions are added using two CustomBondForce entries.
        if not self.customnb:
            return
        vdwname = "VanderWaalsExclusionCorrection"
        vdw_excl_corr = openmm.CustomBondForce(self.vdw_expression)
        vdw_excl_corr.setName(vdwname)
        if self.vdw in [ VdW.WBH, VdW.LJ14_7 ]:
            vdw_excl_corr.addPerBondParameter("sigma")
        if self.vdw == VdW.GBHAM:
            vdw_excl_corr.addPerBondParameter("rmin")
        vdw_excl_corr.addPerBondParameter("epsilon")
        vdw_excl_corr.addPerBondParameter("gamma")
        vdw_excl_corr.addPerBondParameter("vdW")
        if self.vdw in [ VdW.GBHAM, VdW.LJ14_7 ]:
            vdw_excl_corr.addPerBondParameter("delta")
        qq_excl_corr = openmm.CustomBondForce(self.qq_expression)
        qq_excl_corr.setName("CoulombExclusionCorrection")
        qq_excl_corr.addPerBondParameter("charge1")
        qq_excl_corr.addPerBondParameter("charge2")
        qq_excl_corr.addPerBondParameter("zeta")

        nexclvdw = self.sim_params.getInt("nexclvdw")
        nexclqq  = self.sim_params.getInt("nexclqq")
        if self.vdw == VdW.WBH:
            csigma, cepsilon, cgamma = self.comb.combStrings(self.vdw)
        elif self.vdw == VdW.GBHAM:
            crmin, cepsilon, cgamma, cdelta = self.comb.combStrings(self.vdw)
        elif self.vdw == VdW.LJ14_7:
            csigma, cepsilon, cgamma, cdelta = self.comb.combStrings(self.vdw)
        else:
            sys.exit("Do not know how to treat Van der Waals function '%s'" % dictVdW[self.vdw])

        if self.args.debug:
            if self.vdw == VdW.GBHAM:
                print("crmin   = %s" % crmin)
            else:
                print("csigma   = %s" % csigma)
                print("cepsilon = %s" % cepsilon)
                if not self.args.vdw in [ VdW.LJ12_6 ]:
                    print("cgamma   = %s" % cgamma)
            if self.vdw in [ VdW.GBHAM, VdW.LJ14_7 ]:
                print("cdelta   = %s" % cdelta)

        if self.vdw == VdW.WBH:

            for index in range(self.nonbondedforce.getNumExceptions()):
                # Just get the excluded atoms from the regular NB force
                [iatom, jatom, chargeprod_except, sigma_except, epsilon_except] = self.nonbondedforce.getExceptionParameters(index)
                if self.args.debug:
                    print("iatom %d jatom %d" % ( iatom, jatom ))
                # Check for shell exclusions first
                if (self.args.polarizable and 
                    ((jatom,iatom) in self.core_shell or ((iatom,jatom) in self.core_shell))):
                    continue
                # And get the parameters from the Custom NB force
                [vdW1, sigma1, epsilon1, gamma1, charge1, zeta1] = self.customnb.getParticleParameters(iatom)
                [vdW2, sigma2, epsilon2, gamma2, charge2, zeta2] = self.customnb.getParticleParameters(jatom)
                if self.args.debug:
                    print(f" custom nonbonded force i {self.customnb.getParticleParameters(iatom)}")
                    print(f" custom nonbonded force j {self.customnb.getParticleParameters(jatom)}")

                # Coulomb part
                if not self.real_exclusion(nexclqq, iatom, jatom):
                    zeta12 = ((zeta1 * zeta2)/(math.sqrt(zeta1**2 + zeta2**2)))
                    qq_excl_corr.addBond(iatom, jatom, [charge1, charge2, zeta12])
                    if self.args.debug:
                        print("Adding Coul excl corr i %d j %d q1 %g q2 %g zeta12 %g" % ( iatom, jatom, charge1, charge2, zeta12))
                
                # Van der Waals part
                if (not self.real_exclusion(nexclvdw, iatom, jatom) and
                    epsilon1 > 0 and epsilon2 > 0):
                    gamma   = eval(cgamma)
                    print("gamma_ij = %g gamma_i = %g gamma_j = %g" % (gamma, gamma1, gamma2))
                    epsilon = eval(cepsilon)
                    print("epsilon_ij = %g epsilon_i = %g epsilon_j = %g" % (epsilon, epsilon1, epsilon2))
                    sigma   = eval(csigma.replace("^", "**"))
                    print("sigma_ij = %g sigma_i = %g sigma_j = %g" % (sigma, sigma1, sigma2))

                    vdW     = vdW1*vdW2
                    if vdW != 0:
                        vdw_excl_corr.addBond(iatom, jatom, [sigma, epsilon, gamma, vdW])
                        if self.args.debug:
                            print("Adding VDW excl corr i %d j %d sigma %g epsilon %g gamma %g" % ( iatom, jatom, sigma, epsilon, gamma))
###########################
        elif self.args.vdw == "GWBH":
            for index in range(self.nonbondedforce.getNumExceptions()):
                # Just get the excluded atoms from the regular NB force
                [iatom, jatom, chargeprod_except, sigma_except, epsilon_except] = self.nonbondedforce.getExceptionParameters(index)
                # Check for shell exclusions first
                if (self.args.polarizable and
                    ((jatom,iatom) in self.core_shell or ((iatom,jatom) in self.core_shell))):
                    continue
                # And get the parameters from the Custom NB force
                [vdW1, rmin1, epsilon1, gamma1, delta1, charge1, zeta1] = self.customnb.getParticleParameters(iatom)
                [vdW2, rmin2, epsilon2, gamma2, delta1, charge2, zeta2] = self.customnb.getParticleParameters(jatom)
                if self.args.debug:
                    print(f" custom nonbonded force i {self.customnb.getParticleParameters(iatom)}")
                    print(f" custom nonbonded force j {self.customnb.getParticleParameters(jatom)}")

                # Coulomb part
                if not self.real_exclusion(nexclqq, iatom, jatom):
                    zeta = ((zeta1 * zeta2)/(np.sqrt(zeta1**2 + zeta2**2)))
                    qq_excl_corr.addBond(iatom, jatom, [charge1, charge2, zeta])
                    if self.args.debug:
                        print("Adding Coul excl i %d j %d q1 %g q2 %g zeta %g" % ( iatom, jatom, charge1, charge2, zeta))

                # Van der Waals part
                if (not self.real_exclusion(nexclvdw, iatom, jatom) and
                    epsilon1 > 0 and epsilon2 > 0):
                    gamma   = eval(cgamma)
                    epsilon = eval(cepsilon)
                    rmin   = eval(crmin.replace("^", "**"))
                    delta   = eval(cdelta)
                    vdW     = vdW1*vdW2
                    if vdW != 0:
                        vdw_excl_corr.addBond(iatom, jatom, [rmin, epsilon, gamma, delta, vdW])
                        if self.args.debug:
                            print("Adding VDW excl i %d j %d sigma %g epsilon %g gamma %g delta %g" % ( iatom, jatom, rmin, epsilon, gamma, delta))
        elif self.args.vdw == "14_7":
            for index in range(self.nonbondedforce.getNumExceptions()):
                # Just get the excluded atoms from the regular NB force
                [iatom, jatom, chargeprod_except, sigma_except, epsilon_except] = self.nonbondedforce.getExceptionParameters(index)
                # Check for shell exclusions first
                if (self.args.polarizable and
                    ((jatom,iatom) in self.core_shell or ((iatom,jatom) in self.core_shell))):
                    continue
                # And get the parameters from the Custom NB force
                [vdW1, sigma1, epsilon1, gamma1, delta1, charge1, zeta1] = self.nb_correction.getParticleParameters(iatom)
                [vdW2, sigma2, epsilon2, gamma2, delta1, charge2, zeta2] = self.nb_correction.getParticleParameters(jatom)
                if self.args.debug:
                    print(f" custom nonbonded force i {self.nb_correction.getParticleParameters(iatom)}")
                    print(f" custom nonbonded force j {self.nb_correction.getParticleParameters(jatom)}")

                # Coulomb part
                if not self.real_exclusion(nexclqq, iatom, jatom):
                    zeta = ((zeta1 * zeta2)/(np.sqrt(zeta1**2 + zeta2**2)))
                    qq_excl_corr.addBond(iatom, jatom, [charge1, charge2, zeta])
                    if self.args.debug:
                        print("Adding Coul excl i %d j %d q1 %g q2 %g zeta %g" % ( iatom, jatom, charge1, charge2, zeta))

                # Van der Waals part
                if (not self.real_exclusion(nexclvdw, iatom, jatom) and
                    epsilon1 > 0 and epsilon2 > 0):
                    gamma   = eval(cgamma)
                    epsilon = eval(cepsilon)
                    sigma   = eval(csigma)
                    delta   = eval(cdelta)
                    vdW     = vdW1*vdW2
                    if vdW != 0:
                        vdw_excl_corr.addBond(iatom, jatom, [sigma, epsilon, gamma, delta, vdW])
                        if self.args.debug:
                            print("Adding VDW excl i %d j %d sigma %g epsilon %g gamma %g delta %g" % ( iatom, jatom, sigma, epsilon, gamma, delta))
        else:
            print("Unsupported Van der Waals potential %s" % self.args.vdw)


        # Finish off. Did we add any exclusion correction?
        if 0 < qq_excl_corr.getNumBonds():        
            self.add_force_group(qq_excl_corr, False, True)
            self.system.addForce(qq_excl_corr)
            self.count_forces("Excl corr 1")
        if 0 < vdw_excl_corr.getNumBonds():
            self.add_force_group(vdw_excl_corr, False, True)
            self.system.addForce(vdw_excl_corr)
            self.count_forces("Excl corr 2")
        self.count_forces("Excl corr 3")

####################


        if self.vdw == VdW.LJ14_7:
            for index in range(self.nonbondedforce.getNumExceptions()):
                # Just get the excluded atoms from the regular NB force
                [iatom, jatom, chargeprod_except, sigma_except, epsilon_except] = self.nonbondedforce.getExceptionParameters(index)
                # Check for shell exclusions first
                if (self.args.polarizable and
                    ((jatom,iatom) in self.core_shell or ((iatom,jatom) in self.core_shell))):
                    continue
                # And get the parameters from the Custom NB force
                [vdW1, sigma1, epsilon1, gamma1, delta1, charge1, zeta1] = self.nb_correction.getParticleParameters(iatom)
                [vdW2, sigma2, epsilon2, gamma2, delta1, charge2, zeta2] = self.nb_correction.getParticleParameters(jatom)
                if self.args.debug:
                    print(f" custom nonbonded force i {self.nb_correction.getParticleParameters(iatom)}")
                    print(f" custom nonbonded force j {self.nb_correction.getParticleParameters(jatom)}")

                # Coulomb part
                if not self.real_exclusion(nexclqq, iatom, jatom):
                    zeta = ((zeta1 * zeta2)/(np.sqrt(zeta1**2 + zeta2**2)))
                    qq_excl_corr.addBond(iatom, jatom, [charge1, charge2, zeta])
                    if self.args.debug:
                        print("Adding Coul excl i %d j %d q1 %g q2 %g zeta %g" % ( iatom, jatom, charge1, charge2, zeta))

                # Van der Waals part
                if (not self.real_exclusion(nexclvdw, iatom, jatom) and
                    epsilon1 > 0 and epsilon2 > 0):
                    gamma   = eval(cgamma)
                    epsilon = eval(cepsilon)
                    sigma   = eval(csigma)
                    delta   = eval(cdelta)
                    vdW     = vdW1*vdW2
                    if vdW != 0:
                        vdw_excl_corr.addBond(iatom, jatom, [sigma, epsilon, gamma, delta, vdW])
                        if self.args.debug:
                            print("Adding VDW excl i %d j %d sigma %g epsilon %g gamma %g delta %g" % ( iatom, jatom, sigma, epsilon, gamma, delta))



####################
          
#        # Now we do not need the original CustomNonbondedForce anymore
#        self.del_force(self.nb_correction)
   
    def add_bonded_forces(self):
        self.bonds    = []
        self.cb_force = None
        for cb_force in self.system.getForces():
            if 'CustomBondForce' == cb_force.getName():
                if self.args.verbose:
                    print("Found CustomBondForce")
                cb_force.setName("AlexandriaBonds")
                self.count_forces("Add Bondeds")
                self.add_force_group(cb_force, False, False)
                for bond_index in range(cb_force.getNumBonds()):
                    # Retrieve atoms (and parameters but we just want the bonds now).
                    [iatom, jatom, params ] = cb_force.getBondParameters(bond_index)
                    self.bonds.append((iatom, jatom))
                self.cb_force = cb_force
        if self.args.debug:
            print(self.bonds)

    def make_forces(self):
        # Create a new CustomNonbondedForce to mimic the direct space 
        self.add_direct_space_force()
        self.add_bonded_forces()
        self.add_excl_correction()
        for force in self.system.getForces():
            if (force.getName() in [ "CustomAngleForce", "HarmonicAngleForce" ] and
                0 == force.getNumAngles()):
                self.del_force(force)
            elif (force.getName() in [ "RBTorsionForce", "PeriodicTorsionForce" ] and
                  0 == force.getNumTorsions()):
                self.del_force(force)
            elif (force.getName() == "CMMotionRemover" and force.getFrequency() <= 0):
                self.del_force(force)

    def print_force_settings(self):
        for force in self.system.getForces():
            print("----------------------------")
            print("%s Group: %d, PBC: %s" % ( force.getName(), 
                                              force.getForceGroup(),
                                              str(force.usesPeriodicBoundaryConditions())))
            if self.customnb and force.getName() == self.customnb.getName():
                print('"Cutoff?" {0}'.format(force.getCutoffDistance()))
                print('"SwitchingDistance?" {0}'.format(force.getSwitchingDistance ()))
                print('"CustomNonbondedMethod?" {0}'.format(force.getNonbondedMethod()))
                print('"SwitchingFunction?" {0}'.format(force.getUseSwitchingFunction()))
            elif force.getName() == self.nonbondedforce.getName():
                print('"Cutoff?" {0}'.format(force.getCutoffDistance()))
                print('"SwitchingDistance?" {0}'.format(force.getSwitchingDistance ()))
                print('"NonbondedMethod?" {0}'.format(force.getNonbondedMethod()))
                print('"SwitchingFunction?" {0}'.format(force.getUseSwitchingFunction()))
                print('"Disp. Corr.?" {0}'.format(force.getUseDispersionCorrection()))
                print('"Reciprocal Force Group?" {0}'.format(force.getReciprocalSpaceForceGroup()))
            elif force.getName() in [ "CustomBondForce", "AlexandriaBonds" ]:
                print("Number of bonds/pairs %d" % ( force.getNumBonds() ) )
                if self.args.debug:
                    for bond_index in range(force.getNumBonds()):
                        # Print atoms and parameters.
                        print(force.getBondParameters(bond_index))
            elif force.getName() in [ "CustomNonbondedForce", "DrudeForce", "CoulombCorrection", "VanderWaalsCorrection" ]:
                print("Number of particles %d" % force.getNumParticles())
            elif force.getName() in [ "CustomAngleForce", "HarmonicAngleForce" ]:
                print("Angle force %s with %d angles" % (force.getName(), force.getNumAngles()))
            
               
        print("----------------------------")

    def set_algorithms(self):
        #### ethermostat / Barostat ####
        if self.nonbondedMethod != NoCutoff:
            if self.sim_params.getBool('useMonteCarloBarostat'):
                if self.args.verbose:
                    self.txt.write("Monte Carlo Barostat will be used.\n")
                self.system.addForce(MonteCarloBarostat(self.sim_params.getFloat('pressure'),
                                                        self.temperature_c,
                                                        self.sim_params.getInt('barostatInterval')))
            elif self.sim_params.getBool('useMonteCarloAnisotropicBarostat'):
                self.system.addForce(MonteCarloAnisotropicBarostat(self.pressvec,self.temperature_c,self.scaleX,self.scaleY,self.scaleZ,self.sim_params.getInt('barostatInterval'))) 
                if self.args.verbose:
                    self.txt.write(f"Monte Carlo ANISOTROPIC Barostat will be used. The dimensions that can change are: X = {self.scaleX} Y = {self.scaleY} Z = {self.scaleZ}\n")
        if self.useAndersenThermostat:
            self.system.addForce(AndersenThermostat(self.temperature_c, self.col_freq))
            if self.args.verbose:
                self.txt.write(f"Andersen Thermostat will be used with temperature {self.temperature_c}\n")

        #### Integrator ####
        friction_c    = self.sim_params.getFloat('friction_c')
        temperature_s = self.sim_params.getFloat('temperature_s')
        integrator    = self.sim_params.getStr('integrator')
        if self.args.polarizable:
            if "DrudeLangevinIntegrator" == integrator:
                self.integrator = DrudeLangevinIntegrator(self.temperature_c, friction_c, temperature_s, 
                                                          self.sim_params.getFloat('friction_s'), self.dt)
                print(DrudeLangevinIntegrator.getTemperature(self.integrator))

            elif "DrudeNoseHooverIntegrator" == integrator:
                self.integrator = DrudeNoseHooverIntegrator(self.temperature_c, friction_c, temperature_s, 
                                                            self.sim_params.getFloat('friction_s'), self.dt)
            elif "DrudeSCFIntegrator" == integrator:
                self.integrator = DrudeSCFIntegrator(self.dt)
                self.integrator.setDrudeTemperature(temperature_s)
            else:
                sys.exit("Unknown integrator %s for polarizable system" % integrator)
            if self.useAndersenThermostat and not "DrudeSCFIntegrator" == integrator:
                print("Andersen thermostat will be turned off since %s contains a built-in thermostat." % self.integrator)
                self.useAndersenThermostat = False
            self.integrator.setMaxDrudeDistance(self.maxDrudeDist)
        else:
            nhi = "NoseHooverIntegrator"
            if nhi != integrator:
                print("Unsupported integrator %s for non-polarizable system, will use %s instead" % ( integrator, nhi ))
            self.integrator = NoseHooverIntegrator(self.temperature_c, friction_c, self.dt)

        # Print some stuff yey.
        if self.args.verbose:
            print("Core Temperature %g" % self.temperature_c)
            if self.args.polarizable:
                print("Drude Temperature %g" % self.integrator.getDrudeTemperature()._value)
            print("Step size %g" % self.integrator.getStepSize()._value)

    def compute_dipole(self)->list:
        positions = self.simulation.context.getState(getPositions=True).getPositions()
        dip = [ 0, 0, 0 ]
        enm2Debye = 48.0321
        for index in range(self.system.getNumParticles()):
            for m in range(3):
                dip[m] += positions[index][m]._value * self.charges[index] * enm2Debye
        self.txt.write("\nDipole [ %g %g %g ] total %g\n" % ( dip[0], dip[1], dip[2], 
                                                              math.sqrt(dip[0]**2+dip[1]**2+dip[2]**2)))
        return dip

    def init_simulation(self):
        #### Simulation setup ####
        self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform)
        self.simulation.context.setPositions(self.positions)

        #### Set positions of shell system to almost zero) ####
        #### the shell displacement is necessary for the LJPME to work, otherwise an error is thrown:
        #### simtk.openmm.OpenMMException: Particle coordinate is nan
        positions = self.simulation.context.getState(getPositions=True).getPositions()
        new_pos = []
        for index in range(self.system.getNumParticles()):
            if (not self.args.polarizable or not index in self.shells):
                new_pos_x = positions[index][0]
                new_pos.append((new_pos_x,positions[index][1],positions[index][2]))
            if (self.args.polarizable and index in self.shells):
                new_pos_x = positions[index][0]+0.001*nanometer
                new_pos_y = positions[index][1]+0.001*nanometer
                new_pos_z = positions[index][2]+0.001*nanometer
                new_pos.append((new_pos_x,new_pos_y,new_pos_z))

        self.simulation.context.setPositions(new_pos)
        if self.args.debug:
            print(f"number of particles (incl. drudes):  {self.system.getNumParticles()}")
            for np in new_pos:
                print("%10.5f  %10.5f  %10.5f" % ( np[0]._value, np[1]._value, np[2]._value ))
        if self.customnb:
            self.qq_correction.updateParametersInContext(self.simulation.context)
            self.vdw_correction.updateParametersInContext(self.simulation.context)
        if self.cb_force:
            # Make sure the name change trickles up in the system
            self.cb_force.updateParametersInContext(self.simulation.context)

        self.nonbondedforce.updateParametersInContext(self.simulation.context)
        if self.nonbondedMethod == NoCutoff:
            # Remove the default Non-Bonded with OpenMM
            self.txt.write("Will remove standard NonBonded forces\n")
            self.del_force(self.nonbondedforce, False)
        # TODO check whether this if statement should be flipped.
        if self.vdw != VdW.LJ12_6:
            self.del_force(self.customnb)
        
    def dhvap(self, epot:float)->float:
        nmol    = self.topology.getNumResidues()
        relener = epot/nmol - self.args.emonomer
        kB      = 1.380649e-23 * 6.02214e23 / 1000
        return kB*self.temperature_c - relener
    
    def print_energy(self, title:str):
        self.txt.write("\n%s:\n" % title)
        etot = 0.0
        self.count_forces("Print energy")
        if self.args.verbose:
            for myforce in self.system.getForces():
                self.text.write("%s\n" % myforce.getName())
        for group in self.force_group:
            eterm = self.simulation.context.getState(getEnergy=True, groups=(1 << group)).getPotentialEnergy()/unit.kilojoule_per_mole
            etot += eterm
            self.txt.write('%-40s %2d %16.4f kJ/mol\n' % (self.force_group[group], group, eterm))
        potE = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()/unit.kilojoule_per_mole
        self.txt.write('Potential energy = %.2f kJ/mol. potE-etot %.2f\n' % (potE, potE-etot))
        if None != self.args.emonomer:
            nmol = self.topology.getNumResidues()
            relener = potE/nmol - self.args.emonomer
            self.txt.write('Interaction energy for %d-mer %g\n' % ( nmol, relener ))
            self.txt.write('Delta H vap %g kJ/mol\n' % ( self.dhvap(potE) ) )
        if abs(potE-etot) > 1e-3:
            self.txt.write("sum of the above %.2f\n" % (etot))
        
    def minimize_energy(self, maxIter:int)->float:
        #### Minimize and Equilibrate ####
        self.txt.write('\nPerforming energy minimization using maxIter = %d.\n' % maxIter)
        enertol = Quantity(value=1e-8, unit=kilojoule/mole)
        self.simulation.minimizeEnergy(tolerance=enertol, maxIterations=maxIter)
        return self.simulation.context.getState(getEnergy=True).getPotentialEnergy()/unit.kilojoule_per_mole

    def equilibrate(self):
        self.txt.write('\nEquilibrating for %d steps at T = %g K.\n' % ( self.equilibrationSteps, self.temperature_c) )
        self.simulation.context.setVelocitiesToTemperature(self.temperature_c)
        self.simulation.step(self.equilibrationSteps)
    
    def production(self):
        simtime = self.sim_params.getFloat('dt')*self.sim_params.getInt('steps')
        self.txt.write('\nSimulating %g ps at %g K...\n' % (simtime, self.temperature_c ))
        if None != self.dcdReporter:
            self.simulation.reporters.append(self.dcdReporter)
        if None != self.dataReporter:
            self.simulation.reporters.append(self.dataReporter)
        if None != self.pdbReporter:
            self.simulation.reporters.append(self.pdbReporter)
        if None != self.chkReporter:
            self.simulation.reporters.append(self.chkReporter)
        self.simulation.currentStep = 0
        self.simulation.step(self.steps)

    def setup(self):
        self.set_params()
        self.start_output()
        self.make_system()
        self.make_forces()
        if self.args.verbose:
            self.print_force_settings()
        self.set_algorithms()
        self.init_simulation()
        self.print_energy("Initial energies")

    def minimize(self, maxIter:int=0)->float:
        epot = self.minimize_energy(maxIter)
        self.print_energy("After minimization")
        return epot
        
    def write_coordinates(self, outfile:str):
        with open(outfile, "w") as outf:
            vecs = self.simulation.context.getState().getPeriodicBoxVectors()
            self.topology.setPeriodicBoxVectors(vecs)
            self.pdb.writeFile(self.topology,
                               self.simulation.context.getState(getPositions=True, enforcePeriodicBox=True, getParameters=True).getPositions(),
                               outf)

    def run(self):
        self.setup()
        self.minimize(maxIter=100)
        self.equilibrate()
        self.print_energy("After equilibration")
        self.production()
        self.print_energy("After production")

    def log_to_xvg(self, xvg:str, ytargets:list):
        if None == self.args.ene_file or not os.path.exists(self.args.ene_file):
            print("Could not find any log file")
        else:
            xtarget  = "Time (ps)"
            ix = -1
            iy = []
            with open(xvg, "w") as outf:
                outf.write("@ xaxis label \"%s\"\n" % xtarget)
                with open(self.args.ene_file, "r") as inf:
                    for line in inf:
                        words = line.strip().split(";")
                        if line.find("#") >= 0:
                            for i in range(len(words)):
                                if words[i].find(xtarget) >= 0:
                                    ix = i
                                else:
                                    for j in range(len(ytargets)):
                                        if words[i].find(ytargets[j]) >= 0:
                                            iy.append(i)
                        elif ix >= 0 and len(iy) > 0:
                            try:
                                outf.write("%10g" % float(words[ix]))
                                for ii in iy:
                                    outf.write("  %10g" % (float(words[ii])))
                                outf.write("\n")
                            except ValueError:
                                print("Incomprehensible line in ene_file %s" % self.args.ene_file)
                                
    def log_to_average(self, ytargets:dict)->dict:
        if None == self.args.ene_file or not os.path.exists(self.args.ene_file):
            print("Could not find any log file")
            return []
        else:
            myaver  = {}
            for i in ytargets.keys():
                myaver[i] = 0
            naver   = 0
            xtarget = "Time (ps)"
            ix      = -1
            iy      = {}
            iy_rev  = {}
            with open(self.args.ene_file, "r") as inf:
                for line in inf:
                    words = line.strip().split(";")
                    if line.find("#") >= 0:
                        for i in range(len(words)):
                            if words[i].find(xtarget) >= 0:
                                ix = i
                            else:
                                for j in ytargets.keys():
                                    if words[i].find(ytargets[j]) >= 0:
                                        iy[j]     = i
                                        iy_rev[i] = j
                    elif ix >= 0 and len(iy.keys()) > 0:
                        try:
                            for ii in iy.keys():
                                myaver[iy_rev[iy[ii]]] += float(words[iy[ii]])
                            naver += 1
                        except ValueError:
                            print("Incomprehensible line in ene_file %s" % self.args.ene_file)
            if naver > 0:
                for i in myaver.keys():
                    myaver[i] /= naver
            return myaver
