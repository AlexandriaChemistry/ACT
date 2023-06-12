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
import argparse, sys, shutil

ONE_4PI_EPS0 = 138.935456

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

    def check_consist(self):
        if "DrudeLangevinIntegrator" == self.getStr('integrator') and self.getBool('useAndersenThermostat'):
            print("Warning, Andersen thermostat and Drude Langevin Integrators are used together. The Langevin integrator contains an in-built thermostat.")


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
        return eval(self.geomString(str(vara), str(varb)))
        
    def arithmeticString(self, vara:str, varb:str)->str:
        return ("0.5*(%s+%s)" % ( vara, varb ))

    def arithmetic(self, vara:float, varb:float)->float:
        return eval(self.geomString(str(vara), str(varb)))

    def combStrings(self):
        if "Hogervorst" == self.comb:
            cepsilon = "((2 * epsilon1 * epsilon2)/(epsilon1 + epsilon2))"
            cgamma   = "((gamma1 + gamma2)/2)"

            if self.args.vdw == "WBH":
                csigma   = "(((sqrt(((epsilon1*gamma1*sigma1^6)/(gamma1-6)) * ((epsilon2*gamma2*sigma2^6)/(gamma2-6)))*(gamma-6))/(epsilon*gamma))^(1/6))"
                return csigma, cepsilon, cgamma 
            
            elif self.args.vdw == "GWBH":
                crmin   = "(((sqrt(((epsilon1*gamma1*rmin1^6)/(gamma1-6)) * ((epsilon2*gamma2*rmin2^6)/(gamma2-6)))*(gamma-6))/(epsilon*gamma))^(1/6))"
                #This combination rule is likely wrong for GWBH 4 parameter one
                print("WARNING, Hogervorst combination rule and 4-parameter generalized Wang-Buckingham potential has likely no basis for their usage together...")
                cdelta = "(sqrt(delta1*delta2))"
                return crmin, cepsilon, cgamma, cdelta
        elif "Geometric" == self.comb:
            cepsilon = "(sqrt(epsilon1*epsilon2))"
            cgamma   = "(sqrt(gamma1*gamma2))"
            if self.args.vdw == "WBH":
                csigma   = "(sqrt(sigma1*sigma2))"
                return csigma, cepsilon, cgamma
            elif self.args.vdw == "GWBH":
               rmin   = "(sqrt(rmin11*rmin22))"
               cdelta = "(sqrt(delta1*delta2))"
               return crmin, cepsilon, cgamma, cdelta
        else:
            sys.exit("Unknown combination rule '%s'" % self.comb)
            
    def zetaString(self)->str:
        if self.qdist == "Gaussian":
            return ("(zeta1*zeta2/sqrt(zeta1^2+zeta2^2))")
        else:
            sys.exit("No support for charge distribution type %s" % self.qdist)
        
class ActOpenMMSim:
    def __init__(self, pdbfile:str=None, xmlfile:str=None, datfile:str=None, actfile:str=None, polarizable:bool=None):
        self.args        = self.argParser(pdbfile, xmlfile, datfile, actfile, polarizable)
        self.pdb         = PDBFile(self.args.pdb_file)
        self.sim_params  = SimParams(self.args.dat_file)
        self.force_group = {}
        self.fgnumber    = {}
        self.comb        = CombinationRules(self.sim_params.getStr("charge_distribution"),
                                            self.sim_params.getStr("combination_rule"), self.args)
        self.logfile     = self.args.outputdir+'/'+self.args.log_file
        if not os.path.isdir(self.args.outputdir):
            os.makedirs(self.args.outputdir, exist_ok=True)

        self.gen_ff()
        
    def gen_ff(self):
        if None != self.args.xml_file:
            self.forcefield  = ForceField(self.args.xml_file)
        elif None != self.args.act_file:
            # Run alexandria gentop, but first check whether we have alexandria
            if None == shutil.which("alexandria"):
                sys.exit("You provided and ACT force field file, but the alexandria program is not in your PATH")
            xmloutfile = self.args.outputdir+"/"+"act.xml"
            if os.path.exists(xmloutfile):
                os.unlink(xmloutfile)
            mycmd = ("alexandria gentop -ff %s -f %s -openmm %s" % ( self.args.act_file, 
                                                                     self.args.pdb_file,
                                                                     xmloutfile ))
            os.system(mycmd)
            if not os.path.exists(xmloutfile):
                sys.exit("Failed running '%s'" % mycmd)
            if self.args.verbose:
                print("Succesfully generated an OpenMM force field file")
            self.forcefield = ForceField(xmloutfile)
        
    def add_force_group(self, force, nonbond:bool=False):
        fcname  = force.getName()
        fnumber = len(self.fgnumber.keys())
        if nonbond:
            directname = fcname + ' (direct space)'
            self.force_group[fnumber] = directname
            self.fgnumber[directname] = fnumber
            force.setForceGroup(fnumber)
            fnumber = len(self.fgnumber.keys())
            recipname = fcname + ' (reciprocal space)'
            self.force_group[fnumber] = recipname
            self.fgnumber[recipname]  = fnumber
            force.setReciprocalSpaceForceGroup(fnumber)
        else:
            self.force_group[fnumber] = fcname
            self.fgnumber[fcname]     = fnumber
            force.setForceGroup(fnumber)

    def del_force(self, force):
        fcname = force.getName()
        if fcname in self.fgnumber:
            fnumber = self.fgnumber[fcname]
            self.system.removeForce(fnumber)
            del self.fgnumber[fcname]
            del self.force_group[fnumber]
            self.system.removeForce(fnumber)

    def argParser(self, pdbfile:str, xmlfile:str, datfile:str, actfile:str, polarizable:bool):
        mydesc = '''
        Run an OpenMM simulation based on an Alexandria force field generated by the
        Alexandria Chemistry Toolkit. Alternatively, you can reuse an existing OpenMM
        force field file generated previously. You also need to provide a pdb file with
        atomic coordinates and a simulation parameter file. A couple other flags are
        available as shown below.
        '''
        parser = argparse.ArgumentParser(description=mydesc)
        parser.add_argument("-pdb", "--pdb_file", help="Input coordinates .pdb file", default=None)
        parser.add_argument("-act", "--act_file", help="Alexandria force field file, specify either this one or an OpenMM force field file", default=None)
        parser.add_argument("-xml", "--xml_file", help="OpenMM force field .xml file", default=None)
        parser.add_argument("-dat", "--dat_file", help="simulation parameter .dat file", default=None)
        deflog = "log.csv"
        parser.add_argument("-log", "--log_file", help="Log file for energies, default +"+deflog, default=deflog)
        parser.add_argument("-pol", "--polarizable", help="Turn on support for polarization", action="store_true")
        defout = "output"
        parser.add_argument("-odir", "--outputdir", help="Directory to write output to, default: "+defout, type=str, default=defout)
        defvdw = "WBH"
        parser.add_argument("-vdw", "--vdw", help="Van der Waals potential function. LJ for Lennard-Jones, 14_7 for 14-7 potential of Qi 2016, WBH for modified Buckingham (Default), GWBH for generalized 4-parameter Buckingham"+defvdw, type=str, default=defvdw)
        defmDrude = "0.1"
        parser.add_argument("-mDrude", "--Drude_mass", help="The mass of drude particles. Make sure it is consistent with the mass in your forcefield"+defmDrude, type=float, default=defmDrude)
        parser.add_argument("-emonomer", "--emonomer", help="Energy of monomer will be subtracted when printing energies", type=float, default=None)
        parser.add_argument("-v", "--verbose", help="Print more stuff", action="store_true")
        parser.add_argument("-dbg", "--debug", help="Print debug stuff", action="store_true")
        # Parse the arguments
        args = parser.parse_args()
        if args.debug:
            args.verbose = True
        # Do some checking
        if args.vdw != "LJ" and args.vdw != "WBH" and args.vdw != "GWBH" and args.vdw != "14_7":
            sys.exit("Unrecognized form of VDW potential. Use either LJ or WBH or GWBH keywords")

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

    def set_monomer_energy(self, emonomer:float):
        self.args.emonomer = emonomer

    def set_params(self):
        # SET SIMULATION PARAMETERS
        ################################################
        self.dt                 = self.sim_params.getFloat('dt')
        self.equilibrationSteps = self.sim_params.getInt('equilibrationSteps')
        self.steps              = self.sim_params.getInt('steps')
        self.save               = self.sim_params.getInt('save')
        
        nbmethod                        = { 'LJPME':LJPME, 'PME':PME, 'Ewald':Ewald, 
                                            'CutoffPeriodic':CutoffPeriodic, 'NoCutoff':NoCutoff }
        self.nonbondedMethod           = nbmethod[self.sim_params.getStr('nonbondedMethod')]
        self.use_switching_function    = self.sim_params.getBool('use_switching_function')
        self.switch_width              = self.sim_params.getFloat('switch_width')
        self.nonbondedCutoff           = self.sim_params.getFloat('nonbondedCutoff')
        self.use_dispersion_correction = self.sim_params.getBool('use_dispersion_correction')
        self.col_freq                  = self.sim_params.getFloat('collision_frequency', 0.1) 
        self.MaxDrudeDist              = self.sim_params.getFloat('MaxDrudeDistance', 0.02)
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

    def start_output(self):
        # OUTPUT
        ################################################
        save              = self.sim_params.getInt('save')
        self.dcdReporter  = DCDReporter(self.args.outputdir+'/trajectory.dcd', save)
        self.dataReporter = StateDataReporter(self.logfile, save, totalSteps=self.steps,
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
        self.chkReporter = CheckpointReporter(self.args.outputdir+'/checkpnt.chk', save)
        self.pdbReporter = PDBReporter(self.args.outputdir+'/output.pdb', save)
        
    def make_system(self):
        # TOPOLOGY
        ################################################
        topology  = self.pdb.topology
        positions = self.pdb.positions
        self.modeller  = Modeller(topology, positions)
        self.modeller.addExtraParticles(self.forcefield)
        self.topology  = self.modeller.topology
        self.positions = self.modeller.positions
        myDrudeMass    = self.args.Drude_mass
        myEwaldErrorTolerance = self.sim_params.getFloat('ewaldErrorTolerance')
        self.rigidWater = False
        if self.args.verbose:
            print("Using flexible water (if present).")
        if self.args.polarizable:
            self.system = self.forcefield.createSystem(self.topology,
                                                       nonbondedMethod=self.nonbondedMethod,
                                                       nonbondedCutoff=self.nonbondedCutoff,
                                                       ewaldErrorTolerance=myEwaldErrorTolerance,
                                                       rigidWater=self.rigidWater,
                                                       drudeMass=myDrudeMass*unit.amu)
            if self.args.verbose:
                print("The force field is polarizable and the drude mass is %g.\nMake sure it is consistent with your force field file." % self.args.Drude_mass)
        else:
            self.system = self.forcefield.createSystem(self.topology,
                                                       nonbondedMethod=self.nonbondedMethod,
                                                       nonbondedCutoff=self.nonbondedCutoff,
                                                       ewaldErrorTolerance=myEwaldErrorTolerance,
                                                       rigidWater=self.rigidWater)
            if self.args.verbose:
                print("The force field is NOT polarizable.")

        # SETTINGS FOR FORCES
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

    # CODE FOR ALEXANDRIA NONBONDED FORCES
    ################################################
    def add_direct_space_force(self): 
        """
        Create a CustomNonbondedForce to calculate the direct-space force of WBK-LJ and gaussian Coulomb-point Coulomb,
        placing it in specified force group.
        The LJ and point charge is necessary for both the dispersion correction and for the LJPME, and for using PME
        Create a CustomBondForce to calculate the direct space force of WBK and gaussian Coulomb for interactions 
        that are excluded (besides core-shell interactions).
        """
        cnbname      = "CustomNonbondedForce"
        haveCustomNB = False
        forces       = {}
        for force in self.system.getForces():
            fname = force.getName()
            if self.args.debug:
                print("Found force %s" % fname)
            forces[fname] = force
            if cnbname == fname:
                haveCustomNB = True
        # There always is a regular NonbondedForce
        self.nb_openmm  = forces['NonbondedForce']
        self.add_force_group(self.nb_openmm, True)
        if haveCustomNB:
            self.nb_correction = forces[cnbname]
            self.nb_correction.setName("NB_Correction")
        else:
            self.nb_correction = None
            if self.args.debug:
                print("No %s in the system" % cnbname)
        dforce = 'DrudeForce'
        if dforce in forces and not self.args.polarizable:
            sys.exit("There are drudes in the system but you forgot the -pol flag")
        if self.args.verbose:
            print("***************************")
            print(f"Number of particles (incl. drudes):  {self.system.getNumParticles()}")
        if self.args.polarizable:
            pol_force = forces[dforce]
            self.add_force_group(pol_force)
            self.cores = []
            self.shells = []
            self.core_shell = []
            for index in range(pol_force.getNumParticles()):
                if self.args.debug:
                    print(f"Polforce {pol_force.getParticleParameters(index)}")
                [particle, particle1, particle2, particle3, particle4, charge, pol, aniso12, aniso34] = pol_force.getParticleParameters(index)
                self.shells.append(particle) # particle  = shell
                self.cores.append(particle1) # particle1 = core
                self.core_shell.append((particle,particle1))
            if self.args.debug:
                # Checking correct atom/shell pairing
                print(f"cores      {self.cores}")
                print(f"shells     {self.shells}")
                print(f"core_shell {self.core_shell}")
                print("########################")


        # TODO: Only for pbc simulations
        cutoff_distance = self.nb_openmm.getCutoffDistance()
        switch_distance = self.nb_openmm.getSwitchingDistance()
    
        if not self.nb_correction:
            return
            
        # Electrostatics is our screened Coulomb minus the point charge based potential
        expression = 'Coulomb_gauss - Coulomb_point;'
        self.qq_expression = ( "(%f*charge1*charge2*erf(zeta*r)/r)" % ONE_4PI_EPS0 )
        expression += ( 'Coulomb_gauss = %s;' % self.qq_expression )
        expression += ( 'Coulomb_point = (%f*charge1*charge2/r);' % ONE_4PI_EPS0 )
        expression += ( "zeta = %s;" % self.comb.zetaString())
            
        self.qq_correction = openmm.CustomNonbondedForce(expression)
        self.qq_correction.setName("CoulombCorrection")
        self.qq_correction.addPerParticleParameter("charge")
        self.qq_correction.addPerParticleParameter("zeta")
        self.qq_correction.setUseSwitchingFunction(self.use_switching_function)
        if self.nonbondedMethod == NoCutoff:
            self.qq_correction.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
        else:    
            # TODO: why not use self.nonbondedMethod
            self.qq_correction.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        self.qq_correction.setCutoffDistance(cutoff_distance)
        self.qq_correction.setSwitchingDistance(switch_distance)
        # Don't use dispersion correction for coulomb, it does not converge, see:
        # https://github.com/openmm/openmm/issues/3162
        self.qq_correction.setUseLongRangeCorrection(False)
        
        self.charges = []
        for index in range(self.nb_openmm.getNumParticles()):
            if self.args.vdw == "WBH":
                [vdW, sigma, epsilon, gamma, charge, zeta] = self.nb_correction.getParticleParameters(index)
                if self.args.debug:
                    print(f"nonbonded vdw sigma, epsilon, gamma, charge, zeta {self.nb_correction.getParticleParameters(index)}")
            elif self.args.vdw == "GWBH":
                [vdW, rmin, epsilon, gamma, delta, charge, zeta] = self.nb_correction.getParticleParameters(index)
                if self.args.debug:
                    print(f"nonbonded vdw rmin, epsilon, gamma, delta, charge, zeta {self.nb_correction.getParticleParameters(index)}")
            elif self.args.vdw == "14_7":
                [vdW, sigma, epsilon, gamma, delta, charge, zeta] = self.nb_correction.getParticleParameters(index)
                if self.args.debug:
                    print(f"nonbonded vdw sigma, epsilon, gamma, delta, charge, zeta {self.nb_correction.getParticleParameters(index)}")
            else:
                sys.exit("Not implemented what to do")
            self.charges.append(charge)
            self.qq_correction.addParticle([charge,zeta])
            
        for index in range(self.nb_openmm.getNumExceptions()):
            [iatom, jatom, chargeprod, sigma, epsilon] = self.nb_openmm.getExceptionParameters(index)
            self.qq_correction.addExclusion(iatom, jatom)
            if self.args.debug:
                print("Coulomb excl %d iatom %d jatom %d" % ( index, iatom, jatom ))
        self.add_force_group(self.qq_correction)
        self.system.addForce(self.qq_correction)

        print("Using %s for Van der Waals interactions" % self.args.vdw)
        
        # Van der Waals, is our custom potential minus the default LJ.
        LJ_expression = 'U_LJ = 4*epsilon_LJ*((sigma_LJ/r)^12 -(sigma_LJ/r)^6);'
        LJ_expression += ('epsilon_LJ   = %s;' % self.comb.geometricString("epsilon_LJ1", "epsilon_LJ2"))
        LJ_expression += ('sigma_LJ     = %s;' % self.comb.arithmeticString("sigma_LJ1", "sigma_LJ2"))
        LJ_expression += ('sigma_LJ_rec = %s;' % self.comb.geometricString("sigma_LJ1", "sigma_LJ2"))
        if self.args.vdw == "WBH":
            expression = 'U_WKB-U_LJ;'
            expression += LJ_expression
            # TODO: Remove hard-coded combination rules if possible
            self.vdw_expression =('vdW*(((((2*epsilon)/(1-(3/(gamma+3)))) * ((sigma^6)/(sigma^6+r^6))* (((3/(gamma+3))*(exp(gamma*(1-(r/sigma)))))-1))));')
            
            expression += ( 'U_WKB = %s;' % self.vdw_expression )
            csigma, cepsilon, cgamma = self.comb.combStrings()
            # The statements have to be in this order! They are evaluated in the reverse order apparently.
            expression += ( 'sigma    = %s;' % csigma )
            expression += ( 'epsilon  = %s;' % cepsilon )
            expression += ( 'gamma    = %s;' % cgamma )
            expression += 'vdW = vdW1*vdW2;'
            self.vdw_correction = openmm.CustomNonbondedForce(expression)
            self.vdw_correction.setName("VanderWaalsCorrection")
            self.vdw_correction.addPerParticleParameter("sigma")
            self.vdw_correction.addPerParticleParameter("epsilon")
            self.vdw_correction.addPerParticleParameter("gamma")
            self.vdw_correction.addPerParticleParameter("vdW")
            self.vdw_correction.addPerParticleParameter("sigma_LJ")
            self.vdw_correction.addPerParticleParameter("epsilon_LJ")
            self.vdw_correction.setUseSwitchingFunction(self.use_switching_function)
            if self.nonbondedMethod == NoCutoff:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
            else:    
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            self.vdw_correction.setCutoffDistance(cutoff_distance)

            self.vdw_correction.setSwitchingDistance(switch_distance)
            self.vdw_correction.setUseLongRangeCorrection(self.nb_openmm.getUseDispersionCorrection())
            for index in range(self.nb_openmm.getNumParticles()):
                [charge_LJ, sigma_LJ, epsilon_LJ] = self.nb_openmm.getParticleParameters(index)
#                print(self.nb_correction.getParticleParameters(index))
                [vdW, sigma, epsilon, gamma, charge, zeta] = self.nb_correction.getParticleParameters(index)
                self.vdw_correction.addParticle([sigma, epsilon, gamma, vdW, sigma_LJ, epsilon_LJ])
                if self.args.debug:
                    print("index %d sigma %g, epsilon %g, gamma %g, vdW %g, sigma_LJ %g, epsilon_LJ %g" %  (index, sigma, epsilon, gamma, vdW, sigma_LJ._value, epsilon_LJ._value ))
            for index in range(self.nb_openmm.getNumExceptions()):
                [iatom, jatom, chargeprod, sigma, epsilon] = self.nb_openmm.getExceptionParameters(index)
                self.vdw_correction.addExclusion(iatom, jatom)
                if self.args.debug:
                    print("VDW excl %d iatom %d jatom %d" % ( index, iatom, jatom ))
            self.add_force_group(self.vdw_correction)
            self.system.addForce(self.vdw_correction)
#################################################
        elif self.args.vdw == "GWBH":
            expression = 'U_GWKB-U_LJ;'
            expression += LJ_expression

#            self.vdw_expression =('vdW*(((((2*epsilon)/(1-(3/(gamma+3)))) * ((sigma^6)/(sigma^6+r^6))* (((3/(gamma+3))*(exp(gamma*(1-(r/sigma)))))-1))));')
            self.vdw_expression =('vdW*(        epsilon*((delta + 2*gamma + 6)/(2*gamma)) * (1/(1+((r/rmin)^6))) * (  ((6+delta)/(delta + 2*gamma + 6)) * exp(gamma*(1-(r/rmin))) -1 ) - (epsilon/(1+(r/rmin)^delta))           );')
            expression += ( 'U_GWKB = %s;' % self.vdw_expression )
            crmin, cepsilon, cgamma, cdelta = self.comb.combStrings()
            expression += ( 'rmin    = %s;' % crmin )
            expression += ( 'epsilon  = %s;' % cepsilon )
            expression += ( 'gamma    = %s;' % cgamma )
            expression += ( 'delta    = %s;' % cdelta )
            expression += 'vdW = vdW1*vdW2;'
            ############TODO put vdw correction in one block except for unique parameters?
            self.vdw_correction = openmm.CustomNonbondedForce(expression)
            self.vdw_correction.setName("VanderWaalsCorrection")
            self.vdw_correction.addPerParticleParameter("rmin")
            self.vdw_correction.addPerParticleParameter("epsilon")
            self.vdw_correction.addPerParticleParameter("gamma")
            self.vdw_correction.addPerParticleParameter("delta")
            self.vdw_correction.addPerParticleParameter("vdW")
            self.vdw_correction.addPerParticleParameter("sigma_LJ")
            self.vdw_correction.addPerParticleParameter("epsilon_LJ")
            self.vdw_correction.setUseSwitchingFunction(self.use_switching_function)

            if self.nonbondedMethod == NoCutoff:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
            else:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            self.vdw_correction.setCutoffDistance(cutoff_distance)

            self.vdw_correction.setSwitchingDistance(switch_distance)
            self.vdw_correction.setUseLongRangeCorrection(self.nb_openmm.getUseDispersionCorrection())

            for index in range(self.nb_openmm.getNumParticles()):
                [charge_LJ, sigma_LJ, epsilon_LJ] = self.nb_openmm.getParticleParameters(index)
                [vdW, rmin, epsilon, gamma, delta, charge, zeta] = self.nb_correction.getParticleParameters(index)
                self.vdw_correction.addParticle([rmin, epsilon, gamma, delta, vdW, sigma_LJ, epsilon_LJ])
                if self.args.debug:
                    print("index %d rmin %g, epsilon %g, gamma %g, delta %g, vdW %g, sigma_LJ %g, epsilon_LJ %g" %  (index, rmin, epsilon, gamma, delta, vdW, sigma_LJ._value, epsilon_LJ._value ))
            for index in range(self.nb_openmm.getNumExceptions()):
                [iatom, jatom, chargeprod, sigma, epsilon] = self.nb_openmm.getExceptionParameters(index)
                self.vdw_correction.addExclusion(iatom, jatom)
                if self.args.debug:
                    print("excl %d iatom %d jatom %d" % ( index, iatom, jatom ))
            self.add_force_group(self.vdw_correction)
            self.system.addForce(self.vdw_correction)


#################################################



        elif self.args.vdw == "14_7":
            expression = 'U_14_7-U_LJ;'
            expression += LJ_expression

#            self.vdw_expression =('vdW*(((((2*epsilon)/(1-(3/(gamma+3)))) * ((sigma^6)/(sigma^6+r^6))* (((3/(gamma+3))*(exp(gamma*(1-(r/sigma)))))-1))));')
#            self.vdw_expression =('vdW*(        epsilon*((delta + 2*gamma + 6)/(2*gamma)) * (1/(1+((r/rmin)^6))) * (  ((6+delta)/(delta + 2*gamma + 6)) * exp(gamma*(1-(r/rmin))) -1 ) - (epsilon/(1+(r/rmin)^delta))           );')
            self.vdw_expression =( 'vdW*(        epsilon*( ( (1+ delta)/((r/sigma)+ delta))^7 ) * ( ( (1+ gamma)/(((r/sigma)^7) +gamma )  ) -2       ) );'                                            )
            expression += ( 'U_14_7 = %s;' % self.vdw_expression )
            csigma, cepsilon, cgamma, cdelta = self.comb.combStrings()
            expression += ( 'sigma    = %s;' % csigma )
            expression += ( 'epsilon  = %s;' % cepsilon )
            expression += ( 'gamma    = %s;' % cgamma )
            expression += ( 'delta    = %s;' % cdelta )
            expression += 'vdW = vdW1*vdW2;'
            ############TODO put vdw correction in one block except for unique parameters?
            self.vdw_correction = openmm.CustomNonbondedForce(expression)
            self.vdw_correction.setName("VanderWaalsCorrection")
            self.vdw_correction.addPerParticleParameter("sigma")
            self.vdw_correction.addPerParticleParameter("epsilon")
            self.vdw_correction.addPerParticleParameter("gamma")
            self.vdw_correction.addPerParticleParameter("delta")
            self.vdw_correction.addPerParticleParameter("vdW")
            self.vdw_correction.addPerParticleParameter("sigma_LJ")
            self.vdw_correction.addPerParticleParameter("epsilon_LJ")
            self.vdw_correction.setUseSwitchingFunction(self.use_switching_function)

            if self.nonbondedMethod == NoCutoff:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
            else:
                self.vdw_correction.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
            self.vdw_correction.setCutoffDistance(cutoff_distance)

            self.vdw_correction.setSwitchingDistance(switch_distance)
            self.vdw_correction.setUseLongRangeCorrection(self.nb_openmm.getUseDispersionCorrection())

            for index in range(self.nb_openmm.getNumParticles()):
                [charge_LJ, sigma_LJ, epsilon_LJ] = self.nb_openmm.getParticleParameters(index)
                [vdW, sigma, epsilon, gamma, delta, charge, zeta] = self.nb_correction.getParticleParameters(index)
                self.vdw_correction.addParticle([sigma, epsilon, gamma, delta, vdW, sigma_LJ, epsilon_LJ])
                if self.args.debug:
                    print("index %d sigma %g, epsilon %g, gamma %g, delta %g, vdW %g, sigma_LJ %g, epsilon_LJ %g" %  (index, sigma, epsilon, gamma, delta, vdW, sigma_LJ._value, epsilon_LJ._value ))
            for index in range(self.nb_openmm.getNumExceptions()):
                [iatom, jatom, chargeprod, sigma, epsilon] = self.nb_openmm.getExceptionParameters(index)
                self.vdw_correction.addExclusion(iatom, jatom)
                if self.args.debug:
                    print("excl %d iatom %d jatom %d" % ( index, iatom, jatom ))
            self.add_force_group(self.vdw_correction)
            self.system.addForce(self.vdw_correction)









#################################################
    def real_exclusion(self, nexcl:int, iatom:int, jatom:int)->bool:
        if nexcl == 0:
            return False
        elif nexcl == 1:
            return ((iatom,jatom) in self.bonds or (jatom,iatom) in self.bonds)
        else:
            sys.exit("Cannot handle nexcl == %d" % nexcl)
        return False

    def add_excl_correction(self):
        # Add vdW and electrostactics that have been excluded.
        # This has to be done as the number of exclusions is 3 for 
        # nonbonded interactions in OpenMM and it likely less in ACT.
        # These interactions are added using two CustomBondForce entries.
        if not self.nb_correction:
            return
        vdw_excl_corr = openmm.CustomBondForce(self.vdw_expression)
        vdw_excl_corr.setName("VanderWaalsExclusionCorrection")
#        vdw_excl_corr.addPerBondParameter("sigma")
        if self.args.vdw == "WBH" or self.args.vdw == "14_7":
            vdw_excl_corr.addPerBondParameter("sigma")
        if self.args.vdw == "GWBH":
            vdw_excl_corr.addPerBondParameter("rmin")
        vdw_excl_corr.addPerBondParameter("epsilon")
        vdw_excl_corr.addPerBondParameter("gamma")
        vdw_excl_corr.addPerBondParameter("vdW")
        if self.args.vdw == "GWBH" or self.args.vdw == "14_7":
            vdw_excl_corr.addPerBondParameter("delta")
        qq_excl_corr = openmm.CustomBondForce(self.qq_expression)
        qq_excl_corr.setName("CoulombExclusionCorrection")
        qq_excl_corr.addPerBondParameter("charge1")
        qq_excl_corr.addPerBondParameter("charge2")
        qq_excl_corr.addPerBondParameter("zeta")

        nexclvdw = self.sim_params.getInt("nexclvdw")
        nexclqq  = self.sim_params.getInt("nexclqq")
        if self.args.vdw == "WBH":

            csigma, cepsilon, cgamma = self.comb.combStrings()
        if self.args.vdw == "GWBH":
            crmin, cepsilon, cgamma, cdelta = self.comb.combStrings()
        if self.args.vdw == "14_7":
            csigma, cepsilon, cgamma, cdelta = self.comb.combStrings()

        if self.args.debug:
            if self.args.vdw == "GWBH":
                print("crmin   = %s" % crmin)
                print("cdelta   = %s" % cdelta)
            else:
                print("csigma   = %s" % csigma)
            print("cepsilon = %s" % cepsilon)
            print("cgamma   = %s" % cgamma)
            if self.args.vdw == "14_7":
                print("cgamma   = %s" % cdelta)
        if self.args.vdw == "WBH":

            for index in range(self.nb_openmm.getNumExceptions()):
                # Just get the excluded atoms from the regular NB force
                [iatom, jatom, chargeprod_except, sigma_except, epsilon_except] = self.nb_openmm.getExceptionParameters(index)
                # Check for shell exclusions first
                if (self.args.polarizable and 
                    ((jatom,iatom) in self.core_shell or ((iatom,jatom) in self.core_shell))):
                    continue
                # And get the parameters from the Custom NB force
                [vdW1, sigma1, epsilon1, gamma1, charge1, zeta1] = self.nb_correction.getParticleParameters(iatom)
                [vdW2, sigma2, epsilon2, gamma2, charge2, zeta2] = self.nb_correction.getParticleParameters(jatom)
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
                    sigma   = eval(csigma.replace("^", "**"))
                    vdW     = vdW1*vdW2
                    if vdW != 0:
                        vdw_excl_corr.addBond(iatom, jatom, [sigma, epsilon, gamma, vdW])
                        if self.args.debug:
                            print("Adding VDW excl i %d j %d sigma %g epsilon %g gamma %g" % ( iatom, jatom, sigma, epsilon, gamma))
###########################
        if self.args.vdw == "GWBH":
            for index in range(self.nb_openmm.getNumExceptions()):
                # Just get the excluded atoms from the regular NB force
                [iatom, jatom, chargeprod_except, sigma_except, epsilon_except] = self.nb_openmm.getExceptionParameters(index)
                # Check for shell exclusions first
                if (self.args.polarizable and
                    ((jatom,iatom) in self.core_shell or ((iatom,jatom) in self.core_shell))):
                    continue
                # And get the parameters from the Custom NB force
                [vdW1, rmin1, epsilon1, gamma1, delta1, charge1, zeta1] = self.nb_correction.getParticleParameters(iatom)
                [vdW2, rmin2, epsilon2, gamma2, delta1, charge2, zeta2] = self.nb_correction.getParticleParameters(jatom)
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
                    rmin   = eval(crmin.replace("^", "**"))
                    delta   = eval(cdelta)
                    vdW     = vdW1*vdW2
                    if vdW != 0:
                        vdw_excl_corr.addBond(iatom, jatom, [rmin, epsilon, gamma, delta, vdW])
                        if self.args.debug:
                            print("Adding VDW excl i %d j %d sigma %g epsilon %g gamma %g delta %g" % ( iatom, jatom, rmin, epsilon, gamma, delta))


####################


        if self.args.vdw == "14_7":
            for index in range(self.nb_openmm.getNumExceptions()):
                # Just get the excluded atoms from the regular NB force
                [iatom, jatom, chargeprod_except, sigma_except, epsilon_except] = self.nb_openmm.getExceptionParameters(index)
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
                        
        self.add_force_group(qq_excl_corr)
        self.system.addForce(qq_excl_corr)
        self.add_force_group(vdw_excl_corr)
        self.system.addForce(vdw_excl_corr)
        # Now we do not need the original CustomNonbondedForce anymore
        self.del_force(self.nb_correction)
   
    def add_bonded_forces(self):
        self.bonds = []
        for cb_force in self.system.getForces():
            if 'CustomBondForce' == cb_force.getName():
                if self.args.verbose:
                    print("Found CustomBondForce")
                cb_force.setName("AlexandriaBonds")
                self.add_force_group(cb_force)
                for bond_index in range(cb_force.getNumBonds()):
                    # Retrieve atoms (and parameters but we just want the bonds now).
                    [iatom, jatom, params ] = cb_force.getBondParameters(bond_index)
                    self.bonds.append((iatom, jatom))
        if self.args.debug:
            print(self.bonds)

    def make_forces(self):
        # Create a new CustomNonbondedForce to mimic the direct space 
        self.add_direct_space_force()
        self.add_bonded_forces()
        self.add_excl_correction()

    def print_force_settings(self):
        for force in self.system.getForces():
            print("----------------------------")
            print("%s Group: '%s' %d, PBC: %s" % ( force.getName(), 
                                                   self.force_group[force.getForceGroup()],
                                                   force.getForceGroup(),
                                                   str(force.usesPeriodicBoundaryConditions())))
            if self.nb_correction and force.getName() == self.nb_correction.getName():
                print('"Cutoff?" {0}'.format(force.getCutoffDistance()))
                print('"SwitchingDistance?" {0}'.format(force.getSwitchingDistance ()))
                print('"CustomNonbondedMethod?" {0}'.format(force.getNonbondedMethod()))
                print('"SwitchingFunction?" {0}'.format(force.getUseSwitchingFunction()))
            elif force.getName() == self.nb_openmm.getName():
                print('"Cutoff?" {0}'.format(force.getCutoffDistance()))
                print('"SwitchingDistance?" {0}'.format(force.getSwitchingDistance ()))
                print('"NonbondedMethod?" {0}'.format(force.getNonbondedMethod()))
                print('"SwitchingFunction?" {0}'.format(force.getUseSwitchingFunction()))
                print('"Disp. Corr.?" {0}'.format(force.getUseDispersionCorrection()))
                print('"Reciprocal Force Group?" {0}'.format(force.getReciprocalSpaceForceGroup()))
            elif force.getName() in [ "CustomBondForce", "AlexandriaBonds", "QQ", "VDW" ]:
                print("Number of bonds/pairs %d" % ( force.getNumBonds() ) )
                if self.args.debug:
                    for bond_index in range(force.getNumBonds()):
                        # Print atoms and parameters.
                        print(force.getBondParameters(bond_index))
            elif force.getName() in [ "CustomNonbondedForce", "DrudeForce", "CoulombCorrection", "LennardJonesCorrection" ]:
                print("Number of particles %d" % force.getNumParticles())
            
               
        print("----------------------------")

    def set_algorithms(self):
        #### ethermostat / Barostat ####
        if self.nonbondedMethod != NoCutoff:
            if self.sim_params.getBool('useMonteCarloBarostat'):
                if self.args.verbose:
                    print("Monte Carlo Barostat will be used.")
                self.system.addForce(MonteCarloBarostat(self.sim_params.getFloat('pressure'),
                                                        self.temperature_c,
                                                        self.sim_params.getInt('barostatInterval')))
            elif self.sim_params.getBool('useMonteCarloAnisotropicBarostat'):
                self.system.addForce(MonteCarloAnisotropicBarostat(self.pressvec,self.temperature_c,self.scaleX,self.scaleY,self.scaleZ,self.sim_params.getInt('barostatInterval'))) 
                if self.args.verbose:
                    print(f"Monte Carlo ANISOTROPIC Barostat will be used. The dimensions that can change are: X = {self.scaleX} Y = {self.scaleY} Z = {self.scaleZ}")
        if self.sim_params.getBool('useAndersenThermostat'):    
            self.system.addForce(AndersenThermostat(self.temperature_c, self.col_freq))
            if self.args.verbose:
                print(f"Andersen Thermostat will be used with temperature {self.temperature_c}")

        #### Integrator ####
        friction_c    = self.sim_params.getFloat('friction_c')
        temperature_s = self.sim_params.getFloat('temperature_s')
        integrator    = self.sim_params.getStr('integrator')
        if "NoseHooverIntegrator" == integrator:
            self.integrator = NoseHooverIntegrator(self.temperature_c, friction_c, self.dt)
        elif "DrudeLangevinIntegrator" == integrator:
            self.integrator = DrudeLangevinIntegrator(self.temperature_c, friction_c, temperature_s, 
                                                      self.sim_params.getFloat('friction_s'), self.dt)
            print(DrudeLangevinIntegrator.getTemperature(self.integrator))
        elif "DrudeNoseHooverIntegrator" == integrator:
            self.integrator = DrudeNoseHooverIntegrator(self.temperature_c, friction_c, temperature_s, 
                                                        self.sim_params.getFloat('friction_s'), self.dt)
        elif "DrudeSCFIntegrator" == integrator:
            self.integrator = DrudeSCFIntegrator(self.dt)
            self.integrator.setDrudeTemperature(self.temperature_s)
        else:
            sys.exit("Unknown integrator %s" % integrator)
        if self.args.polarizable:
#            self.integrator.setMaxDrudeDistance(0.02*nanometer)
            self.integrator.setMaxDrudeDistance(self.MaxDrudeDist)
            if self.args.verbose:
                print("Core Temperature %g" % self.temperature_c)
                print("Drude Temperature %g" % self.integrator.getDrudeTemperature()._value)
                print("Step size %g" % self.integrator.getStepSize()._value)

    def compute_dipole(self)->list:
        positions = self.simulation.context.getState(getPositions=True).getPositions()
        dip = [ 0, 0, 0 ]
        enm2Debye = 48.0321
        for index in range(self.system.getNumParticles()):
            for m in range(3):
                dip[m] += positions[index][m]._value * self.charges[index] * enm2Debye
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
            if (not self.args.polarizable or index in self.cores):
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
        if self.nb_correction:
            self.qq_correction.updateParametersInContext(self.simulation.context)
            self.vdw_correction.updateParametersInContext(self.simulation.context)
        self.nb_openmm.updateParametersInContext(self.simulation.context)
        
    def print_energy(self, title:str):
        print("%s:" % title)
        etot = 0.0
        for group in self.force_group.keys():
            eterm = self.simulation.context.getState(getEnergy=True, groups=(1 << group)).getPotentialEnergy()/unit.kilojoule_per_mole
            etot += eterm
            print('%-40s : %16.4f kJ/mol' % (self.force_group[group], eterm))
        potE = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()/unit.kilojoule_per_mole
        print('Potential energy = %.2f kJ/mol' % potE)
        if None != self.args.emonomer:
            nmol = self.topology.getNumResidues()
            relener = potE/nmol - self.args.emonomer
            print('Interaction energy for %d-mer %g' % ( nmol, relener ))
            kB = 1.380649e-23 * 6.02214e23 / 1000
            print('Delta H vap %g kJ/mol' % ( kB*self.temperature_c - relener ))
        if abs(potE-etot) > 1e-3:
            print("sum of the above %.2f" % (etot))
        
    def minimize_energy(self, maxIter:int)->float:
        #### Minimize and Equilibrate ####
        if self.args.verbose:
            print('Performing energy minimization...')
        self.simulation.minimizeEnergy(maxIterations=maxIter)
        return self.simulation.context.getState(getEnergy=True).getPotentialEnergy()/unit.kilojoule_per_mole

    def equilibrate(self):
        print('Equilibrating...')
        self.simulation.context.setVelocitiesToTemperature(self.temperature_c)
        self.simulation.step(self.equilibrationSteps)
    
    def production(self):
        if self.args.verbose:
            simtime = self.sim_params.getFloat('dt')*self.sim_params.getInt('steps')
            print('Simulating %g ps...' % simtime )
        self.simulation.reporters.append(self.dcdReporter)
        self.simulation.reporters.append(self.dataReporter)
        self.simulation.reporters.append(self.pdbReporter)
        self.simulation.reporters.append(self.chkReporter)
        self.simulation.currentStep = 0
        self.simulation.step(self.steps)

    def setup(self):
        self.set_params()
        self.sim_params.check_consist()
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
            self.pdb.writeFile(self.topology,
                               self.simulation.context.getState(getPositions=True).getPositions(),
                               outf)
    def run(self):
        self.setup()
        self.minimize()
        self.equilibrate()
        self.print_energy("After equilibration")
        self.production()
        self.print_energy("After production")

    def log_to_xvg(self, xvg:str, ytargets:list):
        if None == self.logfile or not os.path.exists(self.logfile):
            print("Could not find any log file")
        else:
            xtarget  = "Time (ps)"
            ix = -1
            iy = []
            with open(xvg, "w") as outf:
                outf.write("@ xaxis label \"%s\"\n" % xtarget)
                with open(self.logfile, "r") as inf:
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
                                print("Incomprehensible line in logfile")
                                
