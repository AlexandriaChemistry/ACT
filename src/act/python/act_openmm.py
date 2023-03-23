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
import argparse, sys

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
    def __init__(self, qdist:str, comb:str):
        self.qdist = qdist
        self.comb  = comb
        
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
            csigma   = "(((sqrt(((epsilon1*gamma1*sigma1^6)/(gamma1-6)) * ((epsilon2*gamma2*sigma2^6)/(gamma2-6)))*(gamma-6))/(epsilon*gamma))^(1/6))"
            cepsilon = "((2 * epsilon1 * epsilon2)/(epsilon1 + epsilon2))"
            cgamma   = "((gamma1 + gamma2)/2)"
            return csigma, cepsilon, cgamma
        else:
            sys.exit("Unknown combination rule '%s'" % self.comb)
            
    def zetaString(self)->str:
        if self.qdist == "Gaussian":
            return ("(zeta1*zeta2/sqrt(zeta1^2+zeta2^2))")
        else:
            sys.exit("No support for charge distribution type %s" % self.qdist)
        
class ActOpenMMSim:
    def __init__(self):
        self.args        = self.argParser()
        self.pdb         = PDBFile(self.args.pdb_file)
        self.forcefield  = ForceField(self.args.xml_file)
        self.sim_params  = SimParams(self.args.dat_file)
        self.force_group = {}
        self.fgnumber    = {}
        self.comb        = CombinationRules(self.sim_params.getStr("charge_distribution"),
                                            self.sim_params.getStr("combination_rule"))
        self.logfile     = self.args.outputdir+'/'+self.args.log_file
        
    def add_force_group(self, force, fcname:str):
        fnumber = len(self.fgnumber)
        if fcname == 'NonbondedForce':
            directname = 'NonbondedForce (direct space)'
            self.force_group[fnumber] = directname
            self.fgnumber[directname] = fnumber
            force.setForceGroup(fnumber)
            fnumber = len(self.fgnumber)
            recipname = 'NonbondedForce (reciprocal space)'
            self.force_group[fnumber] = recipname
            self.fgnumber[recipname]  = fnumber
            force.setReciprocalSpaceForceGroup(fnumber)
        else:
            self.force_group[fnumber] = fcname
            self.fgnumber[fcname]     = fnumber
            force.setForceGroup(fnumber)
    
    def argParser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("-pdb", "--pdb_file", help="coordinate .pdb file", default=None)
        parser.add_argument("-xml", "--xml_file", help="openMM force field .xml file", default=None)
        parser.add_argument("-dat", "--dat_file", help="simulation parameter .dat file", default=None)
        deflog = "log.csv"
        parser.add_argument("-log", "--log_file", help="Log file for energies, default +"+deflog, default=deflog)
        parser.add_argument("-pol", "--polarizable", help="Turn on support for polarization", action="store_true")
        defout = "output"
        parser.add_argument("-odir", "--outputdir", help="Directory to write output to, default: "+defout, type=str, default=defout)
        defbonded = "morse"
        parser.add_argument("-bonds", "--bonded_potential", help="Either morse, or cubic, default "+defbonded, type=str, default=defbonded)
        defmDrude = "0.1"
        parser.add_argument("-mDrude", "--Drude_mass", help="The mass of drude particles. Make sure it is consistent with the mass in your forcefield"+defmDrude, type=float, default=defmDrude)
        parser.add_argument("-v", "--verbose", help="Print more stuff", action="store_true")
        # Parse the arguments
        args = parser.parse_args()
        # Do some checking
        if None == args.pdb_file or not os.path.exists(args.pdb_file):
            sys.exit("Please pass a correct pdb file")
        if None == args.xml_file or not os.path.exists(args.xml_file):
            sys.exit("Please pass a correct xml force field file")
        if None == args.dat_file or not os.path.exists(args.dat_file):
            sys.exit("Please pass a correct simulation parameter file")
    
        return args

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
        self.useAndersenThermostat     = self.sim_params.getBool('useAndersenThermostat')
        self.temperature_c             = self.sim_params.getFloat('temperature_c')
        self.useMonteCarloBarostat     = self.sim_params.getBool('useMonteCarloBarostat')
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
        os.makedirs(self.args.outputdir, exist_ok=True)
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
        if self.args.polarizable:
            self.rigidWater = False
            self.system    = self.forcefield.createSystem(self.topology,
                                                nonbondedMethod=self.nonbondedMethod,
                                                nonbondedCutoff=self.nonbondedCutoff,
                                                ewaldErrorTolerance=myEwaldErrorTolerance,
                                                rigidWater=self.rigidWater,
                                                drudeMass=myDrudeMass*unit.amu)
            if self.args.verbose:
                print(f"The force field is polarizable and the drude mass is {self.args.Drude_mass}. Make sure it is consistent with your force field file. The water is set to false.")
        else:
            self.system    = self.forcefield.createSystem(self.topology,
                                                      nonbondedMethod=self.nonbondedMethod,
                                                      nonbondedCutoff=self.nonbondedCutoff,
                                                      ewaldErrorTolerance=myEwaldErrorTolerance,
                                                      rigidWater=self.rigidWater)
            if self.args.verbose:
                print(f"The force field is NOT polarizable and the drude mass {self.args.Drude_mass} is not going to be used.")


        if self.args.verbose:
            print(f"We are working with {self.args.bonded_potential} form of bonded potential.\nIf the corresponding parameters are not found in force field file, the code will crash.")
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
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }
        self.reference_nb_force  = forces['NonbondedForce']
        self.reference_cnb_force = forces['CustomNonbondedForce']
        if self.args.verbose:
            print("***************************")
            print(f"Number of particles (incl. drudes):  {self.system.getNumParticles()}")
        if self.args.polarizable:
            reference_pol_force = forces['DrudeForce']
            self.cores = []
            self.shells = []
            self.core_shell = []
            for index in range(reference_pol_force.getNumParticles()):
                if self.args.verbose:
                    print(f"Polforce {reference_pol_force.getParticleParameters(index)}")
                [particle, particle1, particle2, particle3, particle4, charge, pol, aniso12, aniso34] = reference_pol_force.getParticleParameters(index)
                self.shells.append(particle) # particle  = shell
                self.cores.append(particle1) # particle1 = core
                self.core_shell.append((particle,particle1))
            if self.args.verbose:
                # Checking correct atom/shell pairing
                print(f"cores      {self.cores}")
                print(f"shells     {self.shells}")
                print(f"core_shell {self.core_shell}")
                print("########################")


        # TODO: Only for pbc simulations
        [alpha_ewald, nx, ny, nz] = self.reference_nb_force.getPMEParameters()
        if (alpha_ewald/alpha_ewald.unit) == 0.0:
            # If alpha is 0.0, alpha_ewald is computed by OpenMM from from the error tolerance.
            tol = self.reference_nb_force.getEwaldErrorTolerance()
            alpha_ewald = (1.0/self.reference_nb_force.getCutoffDistance()) * np.sqrt(-np.log(2.0*tol))
        cutoff_distance = self.reference_nb_force.getCutoffDistance()
        switch_distance = self.reference_nb_force.getSwitchingDistance()
    
        # Electrostatics
        expression = 'Coulomb_gauss - Coulomb_point;'
        self.qq_expression = ( "(%f*charge1*charge2*erf(zeta*r)/r)" % ONE_4PI_EPS0 )
        expression += ( 'Coulomb_gauss = %s;' % self.qq_expression )
        expression += ( 'Coulomb_point = (%f*charge1*charge2/r);' % ONE_4PI_EPS0 )
        expression += ( "zeta = %s;" % self.comb.zetaString())

        force = openmm.CustomNonbondedForce(expression)
        force.addPerParticleParameter("charge")
        force.addPerParticleParameter("zeta")
        force.setUseSwitchingFunction(self.use_switching_function)
        qqforce = openmm.CustomNonbondedForce(expression)
        qqforce.addPerParticleParameter("charge")
        qqforce.addPerParticleParameter("zeta")
        qqforce.setUseSwitchingFunction(self.use_switching_function)
        if self.nonbondedMethod == NoCutoff:
            qqforce.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
        else:    
            qqforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        qqforce.setCutoffDistance(cutoff_distance)
        qqforce.setSwitchingDistance(switch_distance)
        # Don't use dispersion correction for coulomb, it does not converge: https://github.com/openmm/openmm/issues/3162
        qqforce.setUseLongRangeCorrection(False)

        for index in range(self.reference_nb_force.getNumParticles()):
            [vdW, sigma, epsilon, gamma, charge, zeta] = self.reference_cnb_force.getParticleParameters(index)
            if self.args.verbose:
                print(f"nonbonded vdw sigma, epsilon, gamma, charge, zeta {self.reference_cnb_force.getParticleParameters(index)}")
            force.addParticle([charge,zeta])
            qqforce.addParticle([charge,zeta])
        for index in range(self.reference_nb_force.getNumExceptions()):
            [iatom, jatom, chargeprod, sigma, epsilon] = self.reference_nb_force.getExceptionParameters(index)
            qqforce.addExclusion(iatom, jatom)
        self.add_force_group(qqforce, "Coulomb")
        self.system.addForce(qqforce)
        # vdW
        expression = 'U_WKB- U_LJ;' 

        # TODO: Remove hard-coded combination rules if possible
        expression += 'U_LJ = 4*epsilon_LJ*((sigma_LJ/r)^12 -(sigma_LJ/r)^6);'
        expression += ('epsilon_LJ   = %s;' % self.comb.geometricString("epsilon_LJ1", "epsilon_LJ2"))
        expression += ('sigma_LJ     = %s;' % self.comb.arithmeticString("sigma_LJ1", "sigma_LJ2"))
        expression += ('sigma_LJ_rec = %s;' % self.comb.geometricString("sigma_LJ1", "sigma_LJ2"))
    
        self.vdw_expression =('(((((2*epsilon)/(1-(3/(gamma+3)))) * ((sigma^6)/(sigma^6+r^6))* (((3/(gamma+3))*(exp(gamma*(1-(r/sigma)))))-1))));')

        expression += ( 'U_WKB = %s;' % self.vdw_expression )
        #vdW*(((2*epsilon)/(1-(3/(gamma+3)))) * ((sigma^6)/(sigma^6+r^6))* (((3/(gamma+3))*(exp(gamma*(1-(r/(((sqrt(((epsilon1*gamma1*sigma1^6)/(gamma1-6)) * ((epsilon2*gamma2*sigma2^6)/(gamma2-6)))*(gamma-6))/(epsilon*gamma))^(1/6)))))))-1));'
        csigma, cepsilon, cgamma = self.comb.combStrings()
        # The statements have to be in this order! They are evaluated in the reverse order apparently.
        expression += ( 'sigma    = %s;' % csigma )
        expression += ( 'epsilon  = %s;' % cepsilon )
        expression += ( 'gamma    = %s;' % cgamma )
        expression += 'vdW = vdW1*vdW2;'
        vdwforce = openmm.CustomNonbondedForce(expression)
        vdwforce.addPerParticleParameter("sigma")
        vdwforce.addPerParticleParameter("epsilon")
        vdwforce.addPerParticleParameter("gamma")
        vdwforce.addPerParticleParameter("vdW")
        vdwforce.addPerParticleParameter("sigma_LJ")
        vdwforce.addPerParticleParameter("epsilon_LJ")
        vdwforce.setUseSwitchingFunction(self.use_switching_function)
        if self.nonbondedMethod == NoCutoff:
            vdwforce.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
        else:    
            vdwforce.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        vdwforce.setCutoffDistance(cutoff_distance)
        vdwforce.setSwitchingDistance(switch_distance)
        vdwforce.setUseLongRangeCorrection(self.reference_nb_force.getUseDispersionCorrection())
        for index in range(self.reference_nb_force.getNumParticles()):
            [charge_LJ, sigma_LJ, epsilon_LJ] = self.reference_nb_force.getParticleParameters(index)
            [vdW, sigma, epsilon, gamma, charge, zeta] = self.reference_cnb_force.getParticleParameters(index)
            force.addParticle([sigma, epsilon, gamma, vdW, sigma_LJ, epsilon_LJ])
            vdwforce.addParticle([sigma, epsilon, gamma, vdW, sigma_LJ, epsilon_LJ])
            if self.args.verbose:
                print("index %d sigma %g, epsilon %g, gamma %g, vdW %g, sigma_LJ %g, epsilon_LJ %g" %  (index, sigma, epsilon, gamma, vdW, sigma_LJ._value, epsilon_LJ._value ))
        for index in range(self.reference_nb_force.getNumExceptions()):
            [iatom, jatom, chargeprod, sigma, epsilon] = self.reference_nb_force.getExceptionParameters(index)
            vdwforce.addExclusion(iatom, jatom)
            if self.args.verbose:
                print("excl %d iatom %d jatom %d" % ( index, iatom, jatom ))
        self.add_force_group(vdwforce, "Wang-Buckingham minus Lennard-Jones")
        self.system.addForce(vdwforce)

    def real_exclusion(self, nexcl:int, iatom:int, jatom:int)->bool:
        if nexcl == 0:
            return False
        elif nexcl == 1:
            return ((iatom,jatom) in self.bonds or (jatom,iatom) in self.bonds)
        else:
            sys.exit("Cannot handle nexcl == %d" % nexcl)
        return False

    def add_excl_correction(self):
        # Add vdW and electrostactics that have been excluded (this has to be done as the number of exclusions is 3 for nonbonded interactions in OpenMM)
        # Those interactions are added using two CustomBondForce entries
        vdw_force = openmm.CustomBondForce(self.vdw_expression)
        vdw_force.addPerBondParameter("sigma")
        vdw_force.addPerBondParameter("epsilon")
        vdw_force.addPerBondParameter("gamma")
        
        qq_force = openmm.CustomBondForce(self.qq_expression)
        qq_force.addPerBondParameter("charge1")
        qq_force.addPerBondParameter("charge2")
        qq_force.addPerBondParameter("zeta")

        nexclvdw = self.sim_params.getInt("nexclvdw")
        nexclqq  = self.sim_params.getInt("nexclqq")
        
        for index in range(self.reference_nb_force.getNumExceptions()):
            # Just get the excluded atoms from the regular NB force
            [iatom, jatom, chargeprod_except, sigma_except, epsilon_except] = self.reference_nb_force.getExceptionParameters(index)
            # Check for shell exclusions first
            if (self.args.polarizable and 
                ((jatom,iatom) in self.core_shell or ((iatom,jatom) in self.core_shell))):
                continue
            # And get the parameters from the Custom NB force
            [vdW1, sigma1, epsilon1, gamma1, charge1, zeta1] = self.reference_cnb_force.getParticleParameters(iatom)
            [vdW2, sigma2, epsilon2, gamma2, charge2, zeta2] = self.reference_cnb_force.getParticleParameters(jatom)
            if self.args.verbose:
                print(f" custom nonbonded force i {self.reference_cnb_force.getParticleParameters(iatom)}")
                print(f" custom nonbonded force j {self.reference_cnb_force.getParticleParameters(jatom)}")

            # Coulomb part
            if not self.real_exclusion(nexclqq, iatom, jatom):
                zeta = ((zeta1 * zeta2)/(np.sqrt(zeta1**2 + zeta2**2)))
                qq_force.addBond(iatom, jatom, [charge1, charge2, zeta])
                
            # Van der Waals part
            if not self.real_exclusion(nexclvdw, iatom, jatom):
                # Note that the 6 below is because of the combination rule
                if epsilon1 > 0 and epsilon2 > 0 and gamma1 > 6 and gamma2 > 6:
                    csigma, cepsilon, cgamma = self.comb.combStrings()
                    if self.args.verbose:
                        print(csigma)
                        print(cepsilon)
                        print(cgamma)
                    gamma   = eval(cgamma)
                    epsilon = eval(cepsilon)
                    sigma   = eval(csigma.replace("^", "**"))

                    if self.args.verbose:
                        print("i %d j %d q1 %g q2 %g sigma %g epsilon %g gamma %g zeta1 %g zeta2 %g" % 
                              ( iatom, jatom, charge1, charge2, sigma, epsilon, gamma, zeta1, zeta2 ))
                    if vdW1 != 0 and vdW2 != 0:
                        vdw_force.addBond(iatom, jatom, [sigma, epsilon, gamma])
                        
        self.add_force_group(qq_force, "Coulomb Exclusion Correction")
        self.system.addForce(qq_force)
        self.add_force_group(vdw_force, "Van der Waals Exclusion Correction")
        self.system.addForce(vdw_force)
   
    def add_bonded_forces(self):
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }
        reference_cb_force  = forces['CustomBondForce']
        self.bonds = []
        if self.args.bonded_potential == "morse":
            ### Morse potential ###
            Morse_expression = "(D_e*(1 - exp(-beta*(r-r0)))^2)+D0;"
            Morse_force = openmm.CustomBondForce(Morse_expression)
            Morse_force.addPerBondParameter("beta")
            Morse_force.addPerBondParameter("D_e")
            Morse_force.addPerBondParameter("D0")
            Morse_force.addPerBondParameter("r0")
            for bond_index in range(reference_cb_force.getNumBonds()):
                # Retrieve parameters.
                [iatom, jatom, (beta, D_e, D0, r0)] = reference_cb_force.getBondParameters(bond_index)
                if self.args.verbose:
                    print("Morse: iatom %d jatom %d beta %g D_e %g D0 %g r0 %g" % ( iatom, jatom, beta, D_e, D0, r0 ))
                self.bonds.append((iatom, jatom))
                Morse_force.addBond(iatom, jatom, [beta, D_e, D0, r0])
            self.add_force_group(Morse_force, "Morse bonds")
            self.system.addForce(Morse_force)
            
        ###################################################    
        elif self.args.bonded_potential == "cubic":
            ### Cubic potential ###
            # Maximum in the potential is at r = (2*rmax+bondlength)/3
            # The final value is
            # V((2*rmax+bondlength)/3) = -D_e - (4/27)kb(bondlength -rmax)^2
            Cubic_expression = "(kb*(r-bondlength)^2 * (rmax-r) - D_e)*step(2*rmax+bondlength-3*r) + step(3*r-2*rmax-bondlength)*(-D_e - (4*kb/27)*(bondlength-rmax)^2);" 
            Cubic_force = openmm.CustomBondForce(Cubic_expression)
            Cubic_force.addPerBondParameter("bondlength")
            Cubic_force.addPerBondParameter("rmax")
            Cubic_force.addPerBondParameter("kb")
            Cubic_force.addPerBondParameter("D_e")
            for bond_index in range(reference_cb_force.getNumBonds()):
                # Retrieve parameters.
                [iatom, jatom, (bondlength, rmax, kb, D_e)] = reference_cb_force.getBondParameters(bond_index)
                self.bonds.append((iatom, jatom))
                Cubic_force.addBond(iatom, jatom, [bondlength, rmax, kb, D_e])
            self.add_force_group(Cubic_force, "Cubic bonds")
            self.system.addForce(Cubic_force)
        else:
            sys.exit("Unsupported bonded function %s" % self.args.bonded_potential)

    def make_forces(self):
        # Create a new CustomNonbondedForce to mimic the direct space 
        self.add_direct_space_force()
        self.add_bonded_forces()
        self.add_excl_correction()

        forces = { force.__class__.__name__ : force for force in self.system.getForces() }
        self.Cnbforce   = forces['CustomNonbondedForce']
        self.CBondforce = forces['CustomBondForce']
        self.nbforce    = forces['NonbondedForce']
        self.add_force_group(self.nbforce, "NonbondedForce")
        #### for systems with core shell particles in general; nbforce.getNumParticles() is cores+shells
        #### TODO: get the real number of particle
        self.Numb_particles = self.nbforce.getNumParticles()
        if self.args.polarizable:
            self.add_force_group(forces['DrudeForce'], "Polarization")

    def print_params(self):
        print("----------------------------")
        print("### CustomNonbondedForce ###")
        print('"Cutoff?" {0}'.format(self.Cnbforce.getCutoffDistance()))
        print('"SwitchingDistance?" {0}'.format(self.Cnbforce.getSwitchingDistance ()))
        print('"CustomNonbondedMethod?" {0}'.format(self.Cnbforce.getNonbondedMethod()))
        print('"SwitchingFunction?" {0}'.format(self.Cnbforce.getUseSwitchingFunction()))
        print('"CustomNB PBC?" {0}'.format(self.Cnbforce.usesPeriodicBoundaryConditions()))
        print('"Force Group?" {0}'.format(self.Cnbforce.getForceGroup()))

        print("----------------------------")
        print("### NonbondedForce ###")
        print('"Cutoff?" {0}'.format(self.nbforce.getCutoffDistance()))
        print('"SwitchingDistance?" {0}'.format(self.nbforce.getSwitchingDistance ()))
        print('"NonbondedMethod?" {0}'.format(self.nbforce.getNonbondedMethod()))
        print('"SwitchingFunction?" {0}'.format(self.nbforce.getUseSwitchingFunction()))
        print('"NB PBC?" {0}'.format(self.nbforce.usesPeriodicBoundaryConditions()))
        print('"Disp. Corr.?" {0}'.format(self.nbforce.getUseDispersionCorrection()))
        print('"Force Group?" {0}'.format(self.nbforce.getForceGroup()))
        print('"Reciprocal Force Group?" {0}'.format(self.nbforce.getReciprocalSpaceForceGroup()))
        print("----------------------------")

    def set_algorithms(self):
        #### Thermostat / Barostat ####
        if self.nonbondedMethod != NoCutoff:
            if self.sim_params.getBool('useMonteCarloBarostat'):
                if self.args.verbose:
                    print("Monte Carlo Barostat shall be used, as was asked for.")
                    self.system.addForce(MonteCarloBarostat(self.sim_params.getFloat('pressure'),
                                                    self.temperature_c,
                                                    self.sim_params.getInt('barostatInterval')))
            elif self.sim_params.getBool('useAndersenThermostat'):    
                self.system.addForce(AndersenThermostat(self.temperature_c, self.col_freq))
                if self.args.verbose:
                    print(f"Andersen Thermostat shall be used, as was asked for. With temperature {self.temperature_c}")
            else:
                print("I shall refrain from using a Thermostat...")

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
            self.integrator.setMaxDrudeDistance(0.02*nanometer)
            if self.args.verbose:
                print("Core Temperature %g" % self.temperature_c)
                print("Drude Temperature %g" % self.integrator.getDrudeTemperature()._value)
                print("Step size %g" % self.integrator.getStepSize()._value)

    def init_simulation(self):
        #### Simulation setup ####
        self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform)
        self.simulation.context.setPositions(self.positions)

        #### Set positions of shell systemot zero) ####
        #### the shell displacement is necessary for the LJPME to work, otherwise an error is thrown:
        #### simtk.openmm.OpenMMException: Particle coordinate is nan
        positions = self.simulation.context.getState(getPositions=True).getPositions()
        new_pos = []
        ##print(new_pos)
        for index in range(self.system.getNumParticles()):
            if (not self.args.polarizable or index in self.cores):
                new_pos_x = positions[index][0]+0.000*nanometer
                new_pos.append((new_pos_x,positions[index][1],positions[index][2]))
            if (self.args.polarizable and index in self.shells):
                new_pos_x = positions[index][0]+0.01*nanometer
                new_pos_y = positions[index][1]+0.01*nanometer
                new_pos_z = positions[index][2]+0.01*nanometer
                new_pos.append((new_pos_x,new_pos_y,new_pos_z)) 
        if self.args.verbose:
            print(f"number of particles (incl. drudes):  {self.system.getNumParticles()}")
            for np in new_pos:
                print("%10.5f  %10.5f  %10.5f" % ( np[0]._value, np[1]._value, np[2]._value ))
        self.simulation.context.setPositions(new_pos)
        self.Cnbforce.updateParametersInContext(self.simulation.context)
        self.nbforce.updateParametersInContext(self.simulation.context)
        
    def print_energy(self, title:str):
        print("%s:" % title)
        etot = 0.0
        for group in range(len(self.force_group.keys())):
            eterm = self.simulation.context.getState(getEnergy=True, groups=(1 << group)).getPotentialEnergy()/unit.kilojoule_per_mole
            etot += eterm
            print('%8d : %64s : %16.4f kJ/mol' % (group, self.force_group[group], eterm))
        potE = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()/unit.kilojoule_per_mole
        print('potential energy = %.2f kJ/mol (sum of the above %.2f)\n' % (potE, etot))    
        
    def minimize_energy(self):
        #### Minimize and Equilibrate ####
        if self.args.verbose:
            print('Performing energy minimization...')
        self.simulation.minimizeEnergy()

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

    def run(self):
        self.set_params()
        self.start_output()
        self.make_system()
        self.make_forces()
        self.print_params()
        self.set_algorithms()
        self.init_simulation()
        self.print_energy("Initial energies:")
        self.minimize_energy()
        self.print_energy("After minimization:")
        self.equilibrate()
        self.print_energy("After equilibration:")
        self.production()
        self.print_energy("After production:")

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
                                
