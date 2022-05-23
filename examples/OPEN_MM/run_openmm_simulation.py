#!/usr/bin/env python3

# OpenMM python example script (using DrudeLangevin mintegrator). 
# This script implements a modified Buckingham potential (using Hogervorst combination rules) and Gaussian distributed charges for the nonbonded 
# interactions, and a Morse potential for the bonded interactions.
############################################################ Proceed at your own risk. ############################################################
# Van der Spoel Group
# Author: Marie-Madeleine Walz, Department of Cell and Molecular Biology, Uppsala University, Sweden. marie-madeleine.walz@icm.uu.se
###################################################################################################################################################

from openmm import *
from openmm.app import *
from simtk.unit import *
from simtk import openmm, unit
import numpy as np
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument("-pdb", "--pdb_file", help="coordinate .pdb file", default=None)
parser.add_argument("-xml", "--xml_file", help="openMM force field .xml file", default=None)
parser.add_argument("-dat", "--dat_file", help="simulation parameter .dat file", default=None)
parser.add_argument("-pol", "--polarizable", help="Turn on support for polarization", action="store_true")
defout = "output"
parser.add_argument("-odir", "--outputdir", help="Directory to write output to, default: "+defout, type=str, default=defout)

args = parser.parse_args()

if None == args.pdb_file or not os.path.exists(args.pdb_file):
    sys.exit("Please pass a correct pdb file")
pdb          = PDBFile(args.pdb_file)
if None == args.xml_file or not os.path.exists(args.xml_file):
    sys.exit("Please pass a correct xml force field file")
forcefield   = ForceField(args.xml_file)
if None == args.dat_file or not os.path.exists(args.dat_file):
    sys.exit("Please pass a correct simulation parameter file")
sim_dat_file = args.dat_file


def make_sim_dat_dict(sim_dat_file):
    sim_params = {}
    inFileStream = open(sim_dat_file, 'r')
    for line in inFileStream:
        try:
            key, equal, value =line.split() 
            sim_params[key] = value
        except:
            continue
    inFileStream.close() 
    return sim_params            


def get_value(sim_params,key_x):
    value_x = sim_params[key_x]
    return value_x
                  

def return_bool(value):
    if (value == "True"):
        return True
    else:
        return False     
         
sim_dat_dict = make_sim_dat_dict(sim_dat_file)

sim_dat_dict_internal = {'LJPME':LJPME, 'PME':PME, 'Ewald':Ewald, 'CutoffPeriodic':CutoffPeriodic, 'NoCutoff':NoCutoff, 'HBonds':HBonds, 'HAngles':HAngles, 'None':None}

# SET SIMULATION PARAMETERS
################################################
dt                 = float(get_value(sim_dat_dict, 'dt').split('*')[0])
equilibrationSteps = int(get_value(sim_dat_dict, 'equilibrationSteps'))
steps              = int(get_value(sim_dat_dict, 'steps'))

save               = int(get_value(sim_dat_dict, 'save'))
outStep            = return_bool(get_value(sim_dat_dict, 'outStep'))   
outTime            = return_bool(get_value(sim_dat_dict, 'outTime'))
outSpeed           = return_bool(get_value(sim_dat_dict, 'outSpeed'))
outProgress        = return_bool(get_value(sim_dat_dict, 'outProgress'))
outPotentialEnergy = return_bool(get_value(sim_dat_dict, 'outPotentialEnergy'))
outKineticEnergy   = return_bool(get_value(sim_dat_dict, 'outKineticEnergy'))
outTemperature     = return_bool(get_value(sim_dat_dict, 'outTemperature'))
outVolume          = return_bool(get_value(sim_dat_dict, 'outVolume'))
outDensity         = return_bool(get_value(sim_dat_dict, 'outDensity'))
outSeparator       = get_value(sim_dat_dict, 'outSeparator')

nonbondedMethod           = sim_dat_dict_internal[get_value(sim_dat_dict, 'nonbondedMethod')]
use_switching_function    = return_bool(get_value(sim_dat_dict, 'use_switching_function'))
switch_width              = float(get_value(sim_dat_dict, 'switch_width').split('*')[0])
nonbondedCutoff           = float(get_value(sim_dat_dict, 'nonbondedCutoff').split('*')[0])
use_dispersion_correction = return_bool(get_value(sim_dat_dict, 'use_dispersion_correction'))
ewaldErrorTolerance       = float(get_value(sim_dat_dict, 'ewaldErrorTolerance')) 

useAndersenThermostat     = return_bool(get_value(sim_dat_dict, 'useAndersenThermostat'))
temperature_c             = float(get_value(sim_dat_dict, 'temperature_c').split('*')[0]) 
temperature_s             = float(get_value(sim_dat_dict, 'temperature_s').split('*')[0]) 
useMonteCarloBarostat     = return_bool(get_value(sim_dat_dict, 'useMonteCarloBarostat')) 
friction_c                = float(get_value(sim_dat_dict, 'friction_c').split('/')[0]) 
friction_s                = float(get_value(sim_dat_dict, 'friction_s').split('/')[0])  
pressure                  = float(get_value(sim_dat_dict, 'pressure').split('*')[0])  
barostatInterval          = int(get_value(sim_dat_dict, 'barostatInterval')) 

constraints            = sim_dat_dict_internal[get_value(sim_dat_dict, 'constraints')]
rigidWater             = return_bool(get_value(sim_dat_dict, 'rigidWater'))
constraintTolerance    = float(get_value(sim_dat_dict, 'constraintTolerance'))

usePlatform      = get_value(sim_dat_dict, 'usePlatform') 
usePrecisionCuda = get_value(sim_dat_dict, 'usePrecisionCuda') 




# COMPUTING PLATFORM
################################################
#platform = Platform.getPlatformByName('CUDA')
#properties = {'CudaPrecision': 'single'}
platform = Platform.getPlatformByName(usePlatform)

# OUTPUT
################################################
os.makedirs(args.outputdir, exist_ok=True)
dcdReporter = DCDReporter(args.outputdir+'/trajectory.dcd', save)
dataReporter = StateDataReporter(args.outputdir+'/log.csv', save, totalSteps=steps,
    step=outStep, time= outTime, speed=outSpeed, progress=outProgress, potentialEnergy=outPotentialEnergy, kineticEnergy=outKineticEnergy, temperature=outTemperature, volume=outVolume, density=outDensity, separator=outSeparator)

chkReporter = CheckpointReporter(args.outputdir+'/checkpnt.chk', save)
pdbReporter = PDBReporter(args.outputdir+'/output.pdb', save)


# TOPOLOGY
################################################
topology  = pdb.topology
positions = pdb.positions
modeller  = Modeller(topology, positions)
modeller.addExtraParticles(forcefield)
topology  = modeller.topology
positions = modeller.positions
system    = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,ewaldErrorTolerance=ewaldErrorTolerance, rigidWater=rigidWater)


# SETTINGS FOR FORCES
################################################

for force in system.getForces():
    if hasattr(force, 'setCutoffDistance'):
        force.setCutoffDistance(nonbondedCutoff)
    if hasattr(force, 'setUseSwitchingFunction'):
        force.setUseSwitchingFunction(use_switching_function)
    if hasattr(force, 'setSwitchingDistance'):
        switch_distance = nonbondedCutoff-switch_width
        force.setSwitchingDistance(switch_distance)
    if hasattr(force, 'setEwaldErrorTolerance'):
        force.setEwaldErrorTolerance(ewaldErrorTolerance)
    if hasattr(force, 'setUseDispersionCorrection'):
        force.setUseDispersionCorrection(use_dispersion_correction)


# CODE FOR ALEXANDRIA NONBONDED FORCES
################################################
def add_direct_space_force(system, group_cnb, group_cb): 
    """
    Create a CustomNonbondedForce to calculate the direct-space force of WBK-LJ and gaussian Coulomb-point Coulomb, placing it in specified force group.
    The LJ and point charge is necessary for both the dispersion correction and for the LJPME, and for using PME
    Create a CustomBondForce to calculate the direct space force of WBK and gaussian Coulomb for interactions that are excluded (besides core-shell interactions).

    """
    forces = { force.__class__.__name__ : force for force in system.getForces() }
    reference_nb_force  = forces['NonbondedForce']
    reference_cnb_force = forces['CustomNonbondedForce']
    reference_cb_force  = forces['CustomBondForce']
    if args.polarizable:
        reference_pol_force = forces['DrudeForce']
        cores = []
        shells = []
        core_shell=[]
        for index in range(reference_pol_force.getNumParticles()):
            [particle, particle1, particle2, particle3, particle4, charge, pol, aniso12, aniso34] = reference_pol_force.getParticleParameters(index)
            shells.append(particle) # particle  = shell
            cores.append(particle1) # particle1 = core
            core_shell.append((particle,particle1))

    ONE_4PI_EPS0 = 138.935456
    [alpha_ewald, nx, ny, nz] = reference_nb_force.getPMEParameters()
    if (alpha_ewald/alpha_ewald.unit) == 0.0:
        # If alpha is 0.0, alpha_ewald is computed by OpenMM from from the error tolerance.
        tol = reference_nb_force.getEwaldErrorTolerance()
        alpha_ewald = (1.0/reference_nb_force.getCutoffDistance()) * np.sqrt(-np.log(2.0*tol))
    cutoff_distance = reference_nb_force.getCutoffDistance()
    switch_distance = reference_nb_force.getSwitchingDistance()
    
    # Electrostatics
    expression = 'Coulomb_gauss - Coulomb_point;'
    expression += 'Coulomb_gauss = (ONE_4PI_EPS0*charge1*charge2*erf(beta*r)/r);'
    expression += 'Coulomb_point = (ONE_4PI_EPS0*charge1*charge2/r);'
    expression += 'beta = ((beta1 * beta2)/(sqrt(beta1^2 + beta2^2)));'
    expression += 'ONE_4PI_EPS0 = %.16e;' % (ONE_4PI_EPS0)
    force = openmm.CustomNonbondedForce(expression)
    force.addPerParticleParameter("charge")
    force.addPerParticleParameter("beta")
    force.setUseSwitchingFunction(use_switching_function)
    if nonbondedMethod == NoCutoff:
        force.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
    else:    
        force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    force.setCutoffDistance(cutoff_distance)
    force.setSwitchingDistance(switch_distance)
    force.setUseLongRangeCorrection(False) #Don't use dispersion correction for coulomb, it does not converge: https://github.com/openmm/openmm/issues/3162

    for index in range(reference_nb_force.getNumParticles()):
        [vdW, sigma, epsilon, gamma, charge, beta] = reference_cnb_force.getParticleParameters(index)
        force.addParticle([charge,beta])
    for index in range(reference_nb_force.getNumExceptions()):
        [iatom, jatom, chargeprod, sigma, epsilon] = reference_nb_force.getExceptionParameters(index)
        force.addExclusion(iatom, jatom)
    force.setForceGroup(group_cnb)
    system.addForce(force)

    # vdW
    expression = 'U_WKB- U_LJ;' 

    expression += 'U_LJ = 4*epsilon_LJ*((sigma_LJ/r)^12 -(sigma_LJ/r)^6);'  
    expression += 'epsilon_LJ = sqrt(epsilon_LJ1*epsilon_LJ2);'
    expression += 'sigma_LJ = 0.5*(sigma_LJ1+sigma_LJ2);'
    expression += 'sigma_LJ_rec = sqrt(sigma_LJ1*sigma_LJ2);'
    
    expression += 'U_WKB = (((((2*epsilon)/(1-(3/(gamma+3)))) * (((((sqrt(((epsilon1*gamma1*sigma1^6)/(gamma1-6)) * ((epsilon2*gamma2*sigma2^6)/(gamma2-6)))*(gamma-6))/(epsilon*gamma))^(1/6))^6)/((((sqrt(((epsilon1*gamma1*sigma1^6)/(gamma1-6)) * ((epsilon2*gamma2*sigma2^6)/(gamma2-6)))*(gamma-6))/(epsilon*gamma))^(1/6))^6+r^6))* (((3/(gamma+3))*(exp(gamma*(1-(r/(((sqrt(((epsilon1*gamma1*sigma1^6)/(gamma1-6)) * ((epsilon2*gamma2*sigma2^6)/(gamma2-6)))*(gamma-6))/(epsilon*gamma))^(1/6)))))))-1))*vdW));'
    expression += 'epsilon = ((2 * epsilon1 * epsilon2)/(epsilon1 + epsilon2));'
    expression += 'gamma = ((gamma1 + gamma2)/2);'
    expression += 'vdW = vdW1*vdW2;'
    force = openmm.CustomNonbondedForce(expression)
    force.addPerParticleParameter("sigma")
    force.addPerParticleParameter("epsilon")
    force.addPerParticleParameter("gamma")
    force.addPerParticleParameter("vdW")
    force.addPerParticleParameter("sigma_LJ")
    force.addPerParticleParameter("epsilon_LJ")
    force.setUseSwitchingFunction(use_switching_function)#True
    if nonbondedMethod == NoCutoff:
        force.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
    else:    
        force.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    force.setCutoffDistance(cutoff_distance)
    force.setSwitchingDistance(switch_distance)
    force.setUseLongRangeCorrection(reference_nb_force.getUseDispersionCorrection())
    for index in range(reference_nb_force.getNumParticles()):
        [charge_LJ, sigma_LJ, epsilon_LJ] = reference_nb_force.getParticleParameters(index)
        [vdW, sigma, epsilon, gamma, charge, beta] = reference_cnb_force.getParticleParameters(index)
        force.addParticle([sigma, epsilon, gamma, vdW, sigma_LJ, epsilon_LJ])
    for index in range(reference_nb_force.getNumExceptions()):
        [iatom, jatom, chargeprod, sigma, epsilon] = reference_nb_force.getExceptionParameters(index)
        force.addExclusion(iatom, jatom)
    force.setForceGroup(group_cnb)
    system.addForce(force)


    # Add vdW and electrostactics that have been excluded (this has to be done as the number of exclusions is 3 for nonbonded interactions in OpenMM)
    # Those interactions are added using a CustomBondForce
    bond_expression =('(U_sterics+U_electrostatics);'
                      'U_sterics = (((((2*epsilon)/(1-(3/(gamma+3)))) * ((sigma^6)/(sigma^6+r^6))* (((3/(gamma+3))*(exp(gamma*(1-(r/sigma)))))-1))*vdW));'
                  'U_electrostatics = (ONE_4PI_EPS0*chargeprod* erf(beta*r)/r);'
                  )
    bond_expression += 'ONE_4PI_EPS0 = %.16e;' % (ONE_4PI_EPS0)

    bond_force = openmm.CustomBondForce(bond_expression)
    bond_force.addPerBondParameter("chargeprod")
    bond_force.addPerBondParameter("beta")
    bond_force.addPerBondParameter("sigma")
    bond_force.addPerBondParameter("epsilon")
    bond_force.addPerBondParameter("gamma")
    bond_force.addPerBondParameter("vdW")
    for index in range(reference_nb_force.getNumExceptions()):
        [iatom, jatom, chargeprod_except, sigma_except, epsilon_except] = reference_nb_force.getExceptionParameters(index)
        [vdW1, sigma1, epsilon1, gamma1, charge1, beta1] = reference_cnb_force.getParticleParameters(iatom)
        [vdW2, sigma2, epsilon2, gamma2, charge2, beta2] = reference_cnb_force.getParticleParameters(jatom)
        chargeprod=charge1*charge2
        beta = ((beta1 * beta2)/(np.sqrt(beta1**2 + beta2**2)))
        epsilon = ((2 * epsilon1 * epsilon2)/(epsilon1 + epsilon2))
        gamma = ((gamma1 + gamma2)/2)
        sigma = (((sqrt(((epsilon1*gamma1*sigma1**6)/(gamma1-6)) * ((epsilon2*gamma2*sigma2**6)/(gamma2-6)))*(gamma-6))/(epsilon*gamma))**(1/6))
        vdW = vdW1*vdW2
        if args.polarizable:
            if ((jatom,iatom) not in core_shell) and ((iatom,jatom) not in core_shell): 
                bond_force.addBond(iatom, jatom, [chargeprod, beta, sigma, epsilon, gamma, vdW])
    bond_force.setForceGroup(group_cb)
    system.addForce(bond_force)
    
     

    ### Morse potential ###
    Morse_expression = "(D_e*(1 - exp(-a*(r-r0)))^2)-D_e;"
    Morse_force = openmm.CustomBondForce(Morse_expression)
    Morse_force.addPerBondParameter("D_e")
    Morse_force.addPerBondParameter("a")
    Morse_force.addPerBondParameter("r0")
    
    for bond_index in range(reference_cb_force.getNumBonds()):
        # Retrieve parameters.
        [iatom, jatom, (D_e, a, r0)] = reference_cb_force.getBondParameters(bond_index)  
        Morse_force.addBond(iatom, jatom, [D_e, a, r0])
    Morse_force.setForceGroup(group_cb)
    system.addForce(Morse_force)
    


## Print energy components
force_group_names = {
    0 : 'HarmonicBondForce',
    1 : 'HarmonicAngleForce',
    2 : 'CustomNonbondedForce (direct space)',
    3 : 'NonbondedForce (direct space)',
    4 : 'NonbondedForce (reciprocal space)',
    5 : 'CustomBondForce'
}
if args.polarizable:
    force_group_names[6] = 'DrudeForce'
fgnumbers = {}
for f in force_group_names.keys():
    fgnumbers[force_group_names[f]] = f

# Set force groups:
for force in system.getForces():
    fcname = force.__class__.__name__
    if fcname in fgnumbers:
        force.setForceGroup(fgnumbers[fcname])
    if force.__class__.__name__ == 'NonbondedForce':
       force.setForceGroup(fgnumbers['NonbondedForce (direct space)'])
       force.setReciprocalSpaceForceGroup(fgnumbers['NonbondedForce (reciprocal space)'])


# Create a new CustomNonbondedForce to mimic the direct space 
add_direct_space_force(system, group_cnb=2, group_cb=fgnumbers['CustomBondForce'])

forces = { force.__class__.__name__ : force for force in system.getForces() }
Cnbforce = forces['CustomNonbondedForce']
CBondforce = forces['CustomBondForce']
nbforce  = forces['NonbondedForce']
if args.polarizable:
    polforce = forces['DrudeForce']
    cores = []
    shells = []
    core_shell=[]
    for index in range(polforce.getNumParticles()):
        [particle, particle1, particle2, particle3, particle4, charge, pol, aniso12, aniso34] = polforce.getParticleParameters(index)
        shells.append(particle) # particle  = shell
        cores.append(particle1) # particle1 = core
        core_shell.append([particle,particle1])
  

print("----------------------------")
print("----------------------------")
print("### CustomNonbondedForce ###")
print('"Cutoff?" {0}'.format(Cnbforce.getCutoffDistance()))
print('"SwitchingDistance?" {0}'.format(Cnbforce.getSwitchingDistance ()))
print('"CustomNonbondedMethod?" {0}'.format(Cnbforce.getNonbondedMethod()))
print('"SwitchingFunction?" {0}'.format(Cnbforce.getUseSwitchingFunction()))
print('"CustomNB PBC?" {0}'.format(Cnbforce.usesPeriodicBoundaryConditions()))
print('"Force Group?" {0}'.format(Cnbforce.getForceGroup()))


print("----------------------------")
print("----------------------------")
print("### NonbondedForce ###")
print('"Cutoff?" {0}'.format(nbforce.getCutoffDistance()))
print('"SwitchingDistance?" {0}'.format(nbforce.getSwitchingDistance ()))
print('"NonbondedMethod?" {0}'.format(nbforce.getNonbondedMethod()))
print('"SwitchingFunction?" {0}'.format(nbforce.getUseSwitchingFunction()))
print('"NB PBC?" {0}'.format(nbforce.usesPeriodicBoundaryConditions()))
print('"Disp. Corr.?" {0}'.format(nbforce.getUseDispersionCorrection()))
print('"Force Group?" {0}'.format(nbforce.getForceGroup()))
print('"Reciprocal Force Group?" {0}'.format(nbforce.getReciprocalSpaceForceGroup()))
print("----------------------------")
print("----------------------------")


#### Thermostat / Barostat ####
if nonbondedMethod != NoCutoff:
  system.addForce(MonteCarloBarostat(pressure, temperature_c, barostatInterval))
  system.addForce(AndersenThermostat(temperature_c, 1/picosecond))

#### Integrator ####
if args.polarizable:
    integrator = DrudeLangevinIntegrator(temperature_c, friction_c, temperature_s, friction_s, dt)
else:
    integrator = NoseHooverIntegrator(temperature_c, friction_c, dt)

#### Simulation setup ####
simulation = Simulation(topology, system, integrator, platform)

simulation.context.setPositions(positions)

#### Set positions of shell particles (so that r_ij is not zero) ####
#### the shell displacement is necessary for the LJPME to work, otherwise an error is thrown:
#### simtk.openmm.OpenMMException: Particle coordinate is nan
positions = simulation.context.getState(getPositions=True).getPositions()
new_pos = []
##print(new_pos)
for index in range(system.getNumParticles()):
  if (not args.polarizable or index in cores):
    new_pos_x = positions[index][0]+0.000*nanometer
    new_pos.append((new_pos_x,positions[index][1],positions[index][2]))
  if (args.polarizable and index in shells):
    new_pos_x = positions[index][0]+0.01*nanometer
    new_pos_y = positions[index][1]+0.01*nanometer
    new_pos_z = positions[index][2]+0.01*nanometer
    new_pos.append((new_pos_x,new_pos_y,new_pos_z)) 
simulation.context.setPositions(new_pos)
Cnbforce.updateParametersInContext(simulation.context)


#### for systems with core shell particles in general; nbforce.getNumParticles() is cores+shells  #### TODO: get the real number of particle
Numb_particles = nbforce.getNumParticles()

group_range = range(1, len(force_group_names.keys()))
for group in group_range:
    print('%8d : %64s : %16.4f kJ/mol' % (group, force_group_names[group], simulation.context.getState(getEnergy=True, groups=(1 << group)).getPotentialEnergy()/unit.kilojoule_per_mole))
potE = simulation.context.getState(getEnergy=True).getPotentialEnergy()/unit.kilojoule_per_mole 
print('potential energy = {0:.2f} kJ/mol\n'.format(potE/Numb_particles))    

#### Minimize and Equilibrate ####
print('Performing energy minimization...')
simulation.minimizeEnergy()
potE = simulation.context.getState(getEnergy=True).getPotentialEnergy()/unit.kilojoule_per_mole 
print('potential energy = {0:.2f} kJ/mol\n'.format(potE/Numb_particles)) 

print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature_c)
simulation.step(equilibrationSteps)

# Simulate
for group in group_range:
    print('%8d : %64s : %16.4f kJ/mol' % (group, force_group_names[group], simulation.context.getState(getEnergy=True, groups=(1 << group)).getPotentialEnergy()/unit.kilojoule_per_mole))
potE = simulation.context.getState(getEnergy=True).getPotentialEnergy()/unit.kilojoule_per_mole 
print('potential energy = {0:.2f} kJ/mol\n'.format(potE/Numb_particles))

print('Simulating...')
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(pdbReporter)
simulation.reporters.append(chkReporter)
simulation.currentStep = 0
simulation.step(steps)


for group in group_range:
    print('%8d : %64s : %16.4f kJ/mol' % (group, force_group_names[group], simulation.context.getState(getEnergy=True, groups=(1 << group)).getPotentialEnergy()/unit.kilojoule_per_mole))
potE = simulation.context.getState(getEnergy=True).getPotentialEnergy()/unit.kilojoule_per_mole 
print('potential energy = {0:.2f} kJ/mol\n'.format(potE/Numb_particles)) 

del simulation.context, integrator
