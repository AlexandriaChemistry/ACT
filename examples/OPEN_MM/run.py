#!/usr/bin/env python3

from act_openmm import ActOpenMMSim

sim = ActOpenMMSim(pdbfile="MX_1000.pdb", datfile="simulation.dat", xmlfile="openmm_ff.xml")
sim.run()
sim.log_to_xvg("energy.xvg", [ "Potential Energy (kJ/mole)"])
sim.log_to_xvg("temperature.xvg", [ "Temperature" ])
sim.log_to_xvg("density.xvg", [ "Density" ])
