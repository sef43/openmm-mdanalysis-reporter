from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from mdareporter import MDAReporter


pdb = PDBFile('villin.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy(maxIterations=100)

# output just protein in xyz format
simulation.reporters.append(MDAReporter('traj.xyz',100, enforcePeriodicBox=False, selection="protein"))

simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
        potentialEnergy=True, temperature=True))

simulation.step(1000)
    
