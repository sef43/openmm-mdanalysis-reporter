from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from mdareporter import MDAReporter
from time import perf_counter


def run_openmm_pdbreporter():
    # pure OpenMM
    pdb = PDBFile('villin.pdb')
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.reporters.append(PDBReporter('traj.pdb',10, enforcePeriodicBox=False))
    simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
            potentialEnergy=True, temperature=True))

    simulation.step(1000)

def run_openmm_dcdreporter():
    # pure OpenMM
    pdb = PDBFile('villin.pdb')
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.reporters.append(DCDReporter('traj.pdb',10, enforcePeriodicBox=False))
    simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
            potentialEnergy=True, temperature=True))

    simulation.step(1000)

def run_mdareporter(ext):
    # MDAReporter
    pdb = PDBFile('villin.pdb')
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    simulation.reporters.append(MDAReporter('traj.'+ext,10, enforcePeriodicBox=False))

    simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
            potentialEnergy=True, temperature=True))

    simulation.step(1000)


if __name__ == "__main__":

    t1 = perf_counter()
    run_openmm_pdbreporter()
    t2 = perf_counter()
    print("OpenMM PDBReporter time = ", t2-t1, "s")

    t1 = perf_counter()
    run_openmm_dcdreporter()
    t2 = perf_counter()
    print("OpenMM DCDReporter time = ", t2-t1, "s")

    for ext in ["DCD", "NCDF","PDB", "TRR", "XTC", "XYZ"]:
        t1 = perf_counter()
        run_mdareporter(ext)
        t2 = perf_counter()
        print("MDAReporter format",ext,"time = ", t2-t1, "s")
