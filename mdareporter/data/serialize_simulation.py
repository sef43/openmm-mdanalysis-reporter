from openmm.app import ForceField, Simulation, PME, HBonds
from openmm import LangevinMiddleIntegrator
from openmm.unit import nanometer, picosecond, kelvin

from openmm.app import PDBFile, Topology
from openmm.openmm import XmlSerializer
from mdareporter.data.files import VILLIN_PDB


pdb = PDBFile(VILLIN_PDB)

forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
system = forcefield.createSystem(
    pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1 * nanometer, constraints=HBonds
)
integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picosecond)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()


pdb.writeFile(
    simulation.topology,
    simulation.context.getState(getPositions=True).getPositions(),
    open("villin_simulation/output.pdb", "w"),
)
system_xml = XmlSerializer.serialize(system)
integrator_xml = XmlSerializer.serialize(integrator)
state_xml = XmlSerializer.serialize(
    simulation.context.getState(getPositions=True, getVelocities=True)
)

# write to files
with open("villin_simulation/system.xml", "w") as f:
    f.write(system_xml)
with open("villin_simulation/integrator.xml", "w") as f:
    f.write(integrator_xml)
with open("villin_simulation/state.xml", "w") as f:
    f.write(state_xml)
