from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from mdareporter import MDAReporter
import MDAnalysis as mda
import numpy as np
import pytest
import tempfile


from mdareporter.data.files import VILLIN_PDB

@pytest.mark.parametrize("file_ext", ["DCD", "NCDF","PDB", "TRR", "XTC", "XYZ"])
def test_mdareporter(file_ext):

    with tempfile.TemporaryDirectory() as tempdir:
        traj_name = os.path.join(tempdir, 'test_traj.'+file_ext)
        

        # run a OpenMM simulation
        pdb = PDBFile(os.path.join(VILLIN_PDB))
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(maxIterations=10)

        # output a traj with M frames every N steps
        M=10
        N=10
        valid_positions = []
        simulation.reporters.append(MDAReporter(traj_name, N, enforcePeriodicBox=False))
        simulation.reporters.append(StateDataReporter(stdout, N, step=True,
                potentialEnergy=True, temperature=True))

        for i in range(M):
            simulation.step(N)
            positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(angstrom)
            valid_positions.append(positions)

        simulation.reporters.clear()

        # read in the traj with MDA
        # and compare the positions with the positions extracted direct from OpenMM
        u = mda.Universe(traj_name)

        test_positions = []
        for ts in u.trajectory:
            coords = u.atoms.positions
            test_positions.append(coords)

        assert(len(valid_positions) == len(test_positions))

        max_errors=[]
        for valid, test in zip(valid_positions, test_positions):
            # units are all in Angstrom
            max_errors.append(np.max(np.abs(valid-test)))

            assert(np.allclose(test, valid, atol=0.01))
        
        print("max error in positions is ",np.max(max_errors)," Angstrom")

    
@pytest.mark.parametrize("file_ext", ["DCD", "NCDF","PDB", "TRR", "XTC", "XYZ"])
def test_mdareporter_selection(file_ext):

    with tempfile.TemporaryDirectory() as tempdir:
        traj_name = os.path.join(tempdir, 'test_traj_selection.'+file_ext)

        # run a OpenMM simulation
        pdb = PDBFile(VILLIN_PDB)
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy(maxIterations=10)

        # output a traj with M frames every N steps of just the protein
        M=10
        N=10
        protein_indices = list(range(582))
        valid_positions = []
        simulation.reporters.append(MDAReporter(traj_name, N, enforcePeriodicBox=False, selection="protein"))
        simulation.reporters.append(StateDataReporter(stdout, N, step=True,
                potentialEnergy=True, temperature=True))

        for i in range(M):
            simulation.step(N)
            positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(angstrom)
            valid_positions.append(positions[protein_indices])

        simulation.reporters.clear()

        # read in the traj with MDA
        # and compare the positions with the positions extracted direct from OpenMM
        u = mda.Universe(traj_name)

        test_positions = []
        for ts in u.trajectory:
            coords = u.atoms.positions
            test_positions.append(coords)

        assert(len(valid_positions) == len(test_positions))

        max_errors=[]
        for valid, test in zip(valid_positions, test_positions):

            assert(valid.shape == test.shape)
            # units are all in Angstrom
            max_errors.append(np.max(np.abs(valid-test)))

            assert(np.allclose(test, valid, atol=0.01))
        
        print("max error in positions is ",np.max(max_errors)," Angstrom")


@pytest.mark.parametrize("file_ext", ["DCD", "NCDF","PDB", "TRR", "XTC", "XYZ"])
def test_mdareporter_boxvectors(file_ext):

    with tempfile.TemporaryDirectory() as tempdir:
        traj_name = os.path.join(tempdir, 'test_traj.'+file_ext)
        

        # run a OpenMM simulation
        pdb = PDBFile(os.path.join(VILLIN_PDB))
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                nonbondedCutoff=1*nanometer, constraints=HBonds)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)

        simulation.context.setPeriodicBoxVectors((5.0,0.0,0.0), (0.0,5.0,0.0), (0.0,1.0,5.0))

        simulation.minimizeEnergy(maxIterations=10)

        # output a traj with M frames every N steps
        M=10
        N=10
        valid_positions = []
        simulation.reporters.append(MDAReporter(traj_name, N, enforcePeriodicBox=True))
        simulation.reporters.append(StateDataReporter(stdout, N, step=True,
                potentialEnergy=True, temperature=True))

        for i in range(M):
            simulation.step(N)
            positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(asNumpy=True).value_in_unit(angstrom)
            valid_positions.append(positions)

        simulation.reporters.clear()

        # read in the traj with MDA
        # and compare the positions with the positions extracted direct from OpenMM
        u = mda.Universe(traj_name)

        # not all formats output boxvectors
        if u.dimensions is not None:
            valid_dimensions = mda.lib.mdamath.triclinic_box(*simulation.context.getState().getPeriodicBoxVectors(asNumpy=True).value_in_unit(angstrom))
            assert(np.allclose(u.dimensions, valid_dimensions))

        test_positions = []
        for ts in u.trajectory:
            coords = u.atoms.positions
            test_positions.append(coords)

        assert(len(valid_positions) == len(test_positions))

        max_errors=[]
        for valid, test in zip(valid_positions, test_positions):
            # units are all in Angstrom
            max_errors.append(np.max(np.abs(valid-test)))

            assert(np.allclose(test, valid, atol=0.01))
        
        print("max error in positions is ",np.max(max_errors)," Angstrom")

    