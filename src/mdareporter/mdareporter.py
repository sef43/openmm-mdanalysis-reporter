import MDAnalysis as mda
from openmm import unit

class MDAReporter(object):
    """MDAReporter outputs a series of frames from a Simulation to any file format supported by MDAnalysis.
    To use it, create a MDAReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, enforcePeriodicBox=None, selection:str=None):
        """Create a MDAReporter.
        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        selection : str
            MDAnalysis selection string (https://docs.mdanalysis.org/stable/documentation_pages/selections.html)
            which will be passed to MDAnalysis.Universe.select_atoms. If None (the default), all atoms will we selected.
            
        """
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._filename = file
        self._topology = None
        self._nextModel = 0
        self._mdaUniverse = None
        self._mdaWriter = None
        self._selection = selection
        self._atomGroup = None

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.
        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        Returns
        -------
        tuple
            A six element tuple. The first element is the number of steps
            until the next report. The next four elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.  The final element specifies whether
            positions should be wrapped to lie in a single periodic box.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, True, False, False, False, self._enforcePeriodicBox)

    def report(self, simulation, state):
        """Generate a report.
        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if self._nextModel == 0:
            self._topology = simulation.topology
            self._mdaUniverse = mda.Universe(simulation.topology, simulation,topology_format='OPENMMTOPOLOGY',format='OPENMMSIMULATION', dt=simulation.currentStep)
            if self._selection is not None:
                self._atomGroup = self._mdaUniverse.select_atoms(self._selection)
            else:
                self._atomGroup = self._mdaUniverse.atoms
            print(self._atomGroup)
            self._mdaWriter = mda.Writer(self._filename, n_atoms=len(self._atomGroup))
            self._nextModel += 1

        # update the positions, convert from OpenMM nm to MDAnalysis angstroms
        self._mdaUniverse.atoms.positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

        # write to the trajectory file
        self._mdaWriter.write(self._atomGroup)

        self._nextModel += 1

    def __del__(self):
        if self._mdaWriter:
            self._mdaWriter.close()
