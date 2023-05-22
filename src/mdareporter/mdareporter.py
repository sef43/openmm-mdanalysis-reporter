"""
mdareporter.py Outputs OpenMM simulation trajectories using 
MDAnalysis in any format MDAnalysis can write
Authors: Stephen Farr

OpenMM: https://openmm.org/
MDAnalysis: https://www.mdanalysis.org/

Modified from pdbreporter.py in OpenMM

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""



import MDAnalysis as mda
from MDAnalysis.lib.mdamath import triclinic_box
from openmm import unit
import numpy as np

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
            dt = simulation.integrator.getStepSize()*self._reportInterval
            self._mdaUniverse = mda.Universe(simulation.topology, simulation,topology_format='OPENMMTOPOLOGY',format='OPENMMSIMULATION', dt=dt)
            if self._selection is not None:
                self._atomGroup = self._mdaUniverse.select_atoms(self._selection)
            else:
                self._atomGroup = self._mdaUniverse.atoms
            self._mdaWriter = mda.Writer(self._filename, n_atoms=len(self._atomGroup))
            self._nextModel += 1

        # update the positions, convert from OpenMM nm to MDAnalysis angstroms
        self._mdaUniverse.atoms.positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

        # update box vectors
        boxVectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
        self._mdaUniverse.dimensions = triclinic_box(*boxVectors)
        

        # write to the trajectory file
        self._mdaWriter.write(self._atomGroup)

        self._nextModel += 1

    def __del__(self):
        if self._mdaWriter:
            self._mdaWriter.close()


def _sanitize_box_angles(angles):
    """ Ensure box angles correspond to first quadrant

    See `discussion on unitcell angles <https://github.com/MDAnalysis/mdanalysis/pull/2917/files#r620558575>`_
    """
    inverted = 180 - angles

    return np.min(np.array([angles, inverted]), axis=0)