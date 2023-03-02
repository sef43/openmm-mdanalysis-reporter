# OpenMM MDAnalysis Reporter

*Currently a WIP but it should work* 

This is an [OpenMM](https://openmm.org/) Reporter class that uses [MDAnalysis](https://www.mdanalysis.org/) for output. This means it can use any [format supported by MDAnalysis](https://userguide.mdanalysis.org/stable/formats/index.html).
It also supports MDAnalysis [selection strings](https://docs.mdanalysis.org/stable/documentation_pages/selections.html).

The reporter is called `MDAReporter`, once this package is installed it can be imported as:
```python
from mdareporter import MDAReporter
```

# Installation
Pip:
```bash
pip install git+https://github.com/sef43/openmm-mdanalysis-reporter
```
From source:
```bash
git clone https://github.com/sef43/openmm-mdanalysis-reporter.git
cd openmm-mdanalysis-reporter
pip install .
```


# Usage
It is used like the existing OpenMM reporters, where the file format is read from the suffix of the output file by MDAnalysis, e.g. to output in XYZ format:
```python

from mdareporter import MDAReporter

...

simulation.reporters.append(MDAReporter('traj.xyz',100))
...

```

Additionally it supports the [MDAnalysis selection syntax](https://docs.mdanalysis.org/stable/documentation_pages/selections.html)
```python

# using MDAnalysis selection string to output just Carbon atoms

simulation.reporters.append(MDAReporter('traj.xyz',100, selection='name is C'))
```

# Example

Full example of using MDAReporter to output just the protein (script and data file in `cd` [example](./example)):
```python
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
    
```