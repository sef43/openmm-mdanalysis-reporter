# OpenMM MDAnalysis Reporter

[![CI](https://github.com/sef43/openmm-mdanalysis-reporter/actions/workflows/CI.yml/badge.svg)](https://github.com/sef43/openmm-mdanalysis-reporter/actions/workflows/CI.yml)

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

simulation.reporters.append(MDAReporter('traj.xyz',100, selection='name C'))
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

# Testing
testsuite can be run using `pytest`
```bash
cd mdareporter/tests
pytest
```

# Benchmarks
A benchmarking script which writes trajectory snapshots every 10 steps for 1000 steps can be found at [example/benchmark.py](example/benchmark.py)

The output on a M2 Macbook is:
```
OpenMM PDBReporter time =  7.989748208987294 s
OpenMM DCDReporter time =  5.6639587499958 s
MDAReporter format DCD time =  3.5682871250028256 s
MDAReporter format NCDF time =  3.609358207992045 s
MDAReporter format PDB time =  11.491382707987214 s
MDAReporter format TRR time =  4.7894440419913735 s
MDAReporter format XTC time =  4.086603666975861 s
MDAReporter format XYZ time =  5.835725833981996 s
```

Excluding MDAReporter using PDB format the MDAReporter formats are all as fast, or faster than the OpenMM formats.
