openmm-mdanalysis-reporter
==============================
[//]: # (Badges)

| **Latest release** | [![Last release tag](https://img.shields.io/github/release-pre/sef43/openmm-mdanalysis-reporter.svg)](https://github.com/sef43/openmm-mdanalysis-reporter/releases) ![GitHub commits since latest release (by date) for a branch](https://img.shields.io/github/commits-since/sef43/openmm-mdanalysis-reporter/latest)  |
| :------ | :------- |
| **Status** | [![GH Actions Status](https://github.com/sef43/openmm-mdanalysis-reporter/actions/workflows/gh-ci.yaml/badge.svg)](https://github.com/sef43/openmm-mdanalysis-reporter/actions?query=branch%3Amain+workflow%3Agh-ci) [![codecov](https://codecov.io/gh/sef43/openmm-mdanalysis-reporter/branch/main/graph/badge.svg)](https://codecov.io/gh/sef43/openmm-mdanalysis-reporter/branch/main) |
| **Community** | [![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)|

MDAnalysis based reporter for OpenMM

openmm-mdanalysis-reporter is bound by a [Code of Conduct](https://github.com/sef43/openmm-mdanalysis-reporter/blob/main/CODE_OF_CONDUCT.md).

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
```
cd mdareporter/tests
pytest
```
or
```
pytest --pyargs mdareporter.tests   
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


#### Acknowledgements
 
Project based on the 
[MDAnalysis Cookiecutter](https://github.com/MDAnalysis/cookiecutter-mda) version 0.1.
Please cite [MDAnalysis](https://github.com/MDAnalysis/mdanalysis#citation) when using openmm-mdanalysis-reporter in published work.
