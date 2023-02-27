# Openmm MDAnalysis Reporter

This is an [OpenMM](https://openmm.org/) Reporter class that uses [MDAnalysis](https://www.mdanalysis.org/) for output. This means it can use any format supported by MDAnalysis.
It also supports MDAnalysis selection strings.

The reporter is called 'MDAReporter', oncee this package is installed it can be imported as:
```python
from mdareporter import MDAReporter
```

# Installation
Pip:
```bash
pip install git+
```
From source:
```bash
git clone 
cd 
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