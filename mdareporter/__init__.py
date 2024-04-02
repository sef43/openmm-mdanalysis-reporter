"""
openmm-mdanalysis-reporter
MDAnalysis based reporter for OpenMM
"""

# # Add imports here
from .mdareporter import MDAReporter

# Handle version
from importlib.metadata import version

__version__ = version("mdareporter")
