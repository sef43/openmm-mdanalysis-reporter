"""
openmm-mdanalysis-reporter
MDAnalysis based reporter for OpenMM
"""

# # Add imports here
from .mdareporter import MDAReporter

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
