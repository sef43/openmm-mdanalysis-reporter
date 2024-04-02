"""
Location of data files
======================

Use as ::

    from mdareporter.data.files import *

"""

__all__ = ["MDANALYSIS_LOGO", "VILLIN_PDB"]  # example file of MDAnalysis logo

from pkg_resources import resource_filename

MDANALYSIS_LOGO = resource_filename(__name__, "mda.txt")
VILLIN_PDB = resource_filename(__name__, "villin.pdb")
