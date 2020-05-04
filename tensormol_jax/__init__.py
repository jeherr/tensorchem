"""
TensorMol-Jax
A reimplementation of TensorMol in Jax
"""

# Add imports here
from .tensormol_jax import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
