"""
Unit and regression test for the tensormol_jax package.
"""

# Import package, test suite, and other packages as needed
import tensormol_jax
import pytest
import sys

def test_tensormol_jax_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "tensormol_jax" in sys.modules
