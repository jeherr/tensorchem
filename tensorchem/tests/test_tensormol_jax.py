"""
Unit and regression test for the tensorchem package.
"""

# Import package, test suite, and other packages as needed
import tensorchem
import pytest
import sys

def test_tensorchem_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "tensorchem" in sys.modules
