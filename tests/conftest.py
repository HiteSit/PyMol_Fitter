"""
Global pytest configuration file for PyMol Fitter tests.
"""
import os
import sys
import tempfile
import pytest
from pathlib import Path

# Add the project root to the Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Test data directory
@pytest.fixture(scope="session")
def data_dir() -> Path:
    """Return the path to the test data directory."""
    return Path(__file__).parent / "data"

@pytest.fixture(scope="session")
def example_dir() -> Path:
    """Return the path to the examples directory."""
    return Path(__file__).parent.parent / "examples"

@pytest.fixture(scope="function")
def temp_dir():
    """Create a temporary directory for tests that need file I/O."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)
