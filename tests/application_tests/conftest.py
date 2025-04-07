"""
Pytest configuration for application tests.
"""
import os
import sys
import tempfile
import pytest
import base64
from pathlib import Path
from flask import Flask
import requests
from unittest.mock import patch

# Add the project root to the Python path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from pymol_fitter_server.app import app as flask_app


@pytest.fixture(scope="session")
def app():
    """Create a Flask app for testing."""
    # Configure app for testing
    flask_app.config.update({
        "TESTING": True,
    })
    yield flask_app


@pytest.fixture(scope="session")
def client(app):
    """Create a test client for the Flask app."""
    return app.test_client()


@pytest.fixture(scope="session")
def docker_endpoint():
    """Return the Docker endpoint for testing.
    
    This is the URL for the Docker container running the Flask API.
    """
    # Default to localhost:5000, but could be overridden via environment variable
    return os.environ.get("PYMOL_FITTER_DOCKER_URL", "http://localhost:5000")


@pytest.fixture(scope="session")
def example_files(example_dir):
    """Return paths to sample files for testing."""
    return {
        "protein": example_dir / "8gcy.pdb",
        "ligand_sdf": example_dir / "Ligand.sdf",
        "crystal_sdf": example_dir / "Crystal.sdf",
    }


@pytest.fixture(scope="function")
def encoded_test_files(example_files):
    """Return base64 encoded content of example files."""
    encoded = {}
    for key, path in example_files.items():
        if path.exists():
            with open(path, "rb") as f:
                encoded[key] = base64.b64encode(f.read()).decode("utf-8")
        else:
            pytest.skip(f"Required example file not found: {path}")
    return encoded


@pytest.fixture(scope="function")
def mock_docking_engine():
    """Mock the Docking_Engine module for faster tests."""
    with patch("pymol_fitter_server.pymol_fitter_src.Docking_Engine.Pymol_Docking") as mock_docking:
        # Configure mock to return predictable values
        instance = mock_docking.return_value
        
        # Mock run_smina_docking to create some dummy output files
        def mock_run_docking(mode, output_name):
            # Create temporary output files
            temp_dir = tempfile.gettempdir()
            smina_output = Path(temp_dir) / f"{output_name}.sdf"
            protein_prep = Path(temp_dir) / f"{output_name}_protein.pdb"
            log_path = Path(temp_dir) / f"{output_name}.log"
            
            # Write some dummy content
            with open(smina_output, "w") as f:
                f.write("MOCK DOCKING RESULT")
            with open(protein_prep, "w") as f:
                f.write("MOCK PROTEIN PREPARATION")
            with open(log_path, "w") as f:
                f.write("MOCK DOCKING LOG")
            
            return smina_output, protein_prep, log_path
        
        instance.run_smina_docking.side_effect = mock_run_docking
        
        # Mock run_complex_minimization
        def mock_minimization(protein_path, ligand_path=None):
            temp_dir = tempfile.gettempdir()
            complex_path = Path(temp_dir) / "minimized_complex.pdb"
            with open(complex_path, "w") as f:
                f.write("MOCK MINIMIZED COMPLEX")
            return complex_path
        
        instance.run_complex_minimization.side_effect = mock_minimization
        
        yield mock_docking


@pytest.fixture
def docker_health_check(docker_endpoint):
    """Skip tests if Docker container is not running."""
    try:
        response = requests.get(f"{docker_endpoint}/health", timeout=2)
        if response.status_code != 200:
            pytest.skip("Docker container not healthy")
    except requests.RequestException:
        pytest.skip("Docker container not accessible")
