"""
Tests for the PyMOL Docking Client in pymol_fitter_plugin.client.

These tests focus on testing the client's functionality and its interaction
with the server API.
"""
import os
import json
import pytest
import tempfile
import requests
from pathlib import Path
from unittest.mock import patch, MagicMock

from pymol_fitter_plugin.client import PyMOLDockingClient


class TestPyMOLDockingClient:
    """Test suite for the PyMOL Docking Client."""
    
    @pytest.fixture(scope="class")
    def base_url(self):
        """Return the base URL for testing."""
        return "http://localhost:5000"
    
    @pytest.fixture(scope="function")
    def client(self, base_url):
        """Create a PyMOLDockingClient for testing."""
        return PyMOLDockingClient(base_url=base_url)
    
    @pytest.fixture(scope="class")
    def example_files(self, example_dir):
        """Return paths to sample files for testing."""
        return {
            "protein": example_dir / "8gcy.pdb",
            "ligand_sdf": example_dir / "Ligand.sdf",
            "crystal_sdf": example_dir / "Crystal.sdf",
        }
    
    @pytest.fixture(scope="function")
    def temp_output_dir(self):
        """Create a temporary directory for output files."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            yield Path(tmp_dir)
    
    def test_initialization(self, base_url):
        """Test that the client initializes correctly."""
        client = PyMOLDockingClient(base_url=base_url)
        assert client.base_url == base_url
    
    def test_encode_file(self, client, example_files):
        """Test the _encode_file method."""
        encoded = client._encode_file(example_files["protein"])
        assert isinstance(encoded, str)
        assert len(encoded) > 0
    
    def test_encode_file_nonexistent(self, client):
        """Test _encode_file with non-existent file."""
        with pytest.raises(FileNotFoundError):
            client._encode_file(Path("nonexistent.pdb"))
    
    def test_decode_and_save(self, client, temp_output_dir):
        """Test the _decode_and_save method."""
        # Create a simple test string
        test_content = b"TEST DATA"
        encoded = "VEVTVCBEQVRB"  # Base64 for "TEST DATA"
        
        output_path = temp_output_dir / "decoded.txt"
        result_path = client._decode_and_save(encoded, output_path)
        
        assert result_path.exists()
        assert result_path == output_path
        with open(result_path, "rb") as f:
            assert f.read() == test_content
    
    def test_check_health_success(self, client):
        """Test the check_health method when server is healthy."""
        with patch("requests.get") as mock_get:
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.json.return_value = {"status": "healthy"}
            mock_get.return_value = mock_response
            
            result = client.check_health()
            assert result is True
            mock_get.assert_called_once_with(f"{client.base_url}/health", timeout=5)
    
    def test_check_health_failure(self, client):
        """Test the check_health method when server is unhealthy."""
        with patch("requests.get") as mock_get:
            # Test server error
            mock_response = MagicMock()
            mock_response.status_code = 500
            mock_get.return_value = mock_response
            
            result = client.check_health()
            assert result is False
            
            # Test connection error
            mock_get.side_effect = requests.RequestException()
            result = client.check_health()
            assert result is False
    
    def test_dock_minimize_validation(self, client, example_files):
        """Test input validation in dock_minimize method."""
        # Test with invalid dock mode
        with pytest.raises(ValueError, match="dock_mode must be either 'Dock' or 'Minimize'"):
            with patch("pymol_fitter_plugin.client.PyMOLDockingClient.check_health", return_value=True):
                client.dock_minimize(
                    protein_file=example_files["protein"],
                    ligand=example_files["ligand_sdf"],
                    crystal_file=example_files["crystal_sdf"],
                    dock_mode="InvalidMode"
                )
    
    def test_dock_minimize_success(self, client, example_files, temp_output_dir):
        """Test successful docking and minimization."""
        # Mock the API response
        with patch("pymol_fitter_plugin.client.PyMOLDockingClient.check_health", return_value=True), \
             patch("requests.post") as mock_post:
            
            # Create mock response
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.json.return_value = {
                "success": True,
                "message": "Docking completed successfully with mode: Dock",
                "log": "TU9DSyBMT0c=",  # Base64 for "MOCK LOG"
                "docked_ligand": "TU9DSyBMSUdBTkQ=",  # Base64 for "MOCK LIGAND"
                "prepared_protein": "TU9DSyBQUk9URUlO",  # Base64 for "MOCK PROTEIN"
                "minimized_complex": "TU9DSyBDT01QTEVYIA=="  # Base64 for "MOCK COMPLEX "
            }
            mock_post.return_value = mock_response
            
            # Call the method
            result = client.dock_minimize(
                protein_file=example_files["protein"],
                ligand=example_files["ligand_sdf"],
                crystal_file=example_files["crystal_sdf"],
                minimize=True,
                dock_mode="Dock",
                output_name="test_dock",
                output_dir=temp_output_dir
            )
            
            # Check results
            assert isinstance(result, dict)
            assert "log" in result
            assert "docked_ligand" in result
            assert "prepared_protein" in result
            assert "minimized_complex" in result
            
            # Check that files were created
            assert result["log"].exists()
            assert result["docked_ligand"].exists()
            assert result["prepared_protein"].exists()
            assert result["minimized_complex"].exists()
            
            # Check file content
            with open(result["log"], "r") as f:
                assert f.read() == "MOCK LOG"
    
    def test_dock_minimize_server_error(self, client, example_files):
        """Test handling of server errors."""
        with patch("pymol_fitter_plugin.client.PyMOLDockingClient.check_health", return_value=True), \
             patch("requests.post") as mock_post:
            
            # Create mock response for server error
            mock_response = MagicMock()
            mock_response.status_code = 500
            mock_response.json.return_value = {
                "success": False,
                "message": "Server error"
            }
            mock_post.return_value = mock_response
            
            # Check that the error is properly raised
            with pytest.raises(RuntimeError, match="API call failed: Server error"):
                client.dock_minimize(
                    protein_file=example_files["protein"],
                    ligand=example_files["ligand_sdf"],
                    crystal_file=example_files["crystal_sdf"],
                    dock_mode="Dock"
                )
    
    def test_dock_minimize_connection_error(self, client, example_files):
        """Test handling of connection errors."""
        with patch("pymol_fitter_plugin.client.PyMOLDockingClient.check_health", return_value=True), \
             patch("requests.post") as mock_post:
            
            # Simulate connection error
            mock_post.side_effect = requests.RequestException("Connection error")
            
            # Check that the error is properly raised
            with pytest.raises(ConnectionError):
                client.dock_minimize(
                    protein_file=example_files["protein"],
                    ligand=example_files["ligand_sdf"],
                    crystal_file=example_files["crystal_sdf"],
                    dock_mode="Dock"
                )
    
    def test_inplace_minimization_success(self, client, example_files, temp_output_dir):
        """Test successful in-place minimization."""
        # Mock the API response
        with patch("pymol_fitter_plugin.client.PyMOLDockingClient.check_health", return_value=True), \
             patch("requests.post") as mock_post:
            
            # Create mock response
            mock_response = MagicMock()
            mock_response.status_code = 200
            mock_response.json.return_value = {
                "success": True,
                "message": "Minimization completed successfully",
                "minimized_complex": "TU9DSyBDT01QTEVYIA=="  # Base64 for "MOCK COMPLEX "
            }
            mock_post.return_value = mock_response
            
            # Call the method
            result = client.inplace_minimization(
                protein_file=example_files["protein"],
                ligand=example_files["ligand_sdf"],
                output_name="test_minimization",
                output_dir=temp_output_dir
            )
            
            # Check results
            assert isinstance(result, dict)
            assert "minimized_complex" in result
            
            # Check that files were created
            assert result["minimized_complex"].exists()
            
            # Check file content
            with open(result["minimized_complex"], "rb") as f:
                assert f.read() == b"MOCK COMPLEX "
