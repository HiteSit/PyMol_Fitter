"""
Tests for the Flask API endpoints in pymol_fitter_server.app.

These tests focus on testing the API endpoints and their integration with
the scientific components.
"""
import json
import pytest
from pathlib import Path

from pymol_fitter_server.app import app


class TestFlaskAPI:
    """Test suite for the Flask API endpoints."""
    
    def test_health_check(self, client):
        """Test the health check endpoint."""
        response = client.get("/health")
        assert response.status_code == 200
        assert response.json["status"] == "healthy"
    
    def test_dock_empty_request(self, client):
        """Test the dock endpoint with empty request body."""
        response = client.post("/dock", json={})
        assert response.status_code == 400
        assert "success" in response.json
        assert response.json["success"] is False
        assert "message" in response.json
        assert "No protein data provided" in response.json["message"]
    
    def test_dock_missing_ligand(self, client, encoded_test_files):
        """Test the dock endpoint with missing ligand data."""
        data = {
            "protein": encoded_test_files["protein"],
            "crystal": encoded_test_files["crystal_sdf"],
            "is_smiles": False,
            "minimize": False,
            "dock_mode": "Dock",
            "output_name": "test_dock"
        }
        response = client.post("/dock", json=data)
        assert response.status_code == 400
        assert "success" in response.json
        assert response.json["success"] is False
        assert "message" in response.json
        assert "No ligand data provided" in response.json["message"]
    
    def test_dock_missing_crystal(self, client, encoded_test_files):
        """Test the dock endpoint with missing crystal data."""
        data = {
            "protein": encoded_test_files["protein"],
            "ligand": encoded_test_files["ligand_sdf"],
            "is_smiles": False,
            "minimize": False,
            "dock_mode": "Dock",
            "output_name": "test_dock"
        }
        response = client.post("/dock", json=data)
        assert response.status_code == 400
        assert "success" in response.json
        assert response.json["success"] is False
        assert "message" in response.json
        assert "No crystal data provided" in response.json["message"]
    
    @pytest.mark.parametrize("minimize", [True, False])
    def test_dock_valid_request_mock(self, client, encoded_test_files, mock_docking_engine, minimize):
        """Test the dock endpoint with a valid request using mocked docking engine."""
        data = {
            "protein": encoded_test_files["protein"],
            "ligand": encoded_test_files["ligand_sdf"],
            "crystal": encoded_test_files["crystal_sdf"],
            "is_smiles": False,
            "minimize": minimize,
            "dock_mode": "Dock",
            "output_name": "test_dock"
        }
        response = client.post("/dock", json=data)
        assert response.status_code == 200
        assert response.json["success"] is True
        assert "message" in response.json
        assert "Docking completed successfully" in response.json["message"]
        assert "log" in response.json
        assert "docked_ligand" in response.json
        assert "prepared_protein" in response.json
        
        if minimize:
            assert "minimized_complex" in response.json
    
    def test_minimize_empty_request(self, client):
        """Test the minimize endpoint with empty request body."""
        response = client.post("/minimize", json={})
        assert response.status_code == 400
        assert "success" in response.json
        assert response.json["success"] is False
        assert "message" in response.json
        assert "No protein data provided" in response.json["message"]
    
    def test_minimize_valid_request_mock(self, client, encoded_test_files, mock_docking_engine):
        """Test the minimize endpoint with a valid request using mocked minimization."""
        with app.test_client() as client:
            data = {
                "protein": encoded_test_files["protein"],
                "ligand": encoded_test_files["ligand_sdf"],
                "output_name": "test_minimization"
            }
            
            # Use the mocked function
            response = client.post("/minimize", json=data)
            
            assert response.status_code == 200
            assert response.json["success"] is True
            assert "message" in response.json
            assert "Minimization completed successfully" in response.json["message"]
            assert "minimized_complex" in response.json
    
    def test_minimize_protein_only_mock(self, client, encoded_test_files, mock_docking_engine):
        """Test the minimize endpoint with protein-only minimization using mocked minimization."""
        with app.test_client() as client:
            data = {
                "protein": encoded_test_files["protein"],
                "output_name": "test_minimization_protein_only"
            }
            
            # Use the mocked function
            response = client.post("/minimize", json=data)
            
            assert response.status_code == 200
            assert response.json["success"] is True
            assert "message" in response.json
            assert "Minimization completed successfully" in response.json["message"]
            assert "minimized_complex" in response.json
