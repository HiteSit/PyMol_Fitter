"""
Tests for the CDPK_Utils module in pymol_fitter_src.

These tests focus on the scientific functionality of the CDPK utilities,
independent of the Flask API layer.
"""
import os
import pytest
from pathlib import Path

from pymol_fitter_server.pymol_fitter_src.CDPK_Utils import CDPK_Runner


class TestCDPKRunner:
    """Test suite for the CDPK_Runner class."""
    
    @pytest.fixture(scope="class")
    def example_files(self, example_dir):
        """Return paths to sample files from the examples directory."""
        return {
            "ligand": example_dir / "Ligand.sdf"
        }
    
    @pytest.fixture(scope="function")
    def output_path(self, temp_dir):
        """Return a temporary output path."""
        return temp_dir / "output_ligand.sdf"
    
    def test_initialization_default(self):
        """Test CDPK_Runner initialization with default parameters."""
        cdpk = CDPK_Runner()
        assert cdpk is not None
        assert cdpk.standardize is True
        assert cdpk.protonate is True
        assert cdpk.gen3d is True
    
    def test_initialization_custom(self):
        """Test CDPK_Runner initialization with custom parameters."""
        cdpk = CDPK_Runner(standardize=False, protonate=True, gen3d=False)
        assert cdpk is not None
        assert cdpk.standardize is False
        assert cdpk.protonate is True
        assert cdpk.gen3d is False
    
    @pytest.mark.slow
    def test_prepare_ligands(self, example_files, output_path):
        """Test ligand preparation with CDPK (slow)."""
        try:
            cdpk = CDPK_Runner()
            cdpk.prepare_ligands(example_files["ligand"], output_path)
            
            assert output_path.exists()
            assert output_path.stat().st_size > 0
        except (ImportError, ModuleNotFoundError) as e:
            pytest.skip(f"Required CDPL module not found: {e}. Skipping test.")
    
    @pytest.mark.slow
    def test_prepare_ligands_no_3d(self, example_files, output_path):
        """Test ligand preparation without 3D generation (slow)."""
        try:
            cdpk = CDPK_Runner(gen3d=False)
            cdpk.prepare_ligands(example_files["ligand"], output_path)
            
            assert output_path.exists()
            assert output_path.stat().st_size > 0
        except (ImportError, ModuleNotFoundError) as e:
            pytest.skip(f"Required CDPL module not found: {e}. Skipping test.")
    
    @pytest.mark.slow
    def test_prepare_ligands_no_standardize(self, example_files, output_path):
        """Test ligand preparation without standardization (slow)."""
        try:
            cdpk = CDPK_Runner(standardize=False)
            cdpk.prepare_ligands(example_files["ligand"], output_path)
            
            assert output_path.exists()
            assert output_path.stat().st_size > 0
        except (ImportError, ModuleNotFoundError) as e:
            pytest.skip(f"Required CDPL module not found: {e}. Skipping test.")
    
    @pytest.mark.slow
    def test_prepare_ligands_no_protonate(self, example_files, output_path):
        """Test ligand preparation without protonation (slow)."""
        try:
            cdpk = CDPK_Runner(protonate=False)
            cdpk.prepare_ligands(example_files["ligand"], output_path)
            
            assert output_path.exists()
            assert output_path.stat().st_size > 0
        except (ImportError, ModuleNotFoundError) as e:
            pytest.skip(f"Required CDPL module not found: {e}. Skipping test.")
