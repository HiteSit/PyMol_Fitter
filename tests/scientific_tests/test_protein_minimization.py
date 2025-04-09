"""
Tests for the Protein Minimization module in pymol_fitter_src.

These tests focus on the scientific functionality of the protein minimization,
independent of the Flask API layer.
"""
import os
import pytest
from pathlib import Path
import tempfile

import rdkit.Chem as rdChem
import datamol as dm

from pymol_fitter_server.pymol_fitter_src.Protein_Minimization import (
    minimize_complex
)


class TestProteinMinimization:
    """Test suite for the protein minimization module."""
    
    @pytest.fixture(
        scope="class",
        params=[
            {"protein": "LAC3.pdb", "ligand": "Lig_Min.sdf"},
            # {"protein": "8gcy.pdb", "ligand": "8gcy_Crystal.sdf"},
            # {"protein": "GLP-R.pdb", "ligand": "Complex_Ligand.sdf"}
        ],
        ids=["glpr"]
    )
    def example_files(self, request, example_dir):
        """Return paths to sample files from the examples directory.
        
        This fixture is parameterized to run tests with multiple sets of example files.
        """
        file_paths = {
            "protein": example_dir / request.param["protein"],
            "ligand": example_dir / request.param["ligand"]
        }
        
        # Skip if files don't exist
        if not file_paths["protein"].exists():
            pytest.skip(f"Protein file not found: {file_paths['protein']}")
        if not file_paths["ligand"].exists():
            pytest.skip(f"Ligand file not found: {file_paths['ligand']}")
            
        return file_paths
    
    @pytest.fixture(scope="function")
    def ligand_mol(self, example_files):
        """Return an RDKit molecule object for the example ligand."""
        mols = dm.read_sdf(example_files["ligand"])
        if not mols:
            pytest.skip("Could not read ligand SDF file")
        return mols[0]
    
    def test_minimize_complex_protein_only(self, example_files):
        """Test minimization with protein only (no ligand)."""
        try:
            result = minimize_complex(example_files["protein"], None)
            
            # Check the returned data structure
            assert isinstance(result, dict)
            assert "PDB_BEFORE" in result
            assert "PDB_AFTER" in result
            
            # Basic content checks on the PDB strings
            assert len(result["PDB_BEFORE"]) > 0
            assert len(result["PDB_AFTER"]) > 0
            assert "ATOM" in result["PDB_BEFORE"]
            assert "ATOM" in result["PDB_AFTER"]
            
            # Optional energy checks if available in the result
            if "E_INITIAL" in result and "E_FINAL" in result:
                # Final energy should be lower than initial energy
                assert result["E_FINAL"] <= result["E_INITIAL"]
        except (ImportError, ModuleNotFoundError) as e:
            pytest.skip(f"Required module not found: {e}. Skipping test.")
    
    @pytest.mark.slow
    def test_minimize_complex_with_ligand(self, example_files, ligand_mol):
        """Test minimization with both protein and ligand."""
        try:
            result = minimize_complex(example_files["protein"], ligand_mol)
            
            # Check the returned data structure
            assert isinstance(result, dict)
            assert "PDB_BEFORE" in result
            assert "PDB_AFTER" in result
            
            # Basic content checks on the PDB strings
            assert len(result["PDB_BEFORE"]) > 0
            assert len(result["PDB_AFTER"]) > 0
            assert "ATOM" in result["PDB_BEFORE"]
            assert "ATOM" in result["PDB_AFTER"]
            
            # Ligand should be included in the PDB
            assert "HETATM" in result["PDB_BEFORE"]
            assert "HETATM" in result["PDB_AFTER"]
            
            # Optional energy checks if available in the result
            if "E_INITIAL" in result and "E_FINAL" in result:
                # Final energy should be lower than initial energy
                assert result["E_FINAL"] <= result["E_INITIAL"]
        except (ImportError, ModuleNotFoundError) as e:
            pytest.skip(f"Required module not found: {e}. Skipping test.")
    
    def test_minimize_complex_nonexistent_protein(self, ligand_mol):
        """Test handling of non-existent protein file."""
        with pytest.raises(FileNotFoundError):
            minimize_complex(Path("nonexistent.pdb"), ligand_mol)
