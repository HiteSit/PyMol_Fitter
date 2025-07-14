"""
Tests for the Docking Engine module in pymol_fitter_src.

These tests focus on the scientific functionality of the docking engine,
independent of the Flask API layer.
"""
import os
import shutil
import sys
import pytest
from pathlib import Path
import tempfile
import numpy as np

from pymol_fitter_server.pymol_fitter_src.Docking_Engine import (
    Pymol_Docking,
    outer_minimization
)

class TestPymolDocking:
    """Test suite for the Pymol_Docking class."""
    
    @pytest.fixture(scope="class")
    def example_files(self, example_dir):
        """Return paths to sample files from the examples directory."""
        return {
            "protein": example_dir / "8gcy.pdb",
            "ligand_sdf": example_dir / "Ligand.sdf",
            "crystal": example_dir / "Crystal.sdf",
            "smiles": "COc1cc(OC)c(cc1OC)C(=O)CCCCCCn1c(=O)c(-c2cc(OC)c(OC)c(OC)c2)cc2cc(OC)c(OC)c(OC)c21"
        }
    
    @pytest.fixture
    def docking_instance(self, example_files):
        """Create a Pymol_Docking instance with sample files."""
        return Pymol_Docking(
            protein_pdb=str(example_files["protein"]),
            input_ligands=str(example_files["ligand_sdf"]),
            crystal_sdf=str(example_files["crystal"]),
            is_smiles=False
        )
    
    @pytest.fixture
    def docking_instance_smiles(self, example_files):
        """Create a Pymol_Docking instance with SMILES input."""
        return Pymol_Docking(
            protein_pdb=str(example_files["protein"]),
            input_ligands=example_files["smiles"],
            crystal_sdf=str(example_files["crystal"]),
            is_smiles=True
        )
    
    def test_initialization(self, docking_instance, example_files):
        """Test that the Pymol_Docking class initializes correctly."""
        assert docking_instance.protein_pdb == Path(example_files["protein"])
        assert docking_instance.ligands_sdf == Path(example_files["ligand_sdf"])
        assert docking_instance.crystal_sdf == Path(example_files["crystal"])
        assert docking_instance.input_mode == "SDF"
    
    def test_initialization_smiles(self, docking_instance_smiles, example_files):
        """Test initialization with SMILES input."""
        assert docking_instance_smiles.protein_pdb == Path(example_files["protein"])
        assert docking_instance_smiles.ligands_smiles == example_files["smiles"]
        assert docking_instance_smiles.crystal_sdf == Path(example_files["crystal"])
        assert docking_instance_smiles.input_mode == "SMILES"
    
    def test_filter_protein(self, docking_instance, example_files, temp_dir):
        """Test the filter_protein static method."""
        # Save original working directory
        original_dir = os.getcwd()
        
        try:
            # Change to temp dir
            os.chdir(temp_dir)
            
            # Verify the method returns a Path
            filtered_path = Pymol_Docking.filter_protein(example_files["protein"])
            assert isinstance(filtered_path, Path)
            assert filtered_path.exists()
            
            # Clean up the temporary file
            if filtered_path.exists():
                os.unlink(filtered_path)
        finally:
            # Return to original directory
            os.chdir(original_dir)
    
    @pytest.mark.parametrize("mode", ["Dock", "Minimize"])
    def test_cdpk_fixer(self, mode, example_files):
        """Test the cdpk_fixer static method with different modes."""
        # Test with valid mode
        fixed_path = Pymol_Docking.cdpk_fixer(example_files["ligand_sdf"], mode)
        assert isinstance(fixed_path, Path)
        assert fixed_path.exists()
        
        # Clean up
        if fixed_path.exists():
            os.unlink(fixed_path)
    
    def test_cdpk_fixer_invalid_mode(self, example_files):
        """Test cdpk_fixer with invalid mode."""
        with pytest.raises(ValueError, match="Mode must be either 'Dock' or 'Minimize'"):
            Pymol_Docking.cdpk_fixer(example_files["ligand_sdf"], "InvalidMode")
    
    def test_cdpk_fixer_nonexistent_file(self):
        """Test cdpk_fixer with non-existent file."""
        with pytest.raises(FileNotFoundError):
            Pymol_Docking.cdpk_fixer(Path("nonexistent.sdf"), "Dock")
    
    # NOTE: The following tests would run the actual docking/minimization
    # They are marked as slow and may be skipped during quick testing
    
    @pytest.mark.slow
    def test_prepare_protein(self, docking_instance):
        """Test protein preparation (slow)."""
        prepared_path = docking_instance.prepare_protein()
        assert prepared_path.exists()
        assert prepared_path.suffix == ".pdb"
        
        # Clean up the temporary file
        if prepared_path.exists():
            os.unlink(prepared_path)
    
    @pytest.mark.slow
    @pytest.mark.parametrize("mode", ["Dock", "Minimize"])
    def test_prepare_ligands(self, docking_instance, mode):
        """Test ligand preparation with SDF input (slow)."""
        fixed_ligands, fixed_crystal = docking_instance.prepare_ligands(mode)
        assert fixed_ligands.exists()
        assert fixed_crystal.exists()
        
        # Clean up the temporary files
        if fixed_ligands.exists() and fixed_ligands != docking_instance.ligands_sdf:
            os.unlink(fixed_ligands)
    
    @pytest.mark.slow
    def test_prepare_ligands_smiles(self, docking_instance_smiles):
        """Test ligand preparation with SMILES input (slow)."""
        fixed_ligands, fixed_crystal = docking_instance_smiles.prepare_ligands("Dock")
        assert fixed_ligands.exists()
        assert fixed_crystal.exists()
        
        # Clean up the temporary file
        if fixed_ligands.exists():
            os.unlink(fixed_ligands)
    
    @pytest.mark.slow
    @pytest.mark.parametrize("mode", ["Dock", "Minimize"])
    def test_run_smina_docking(self, docking_instance, mode, temp_dir):
        """Test running smina docking (very slow)."""
        # Save original working directory
        original_dir = os.getcwd()
        
        try:
            # Change to the temp dir for output files
            os.chdir(temp_dir)
            
            # Create a new docking instance that will use the temp_dir as workdir
            # We need to recreate it after changing directories so its workdir is set to temp_dir
            protein_path = docking_instance.protein_pdb
            if docking_instance.input_mode == "SDF":
                ligand_input = str(docking_instance.ligands_sdf)
            else:
                ligand_input = docking_instance.ligands_smiles
            
            temp_docking_instance = Pymol_Docking(
                protein_pdb=str(protein_path),
                input_ligands=ligand_input,
                crystal_sdf=str(docking_instance.crystal_sdf),
                is_smiles=(docking_instance.input_mode == "SMILES")
            )
            
            # Run the docking
            smina_output, protein_prep, log_path = temp_docking_instance.run_smina_docking(
                mode, f"test_{mode.lower()}"
            )
            
            assert smina_output is not None, f"Smina output path: {smina_output}"
            assert protein_prep is not None, f"Protein prep path: {protein_prep}"
            assert log_path is not None, f"Log path: {log_path}"
            
            assert smina_output.exists()
            assert protein_prep.exists()
            assert log_path.exists()
            
            # Basic content check
            assert smina_output.stat().st_size > 0
            
        finally:
            # Return to original directory
            os.chdir(original_dir)
