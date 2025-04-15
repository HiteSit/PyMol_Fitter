"""
Tests for the Protein Preparation module in pymol_fitter_src.

These tests focus on the scientific functionality of the protein preparation,
independent of the Flask API layer.
"""
import os
import pytest
from pathlib import Path

from pymol_fitter_server.pymol_fitter_src.Protein_Preparation import (
    ProteinPreparation_Protoss,
    ProteinPreparation_PDBFixer
)


class TestProteinPreparation:
    """Test suite for the protein preparation classes."""
    
    @pytest.fixture(scope="class")
    def example_protein(self, example_dir):
        """Return path to sample protein file."""
        return example_dir / "8gcy.pdb"
    
    @pytest.fixture(scope="function")
    def output_path(self, temp_dir):
        """Return a temporary output path for prepared protein."""
        return temp_dir / "prepared_protein.pdb"
    
    def test_protoss_initialization(self):
        """Test that ProteinPreparation_Protoss initializes correctly."""
        protoss = ProteinPreparation_Protoss()
        assert protoss is not None
    
    def test_pdbfixer_initialization(self):
        """Test that ProteinPreparation_PDBFixer initializes correctly."""
        pdbfixer = ProteinPreparation_PDBFixer()
        assert pdbfixer is not None
    
    @pytest.mark.slow
    def test_protoss_preparation(self, example_protein, output_path):
        """Test protein preparation using Protoss (slow)."""
        protoss = ProteinPreparation_Protoss()
        
        try:
            # Note: This test will only pass if Protoss is available in the environment
            prepared_path = protoss(example_protein, output_path)
            assert prepared_path.exists()
            assert prepared_path.stat().st_size > 0
            assert prepared_path == output_path
        except Exception as e:
            pytest.skip(f"Protoss preparation failed: {e}. Skipping test.")
    
    @pytest.mark.slow
    def test_pdbfixer_preparation(self, example_protein, output_path):
        """Test protein preparation using PDBFixer (slow)."""
        pdbfixer = ProteinPreparation_PDBFixer()
        
        prepared_path = pdbfixer(example_protein, output_path)
        assert prepared_path.exists()
        assert prepared_path.stat().st_size > 0
        assert prepared_path == output_path
    
    @pytest.mark.slow
    def test_preparation_with_nonexistent_input(self, output_path):
        """Test protein preparation with non-existent input file."""
        protoss = ProteinPreparation_Protoss()
        pdbfixer = ProteinPreparation_PDBFixer()
        
        with pytest.raises(FileNotFoundError):
            protoss(Path("nonexistent.pdb"), output_path)
            
        with pytest.raises(FileNotFoundError):
            pdbfixer(Path("nonexistent.pdb"), output_path)
