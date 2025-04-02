# Standard library imports
import logging
import os
import subprocess
from pathlib import Path
from tempfile import _TemporaryFileWrapper, gettempdir, NamedTemporaryFile
from typing import Literal, Union, Optional, Dict, Any, Tuple

# Third party imports
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
import biotite.structure as struc
from biotite.structure.io import pdb
import CDPL.Chem as Chem
import CDPL.ConfGen as ConfGen
import datamol as dm
from openbabel import pybel
import numpy as np

import biotite.structure as struc
from biotite.structure.io import pdb

from pymol import cmd
import rdkit.Chem as rdChem
from rdkit.Chem.EnumerateStereoisomers import (
    EnumerateStereoisomers,
    StereoEnumerationOptions
)

# Local imports
from .Protein_Preparation import ProteinPreparation_Protoss, ProteinPreparation_PDBFixer
from .Protein_Minimization import minimize_complex
from .CDPK_Utils import CDPK_Runner

# Set up logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Create logger for this module
logger = logging.getLogger(__name__)

class Pymol_Docking:
    def __init__(self, protein_pdb: str, input_ligands: str, crystal_sdf: Optional[str] = None, is_smiles: bool = False):
        self.workdir: Path = Path(os.getcwd())
        self.protein_pdb: Path = Path(protein_pdb)
        if crystal_sdf is None:
            self.crystal_sdf: Path = Path("Crystal.sdf")
        else:
            self.crystal_sdf: Path = Path(crystal_sdf)
        
        self.protein_preared = None
        
        if is_smiles:
            self.ligands_smiles: str = input_ligands
            self.input_mode = "SMILES"
        else:
            self.ligands_sdf: Path = Path(input_ligands)
            self.input_mode = "SDF"
    
    @staticmethod
    def filter_protein(protein_pdb: Path) -> Path:
        struct_array = pdb.PDBFile.read(protein_pdb).get_structure(model=1)
        struct_array_filter = struct_array[struc.filter_amino_acids(struct_array)]
        
        struct_out = pdb.PDBFile()
        struct_out.set_structure(struct_array_filter)
        
        with NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_file:
            struct_out.write(tmp_file.name)
            return Path(tmp_file.name)
        
    def prepare_protein(self) -> Path:
        try:
            protein_basename: str = self.protein_pdb.stem
            
            protein_FILT = self.filter_protein(self.protein_pdb)
            protein_PREP: Path = self.workdir / f"{protein_basename}_Prep.pdb"
            
            pp = ProteinPreparation_Protoss()
            prepared_protein = pp(protein_FILT, protein_PREP)
            
            self.protein_preared = prepared_protein
            return prepared_protein
        except Exception as e:
            logger.warning("Error in prepare_protein: Protoss failed")
            logger.info("Trying PDBFixer")
            pp = ProteinPreparation_PDBFixer()
            prepared_protein = pp(protein_FILT, protein_PREP)
            
            self.protein_preared = prepared_protein
            return prepared_protein
                
    @staticmethod
    def cdpk_fixer(input_sdf: Path, mode: Literal["Dock", "Minimize"]) -> None:
        TMP_fixed_sdf: _TemporaryFileWrapper[bytes] = NamedTemporaryFile(suffix=".sdf", delete=False)
        
        if mode == "Dock":
            cdpk_runner = CDPK_Runner(standardize=True, protonate=True, gen3d=True)
            cdpk_runner.prepare_ligands(input_sdf, TMP_fixed_sdf.name)
            
            print(TMP_fixed_sdf.name)
            
            return Path(TMP_fixed_sdf.name)
        
        elif mode == "Minimize":
            cdpk_runner = CDPK_Runner(standardize=True, protonate=True, gen3d=False)
            cdpk_runner.prepare_ligands(input_sdf, TMP_fixed_sdf.name)
            
            print(TMP_fixed_sdf.name)
            
            return Path(TMP_fixed_sdf.name)
    
    def prepare_ligands(self, mode: Literal["Dock", "Minimize"]):
        assert mode in ["Dock", "Minimize"], "Mode must be either Dock or Minimize"
        
        def _rdkit_addhs(mol_path: Path):
            another_tmp = NamedTemporaryFile(suffix=".sdf", delete=False)
            mols = dm.read_sdf(mol_path)
            mols_H = []
            for mol in mols:
                mol_H = dm.add_hs(mol, add_coords=True)
                mols_H.append(mol_H)
            
            dm.to_sdf(mols_H, another_tmp.name)
            
            return Path(another_tmp.name)
        
        if mode == "Dock":
            if self.input_mode == "SDF":
                logger.info("Preparing ligands from SDF")
                fixed_crystal = self.crystal_sdf.resolve()
                fixed_ligands = self.cdpk_fixer(self.ligands_sdf, mode)
            
            elif self.input_mode == "SMILES":
                TMP_SMILES_SDF = Path(gettempdir()) / "TMP_SMILES.sdf"
                mol = dm.to_mol(self.ligands_smiles, sanitize=True, kekulize=True)
                mol.SetProp("_Name", "ligand")
                dm.to_sdf(mol, TMP_SMILES_SDF)
                
                fixed_crystal = self.crystal_sdf.resolve()
                fixed_ligands = self.cdpk_fixer(TMP_SMILES_SDF, mode)
                
            return _rdkit_addhs(fixed_ligands), fixed_crystal
            
        elif mode == "Minimize":
            logger.info("Preparing ligands from SDF")
            fixed_crystal = self.crystal_sdf.resolve()
            fixed_ligands = self.cdpk_fixer(self.ligands_sdf, mode)
            
            return _rdkit_addhs(fixed_ligands), fixed_crystal
        
    @staticmethod
    def pose_buster_processer(mol_pred: Path, mol_crystal: Path, mol_prot: Path):
        try:
            from posebusters import PoseBusters
            buster = PoseBusters()
            df = buster.bust(mol_pred, mol_crystal, mol_prot)
            df.reset_index(drop=False, inplace=True)
            bool_columns = df.select_dtypes(include=[bool]).columns
            success_rate_series = df[bool_columns].mean(axis=1) * 100
            success_rate_mean = round(np.mean(success_rate_series), 2)
            
            return mol_pred, success_rate_mean
        except Exception as e:
            logger.warning(f"Error in pose buster processer: {e}")
            return mol_pred, 0
            
    def run_smina_docking(self, mode: Literal["Dock", "Minimize"], docking_basename: str) -> Tuple[Path, Path]:
        # Fix the protein
        print("Running protein preparation")
        protein_PKA: Path = self.prepare_protein()

        # Fix the ligands
        try:
            fixed_ligands, fixed_crystal = self.prepare_ligands(mode)
        except Exception as e:
            logger.critical(f"Error in prepare_ligands: {e}")
            return None, None, None

        # Define the output
        smina_output: Path = self.workdir / f"{docking_basename}.sdf"
        smina_log = smina_output.with_suffix(".log")

        if mode == "Dock":
            smina_lst = [
                "smina",
                "-r", protein_PKA.as_posix(),
                "-l", fixed_ligands.as_posix(),
                "--autobox_ligand", fixed_crystal.as_posix(),
                "-o", str(smina_output),
                "--exhaustiveness", "32"
            ]

            print("Running docking")
            with open(smina_log, "w") as log_file:
                subprocess.run(smina_lst, check=True, stdout=log_file, stderr=log_file)

        elif mode == "Minimize":
            assert self.input_mode == "SDF", "Minimization only works with SDF input"
            smina_lst = [
                "smina",
                "-r", protein_PKA.as_posix(),
                "-l", fixed_ligands.as_posix(),
                "--autobox_ligand", fixed_crystal.as_posix(),
                "--minimize",
                "-o", smina_output,
            ]
            print("Running docking")
            with open(smina_log, "w") as log_file:
                subprocess.run(smina_lst, check=True, stdout=log_file, stderr=log_file)
                
        # Check the size of the output and check if it is > 0 bytes
        if smina_output.stat().st_size == 0:
            logger.critical("Docking failed")
            return None, None, None
        
        logger.info("Running pose buster")
        _, compliance_rate = self.pose_buster_processer(smina_output, fixed_crystal, protein_PKA)
        print(f"\n\nCompliance rate: {round(compliance_rate, 2)}\n\n")
        
        smina_mol = dm.read_sdf(smina_output, sanitize=False)[0]
        smina_mol.SetProp("_Name", "Ligand")
        smina_mol.SetProp("Buster_Compliace", f"{compliance_rate}")
        dm.to_sdf(smina_mol, smina_output)
        
        return smina_output, protein_PKA, smina_log

    def run_complex_minimization(self, protein_prep: Path, docked_ligand: Path):
        # Get the mol object
        mol = dm.read_sdf(docked_ligand)[0]
        
        minimizer = minimize_complex(protein_prep, mol)
        
        pdb_string_before: str | float = minimizer["PDB_BEFORE"]
        pdb_string_after: str | float = minimizer["PDB_AFTER"]
        
        def write_temp_pdb(content: str) -> str:
            """Write PDB content to a temporary file and return the file path."""

            tmp_file = NamedTemporaryFile(delete=False, suffix=".pdb")
            with open(tmp_file.name, "w") as f:
                f.write(content)
            
            return tmp_file.name

        temp_pdb_path_before = write_temp_pdb(pdb_string_before)
        temp_pdb_path_after = write_temp_pdb(pdb_string_after)

        struct_1 = pdb.PDBFile.read(temp_pdb_path_before).get_structure(model=1)
        struct_2 = pdb.PDBFile.read(temp_pdb_path_after).get_structure(model=1)
        
        multistate = struc.stack([struct_1, struct_2])
        pdb_multistate = pdb.PDBFile()
        pdb_multistate.set_structure(multistate)
        
        protein_basename: str = self.protein_pdb.stem
        protein_complex = self.workdir / f"{protein_basename}_Complex.pdb"
        pdb_multistate.write(protein_complex)
        
        return protein_complex