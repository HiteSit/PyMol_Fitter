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
from .Protein_Preparation import ProteinPreparation_Protoss
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

        # if os.path.exists(input_ligands):
        #     self.ligands_sdf: Path = Path(input_ligands)
        #     self.input_mode = "SDF"
        # else:
        #     self.ligands_smiles: str = input_ligands
        #     self.input_mode = "SMILES"
    
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
        protein_basename: str = self.protein_pdb.stem
        
        protein_FILT = self.filter_protein(self.protein_pdb)
        protein_PREP: Path = self.workdir / f"{protein_basename}_Prep.pdb"
        
        pp = ProteinPreparation_Protoss()
        prepared_protein = pp(protein_FILT, protein_PREP)
        
        self.protein_preared = prepared_protein
        return prepared_protein

    @staticmethod
    def cdpk_fixer(input_sdf: Path, TMP_basename: str, mode: Literal["Dock", "Minimize"]) -> None:
        TMP_fixed_sdf = Path(gettempdir()) / "TMP_fixed.sdf"
        
        supplier = rdChem.SDMolSupplier(input_sdf.as_posix())
        writer = rdChem.SDWriter(TMP_fixed_sdf.as_posix())
        for mol in supplier:
            mol.SetProp("_Name", TMP_basename)
            writer.write(mol)
        writer.close()
        
        if mode == "Dock":
            cdpk_runner = CDPK_Runner(standardize=True, protonate=True, gen3d=True)
            cdpk_runner.prepare_ligands(TMP_fixed_sdf, TMP_fixed_sdf)
            
            return TMP_fixed_sdf.resolve()
        elif mode == "Minimize":
            cdpk_runner = CDPK_Runner(standardize=True, protonate=True, gen3d=False)
            cdpk_runner.prepare_ligands(TMP_fixed_sdf, TMP_fixed_sdf)
            
            return TMP_fixed_sdf.resolve()
    
    def prepare_ligands(self, mode: Literal["Dock", "Minimize"]):
        
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
                fixed_ligands = self.cdpk_fixer(self.ligands_sdf, "ligand", mode)
            
            elif self.input_mode == "SMILES":
                TMP_SMILES_SDF = Path(gettempdir()) / "TMP_SMILES.sdf"
                mol = dm.to_mol(self.ligands_smiles, sanitize=True, kekulize=True)
                mol.SetProp("_Name", "ligand")
                dm.to_sdf(mol, TMP_SMILES_SDF)
                
                fixed_crystal = self.crystal_sdf.resolve()
                fixed_ligands = self.cdpk_fixer(TMP_SMILES_SDF, "ligand", mode)
                
            return _rdkit_addhs(fixed_ligands), fixed_crystal
            
        elif mode == "Minimize":
            logger.info("Preparing ligands from SDF")
            fixed_crystal = self.crystal_sdf.resolve()
            fixed_ligands = self.cdpk_fixer(self.ligands_sdf, "ligand", mode)
            
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
            success_rate_list = success_rate_series.tolist()
            success_rate_mean = round(np.mean(success_rate_series), 2)
            
            # # Read all molecules first
            # mol_list = list(rdChem.SDMolSupplier(mol_pred.as_posix()))
            
            # # Create a temporary file path
            # temp_path = mol_pred.parent / f"temp_{mol_pred.name}"
            
            # # Write to temporary file
            # with rdChem.SDWriter(temp_path.as_posix()) as writer:
            #     for pose, success_rate in zip(mol_list, success_rate_list):
            #         pose.SetProp("Compliance", str(success_rate))  # Convert to string
            #         writer.write(pose)
            
            # # Replace original file with temporary file
            # temp_path.replace(mol_pred)
            
            return mol_pred, success_rate_mean
        except Exception as e:
            logger.error(f"Error in pose buster processer: {e}")
            return mol_pred, 0
            
    def run_smina_docking(self, mode: Literal["Dock", "Minimize"], docking_basename: str) -> Tuple[Path, Path]:
        # Fix the protein
        print("Running protein preparation")
        protein_PKA: Path = self.prepare_protein()

        # Fix the ligands
        fixed_ligands, fixed_crystal = self.prepare_ligands(mode)

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
        
        logger.info("Running pose buster")
        smina_bustered, compliance_rate = self.pose_buster_processer(smina_output, fixed_crystal, protein_PKA)
        print(f"\n\nCompliance rate: {round(compliance_rate, 2)}\n\n")
        
        return smina_bustered, protein_PKA, smina_log

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