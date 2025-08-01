# Standard library imports
import logging
import os
import subprocess
from pathlib import Path
from tempfile import _TemporaryFileWrapper, gettempdir, NamedTemporaryFile
from typing import Dict, Any, List, Literal, Optional, Tuple, Union
import csv
import multiprocessing as mp
from functools import partial

# Third party imports
from Bio.PDB import PDBParser
from scrubber import Scrub
from Bio.PDB.PDBIO import PDBIO
import biotite.structure as struc
from biotite.structure.io import pdb
import CDPL.Chem as Chem
import CDPL.ConfGen as ConfGen
import datamol as dm
from openbabel import pybel
import numpy as np
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
    """
    Class for handling molecular docking using PyMOL and external docking tools.

    This class provides methods for preparing proteins, ligands, and performing
    docking operations using the smina docking engine.

    Attributes:
        workdir: Working directory path
        protein_pdb: Path to the input protein PDB file
        crystal_sdf: Path to the crystal structure SDF file
        protein_preared: Path to the prepared protein PDB file (set after preparation)
        input_mode: Mode of input, either "SMILES" or "SDF"
        ligands_smiles: SMILES string of the ligand (if input_mode is "SMILES")
        ligands_sdf: Path to the ligand SDF file (if input_mode is "SDF")
    """

    # authoritative SDF tag used by smina for affinity (see docs & community threads)
    _SMINA_TAG = "minimizedAffinity"

    def __init__(
        self,
        protein_pdb: str,
        input_ligands: str,
        crystal_sdf: Optional[str] = None,
        is_smiles: bool = False,
        virtual_screening: bool = False,
    ):
        """
        Initialize the Pymol_Docking class.

        Args:
            protein_pdb: Path to the protein PDB file
            input_ligands: Either a SMILES string or path to ligand SDF file
            crystal_sdf: Path to the crystal SDF file (optional)
            is_smiles: Whether input_ligands is a SMILES string
            virtual_screening: Whether to enable virtual screening
        """
        self.workdir: Path = Path(os.getcwd())
        self.protein_pdb: Path = Path(protein_pdb)
        if crystal_sdf is None:
            self.crystal_sdf: Path = Path("Crystal.sdf")
        else:
            self.crystal_sdf: Path = Path(crystal_sdf)
        
        self.protein_preared = None
        
        self.virtual_screening = virtual_screening  # VS flag
        
        if is_smiles:
            self.ligands_smiles: str = input_ligands
            self.input_mode = "SMILES"
        else:
            self.ligands_sdf: Path = Path(input_ligands)
            self.input_mode = "SDF"
    
    @staticmethod
    def filter_protein(protein_pdb: Path) -> Path:
        """
        Filter protein structure to keep only amino acids.

        Args:
            protein_pdb: Path to the protein PDB file

        Returns:
            Path to the filtered protein PDB file

        Raises:
            FileNotFoundError: If the protein PDB file doesn't exist
            RuntimeError: If filtering fails
        """
        if not protein_pdb.exists():
            raise FileNotFoundError(f"Protein file not found: {protein_pdb}")
            
        try:
            struct_array = pdb.PDBFile.read(protein_pdb).get_structure(model=1)
            struct_array_filter = struct_array[struc.filter_amino_acids(struct_array)]
            
            struct_out = pdb.PDBFile()
            struct_out.set_structure(struct_array_filter)
            
            with NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_file:
                struct_out.write(tmp_file.name)
                return Path(tmp_file.name)
        except Exception as e:
            raise RuntimeError(f"Error filtering protein structure: {e}") from e
        
    def prepare_protein(self) -> Path:
        """
        Prepare protein structure for docking.

        This method attempts to prepare the protein using Protoss first,
        and falls back to PDBFixer if Protoss fails.

        Returns:
            Path to the prepared protein PDB file

        Raises:
            RuntimeError: If both preparation methods fail
        """
        try:
            protein_basename: str = self.protein_pdb.stem
            
            protein_FILT = self.filter_protein(self.protein_pdb)
            protein_PREP: Path = self.workdir / f"{protein_basename}_Prep.pdb"
            
            pp = ProteinPreparation_Protoss()
            prepared_protein = pp(protein_FILT, protein_PREP)
            
            self.protein_preared = prepared_protein
            return prepared_protein
        except Exception as e:
            logger.warning(f"Error in prepare_protein: Protoss failed - {e}")
            logger.info("Trying PDBFixer")
            
            try:
                pp = ProteinPreparation_PDBFixer()
                prepared_protein = pp(protein_FILT, protein_PREP)
                
                self.protein_preared = prepared_protein
                return prepared_protein
            except Exception as e2:
                raise RuntimeError(f"Both protein preparation methods failed. Last error: {e2}") from e2
                
    @staticmethod
    def cdpk_fixer(input_sdf: Path, mode: Literal["Dock", "Minimize"]) -> Path:
        """
        Fix ligand structures using CDPK utilities.

        Args:
            input_sdf: Path to the input SDF file
            mode: Mode of operation, either "Dock" or "Minimize"

        Returns:
            Path to the fixed SDF file

        Raises:
            ValueError: If mode is not "Dock" or "Minimize"
            FileNotFoundError: If input SDF file doesn't exist
        """
        if mode not in ["Dock", "Minimize"]:
            raise ValueError("Mode must be either 'Dock' or 'Minimize'")
            
        if not input_sdf.exists():
            raise FileNotFoundError(f"Input SDF file not found: {input_sdf}")
        
        try:
            TMP_fixed_sdf: _TemporaryFileWrapper[bytes] = NamedTemporaryFile(suffix=".sdf", delete=False)
            
            if mode == "Dock":
                cdpk_runner = CDPK_Runner(standardize=True, protonate=True, gen3d=True)
                cdpk_runner.prepare_ligands(input_sdf, TMP_fixed_sdf.name)
                
            elif mode == "Minimize":
                cdpk_runner = CDPK_Runner(standardize=True, protonate=True, gen3d=False)
                cdpk_runner.prepare_ligands(input_sdf, TMP_fixed_sdf.name)
                
            logger.info(f"CDPK fixed ligands saved to: {TMP_fixed_sdf.name}")
            return Path(TMP_fixed_sdf.name)
            
        except Exception as e:
            if TMP_fixed_sdf and os.path.exists(TMP_fixed_sdf.name):
                os.unlink(TMP_fixed_sdf.name)
            raise RuntimeError(f"Error in CDPK fixer: {e}") from e
    
    @staticmethod
    def meeko_fixer(input_sdf: Path, mode: Literal["Dock", "Minimize"]) -> Path:
        """
        Fix ligand structures using Meeko/Scrub for protonation.

        Args:
            input_sdf: Path to the input SDF file
            mode: Mode of operation, either "Dock" or "Minimize" (parameter maintained for compatibility)

        Returns:
            Path to the fixed SDF file

        Raises:
            ValueError: If mode is not "Dock" or "Minimize"
            FileNotFoundError: If input SDF file doesn't exist
        """
        if mode not in ["Dock", "Minimize"]:
            raise ValueError("Mode must be either 'Dock' or 'Minimize'")
            
        if not input_sdf.exists():
            raise FileNotFoundError(f"Input SDF file not found: {input_sdf}")
        
        try:
            TMP_fixed_sdf: _TemporaryFileWrapper[bytes] = NamedTemporaryFile(suffix=".sdf", delete=False)
            
            # Mode parameter is kept for API compatibility but doesn't affect functionality
            # in this implementation
            
            scrub = Scrub(
                ph_low=7.2,  # Default pH
                ph_high=7.4,
            )
            
            # Read molecule from SDF
            mols = dm.read_sdf(input_sdf)
            if not mols or len(mols) == 0:
                raise ValueError(f"No molecules found in {input_sdf}")
            
            mol = mols[0]
            if mol is None:
                raise ValueError(f"Could not read molecule from {input_sdf}")
            
            mols_states = []
            for mol_state in scrub(mol):
                mols_states.append(mol_state)
            
            if not mols_states or len(mols_states) == 0:
                raise ValueError(f"Scrub did not generate any molecular states for {input_sdf}")
            
            best_mol_state = mols_states[0]
            dm.to_sdf(best_mol_state, TMP_fixed_sdf.name)
            
            logger.info(f"Meeko fixed ligands saved to: {TMP_fixed_sdf.name}")
            return Path(TMP_fixed_sdf.name)
            
        except Exception as e:
            if TMP_fixed_sdf and os.path.exists(TMP_fixed_sdf.name):
                os.unlink(TMP_fixed_sdf.name)
            raise RuntimeError(f"Error in Meeko fixer: {e}") from e
    
    # def prepare_ligands(self, mode: Literal["Dock", "Minimize"]) -> Tuple[Path, Path]:
    #     """
    #     Prepare ligands for docking or minimization.

    #     Args:
    #         mode: Mode of operation, either "Dock" or "Minimize"

    #     Returns:
    #         Tuple of (fixed_ligands_path, fixed_crystal_path)

    #     Raises:
    #         ValueError: If mode is invalid or input mode is incompatible 
    #         RuntimeError: If ligand preparation fails
    #     """
    #     if mode not in ["Dock", "Minimize"]:
    #         raise ValueError("Mode must be either 'Dock' or 'Minimize'")
        
    #     def _rdkit_addhs(mol_path: Path) -> Path:
    #         """Add hydrogens to molecules using RDKit."""
    #         try:
    #             tmp_file = NamedTemporaryFile(suffix=".sdf", delete=False)
    #             mols = dm.read_sdf(mol_path)
    #             if not mols:
    #                 raise ValueError(f"No molecules found in {mol_path}")
                    
    #             mols_H = []
    #             for mol in mols:
    #                 mol_H = dm.add_hs(mol, add_coords=True)
    #                 mols_H.append(mol_H)
                
    #             dm.to_sdf(mols_H, tmp_file.name)
    #             return Path(tmp_file.name)
    #         except Exception as e:
    #             if os.path.exists(tmp_file.name):
    #                 os.unlink(tmp_file.name)
    #             raise RuntimeError(f"Error adding hydrogens: {e}") from e
        
    #     try:
    #         if mode == "Dock":
    #             if self.input_mode == "SDF":
    #                 logger.info("Preparing ligands from SDF")
    #                 fixed_crystal = self.crystal_sdf.resolve()
    #                 fixed_ligands = self.cdpk_fixer(self.ligands_sdf, mode)
                
    #             elif self.input_mode == "SMILES":
    #                 logger.info("Preparing ligands from SMILES")
    #                 TMP_SMILES_SDF = Path(gettempdir()) / "TMP_SMILES.sdf"
    #                 mol = dm.to_mol(self.ligands_smiles, sanitize=True, kekulize=True)
    #                 if mol is None:
    #                     raise ValueError(f"Failed to convert SMILES to molecule: {self.ligands_smiles}")
                        
    #                 mol.SetProp("_Name", "ligand")
    #                 dm.to_sdf(mol, TMP_SMILES_SDF)
                    
    #                 fixed_crystal = self.crystal_sdf.resolve()
    #                 fixed_ligands = self.cdpk_fixer(TMP_SMILES_SDF, mode)
                    
    #             return _rdkit_addhs(fixed_ligands), fixed_crystal
                
    #         elif mode == "Minimize":
    #             if self.input_mode != "SDF":
    #                 raise ValueError("Minimization only works with SDF input")
                    
    #             logger.info("Preparing ligands for minimization")
    #             fixed_crystal = self.crystal_sdf.resolve()
    #             fixed_ligands = self.cdpk_fixer(self.ligands_sdf, mode)
                
    #             return _rdkit_addhs(fixed_ligands), fixed_crystal
                
    #     except Exception as e:
    #         raise RuntimeError(f"Error in prepare_ligands: {e}") from e
        
    def prepare_ligands(self, mode: Literal["Dock", "Minimize"]) -> Tuple[Path, Path]:
        """
        Prepare ligands for docking or minimization using Meeko fixer.

        Args:
            mode: Mode of operation, either "Dock" or "Minimize"

        Returns:
            Tuple of (fixed_ligands_path, fixed_crystal_path)

        Raises:
            ValueError: If mode is invalid or input mode is incompatible 
            RuntimeError: If ligand preparation fails
        """
        if mode not in ["Dock", "Minimize"]:
            raise ValueError("Mode must be either 'Dock' or 'Minimize'")
        
        def _rdkit_addhs(mol_path: Path) -> Path:
            """Add hydrogens to molecules using RDKit."""
            try:
                tmp_file = NamedTemporaryFile(suffix=".sdf", delete=False)
                mols = dm.read_sdf(mol_path)
                if not mols:
                    raise ValueError(f"No molecules found in {mol_path}")
                    
                mols_H = []
                for mol in mols:
                    mol_H = dm.add_hs(mol, add_coords=True)
                    mols_H.append(mol_H)
                
                dm.to_sdf(mols_H, tmp_file.name)
                return Path(tmp_file.name)
            except Exception as e:
                if os.path.exists(tmp_file.name):
                    os.unlink(tmp_file.name)
                raise RuntimeError(f"Error adding hydrogens: {e}") from e
        
        try:
            if mode == "Dock":
                if self.input_mode == "SDF":
                    logger.info("Preparing ligands from SDF using Meeko")
                    fixed_crystal = self.crystal_sdf.resolve()
                    fixed_ligands = self.meeko_fixer(self.ligands_sdf, mode)
                
                elif self.input_mode == "SMILES":
                    logger.info("Preparing ligands from SMILES using Meeko")
                    TMP_SMILES_SDF = Path(gettempdir()) / "TMP_SMILES.sdf"
                    mol = dm.to_mol(self.ligands_smiles, sanitize=True, kekulize=True)
                    if mol is None:
                        raise ValueError(f"Failed to convert SMILES to molecule: {self.ligands_smiles}")
                        
                    mol.SetProp("_Name", "ligand")
                    dm.to_sdf(mol, TMP_SMILES_SDF)
                    
                    fixed_crystal = self.crystal_sdf.resolve()
                    fixed_ligands = self.meeko_fixer(TMP_SMILES_SDF, mode)
                    
                return _rdkit_addhs(fixed_ligands), fixed_crystal
                
            elif mode == "Minimize":
                if self.input_mode != "SDF":
                    raise ValueError("Minimization only works with SDF input")
                    
                logger.info("Preparing ligands for minimization using Meeko")
                fixed_crystal = self.crystal_sdf.resolve()
                fixed_ligands = self.meeko_fixer(self.ligands_sdf, mode)
                
                return _rdkit_addhs(fixed_ligands), fixed_crystal
                
        except Exception as e:
            raise RuntimeError(f"Error in prepare_ligands (Meeko): {e}") from e
        
    @staticmethod
    def pose_buster_processer(mol_pred: Path, mol_crystal: Path, mol_prot: Path) -> Tuple[Path, float]:
        """
        Process docking results using PoseBusters to calculate pose compliance.

        Args:
            mol_pred: Path to the predicted molecule SDF file
            mol_crystal: Path to the crystal structure SDF file
            mol_prot: Path to the protein PDB file

        Returns:
            Tuple of (molecule_path, compliance_rate)
        """
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
            return mol_pred, 0.0
            
    def run_smina_docking(
        self, 
        mode: Literal["Dock", "Minimize"], 
        docking_basename: str
    ) -> Tuple[Optional[Path], Optional[Path], Optional[Path]]:
        """
        Run smina docking or minimization.

        Args:
            mode: Mode of operation, either "Dock" or "Minimize"
            docking_basename: Base name for output files

        Returns:
            Tuple of (smina_output_path, prepared_protein_path, log_path)
            Returns None values if any step fails

        Raises:
            ValueError: If mode is invalid
            RuntimeError: If a critical error occurs during docking
        """
        if mode not in ["Dock", "Minimize"]:
            raise ValueError("Mode must be either 'Dock' or 'Minimize'")
            
        try:
            # Fix the protein
            logger.info("Running protein preparation")
            protein_PKA: Path = self.prepare_protein()

            # Fix the ligands
            try:
                fixed_ligands, fixed_crystal = self.prepare_ligands(mode)
            except Exception as e:
                logger.critical(f"Error in prepare_ligands: {e}")
                return None, None, None

            # Define the output paths
            smina_output: Path = self.workdir / f"{docking_basename}.sdf"
            smina_log = smina_output.with_suffix(".log")

            if mode == "Dock":
                # Prepare docking command
                smina_lst = [
                    "smina",
                    "-r", protein_PKA.as_posix(),
                    "-l", fixed_ligands.as_posix(),
                    "--autobox_ligand", fixed_crystal.as_posix(),
                    "-o", str(smina_output),
                    "--exhaustiveness", "32"
                ]

                logger.info("Running docking")
                with open(smina_log, "w") as log_file:
                    try:
                        subprocess.run(smina_lst, check=True, stdout=log_file, stderr=log_file)
                    except subprocess.CalledProcessError as e:
                        logger.error(f"Smina docking process failed: {e}")
                        return None, None, None

            elif mode == "Minimize":
                # Verify input mode is SDF
                if self.input_mode != "SDF":
                    raise ValueError("Minimization only works with SDF input")
                    
                # Prepare minimization command
                smina_lst = [
                    "smina",
                    "-r", protein_PKA.as_posix(),
                    "-l", fixed_ligands.as_posix(),
                    "--autobox_ligand", fixed_crystal.as_posix(),
                    "--minimize",
                    "-o", str(smina_output),
                ]
                
                logger.info("Running minimization")
                with open(smina_log, "w") as log_file:
                    try:
                        subprocess.run(smina_lst, check=True, stdout=log_file, stderr=log_file)
                    except subprocess.CalledProcessError as e:
                        logger.error(f"Smina minimization process failed: {e}")
                        return None, None, None
                    
            # Check if output file was created and is not empty
            if not smina_output.exists() or smina_output.stat().st_size == 0:
                logger.critical("Docking failed: output file is empty or does not exist")
                return None, None, None
            
            # Run pose buster to evaluate pose quality
            logger.info("Running pose buster")
            _, compliance_rate = self.pose_buster_processer(smina_output, fixed_crystal, protein_PKA)
            logger.info(f"Compliance rate: {round(compliance_rate, 2)}")
            
            # # Add compliance rate property to output molecule
            # try:
            #     smina_mols = dm.read_sdf(smina_output, sanitize=False)
            #     if not smina_mols or len(smina_mols) == 0:
            #         logger.warning("No molecules found in smina output")
            #     else:
            #         smina_mol = smina_mols[0]
            #         smina_mol.SetProp("_Name", "Ligand")
            #         smina_mol.SetProp("Buster_Compliace", f"{compliance_rate}")
            #         dm.to_sdf(smina_mol, smina_output)
            # except Exception as e:
            #     logger.warning(f"Error updating molecule properties: {e}")
            
            return smina_output, protein_PKA, smina_log
            
        except Exception as e:
            logger.error(f"Error in run_smina_docking: {e}")
            return None, None, None

    def run_complex_minimization(self, protein_prep: Path, docked_ligand: Optional[Path] = None) -> Path:
        """
        Run molecular dynamics minimization on a protein-ligand complex.

        Args:
            protein_prep: Path to the prepared protein PDB file
            docked_ligand: Path to the docked ligand SDF file (optional)

        Returns:
            Path to the minimized complex PDB file

        Raises:
            FileNotFoundError: If input files don't exist
            ValueError: If no molecules are found in the ligand file when provided
            RuntimeError: If minimization fails
        """
        if not protein_prep.exists():
            raise FileNotFoundError(f"Protein file not found: {protein_prep}")
            
        try:
            # Get the mol object if ligand is provided
            mol = None
            if docked_ligand is not None and docked_ligand.exists():
                mols = dm.read_sdf(docked_ligand)
                if not mols or len(mols) == 0:
                    raise ValueError(f"No molecules found in ligand file: {docked_ligand}")
                mol = mols[0]
            
            # Run minimization
            minimizer = minimize_complex(protein_prep, mol)
            
            pdb_string_before: str = minimizer["PDB_BEFORE"]
            pdb_string_after: str = minimizer["PDB_AFTER"]
            
            def write_temp_pdb(content: str) -> str:
                """Write PDB content to a temporary file and return the file path."""
                tmp_file = NamedTemporaryFile(delete=False, suffix=".pdb")
                with open(tmp_file.name, "w") as f:
                    f.write(content)
                return tmp_file.name

            # Create temporary PDB files
            temp_pdb_path_before = write_temp_pdb(pdb_string_before)
            temp_pdb_path_after = write_temp_pdb(pdb_string_after)

            # Create a multi-state PDB with both conformations
            struct_1 = pdb.PDBFile.read(temp_pdb_path_before).get_structure(model=1)
            struct_2 = pdb.PDBFile.read(temp_pdb_path_after).get_structure(model=1)
            
            multistate = struc.stack([struct_1, struct_2])
            pdb_multistate = pdb.PDBFile()
            pdb_multistate.set_structure(multistate)
            
            # Save the complex
            protein_basename: str = self.protein_pdb.stem
            protein_complex = self.workdir / f"{protein_basename}_Complex.pdb"
            pdb_multistate.write(protein_complex)
            
            # Clean up temporary files
            for path in [temp_pdb_path_before, temp_pdb_path_after, protein_prep]:
                if os.path.exists(path):
                    os.unlink(path)
            
            return protein_complex
            
        except Exception as e:
            raise RuntimeError(f"Error in complex minimization: {e}") from e

    # ────────────────── explode to per‑ligand SDFs (VS) ─────────────────

    def _explode_multi_ligands(self, tmpdir: Path) -> List[Tuple[Path, str]]:
        """Return list of (sdf_path, ligand_name) for VS."""
        ligs: List[Tuple[Path, str]] = []
        if self.input_mode == "SDF":
            mols = dm.read_sdf(self.ligands_sdf, sanitize=False)
            assert len(mols) > 1, "virtual_screening=True but SDF has ≤1 record"
            for i, m in enumerate(mols, 1):
                name = m.GetProp("_Name") or f"lig_{i}"
                p = tmpdir / f"{name}.sdf"; dm.to_sdf(m, p); ligs.append((p, name))
        else:
            smi_list = [s for s in self.ligands_smiles.split(".") if s]
            assert len(smi_list) > 1, "virtual_screening=True but only one SMILES"
            for i, smi in enumerate(smi_list, 1):
                mol = dm.to_mol(smi, sanitize=True, kekulize=True)
                mol.SetProp("_Name", f"lig_{i}")
                p = tmpdir / f"lig_{i}.sdf"; dm.to_sdf(mol, p); ligs.append((p, f"lig_{i}"))
        return ligs

    # ────────────────────────── smina wrapper ──────────────────────────

    @staticmethod
    def _smina_single(
        protein: Path,
        crystal: Path,
        lig_tuple: Tuple[Path, str],
        exhaust: int,
        outdir: Path,
    ) -> Dict[str, object]:
        lig_path, lig_name = lig_tuple
        out_base = outdir / lig_name
        out_sdf, out_log = out_base.with_suffix(".sdf"), out_base.with_suffix(".log")

        cmd_lst = [
            "smina", "-r", str(protein), "-l", str(lig_path),
            "--autobox_ligand", str(crystal),
            "--exhaustiveness", str(exhaust), "--cpu", "1",
            "-o", str(out_sdf),
        ]

        with open(out_log, "w") as lf:
            subprocess.run(cmd_lst, check=True, stdout=lf, stderr=lf)

        # fetch affinity if present
        aff: Optional[float] = None
        try:
            mol = dm.read_sdf(out_sdf, sanitize=False)[0]
            if mol.HasProp(Pymol_Docking._SMINA_TAG):
                aff = float(mol.GetProp(Pymol_Docking._SMINA_TAG))
        except Exception:
            pass

        return {"name": lig_name, "sdf": out_sdf, "log": out_log, "aff": aff}

    # ─────────────────────────── public VS API ──────────────────────────

    def run_virtual_screen(self, basename: str) -> Dict[str, object]:
        """Execute docking‑only screen across **multiple** ligands."""
        if not self.virtual_screening:
            raise RuntimeError("Instantiate with virtual_screening=True to run VS")

        protein = self.prepare_protein()
        crystal = self.crystal_sdf.resolve()

        tmpdir = Path(gettempdir()) / f"vs_{os.getpid()}"; tmpdir.mkdir(exist_ok=True)
        ligand_list = self._explode_multi_ligands(tmpdir)

        cpu = max(1, int(os.cpu_count() * 0.9))
        worker = partial(self._smina_single, protein, crystal, exhaust=8, outdir=self.workdir)

        with mp.Pool(cpu) as pool:
            res = pool.map(worker, ligand_list)

        # write CSV summary
        csv_path = self.workdir / f"{basename}_summary.csv"
        with open(csv_path, "w", newline="") as fh:
            wr = csv.writer(fh); wr.writerow(["ligand", Pymol_Docking._SMINA_TAG])
            for r in res:
                wr.writerow([r["name"], r["aff"]])

        return {
            "success": True,
            "docked_ligands": [r["sdf"] for r in res],
            "logs": [r["log"] for r in res],
            "summary": csv_path,
            "prepared_protein": protein,
        }


def outer_minimization(
    protein_prep: Path, 
    docked_ligand: Optional[Path] = None, 
    complex_path: Path = Path("minimized_complex.pdb")
) -> None:
    """
    Perform minimization of a protein-ligand complex and write to an output file.

    Args:
        protein_prep: Path to the prepared protein PDB file
        docked_ligand: Path to the docked ligand SDF file (optional)
        complex_path: Path where the minimized complex will be saved

    Raises:
        FileNotFoundError: If protein file doesn't exist or ligand file doesn't exist when provided
        ValueError: If no molecules are found in the ligand file when provided
        RuntimeError: If minimization fails
    """
    if not protein_prep.exists():
        raise FileNotFoundError(f"Protein file not found: {protein_prep}")
    
    def filter_protein(protein_pdb: Path) -> Path:
        """Filter protein structure to keep only amino acids."""
        struct_array = pdb.PDBFile.read(protein_pdb).get_structure(model=1)
        struct_array_filter = struct_array[struc.filter_amino_acids(struct_array)]
        
        struct_out = pdb.PDBFile()
        struct_out.set_structure(struct_array_filter)
        
        with NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_file:
            struct_out.write(tmp_file.name)
            return Path(tmp_file.name)
    
    try:
        # Get the mol object if ligand is provided
        mol = None
        if docked_ligand is not None and docked_ligand.exists():
            if not docked_ligand.exists():
                raise FileNotFoundError(f"Ligand file not found: {docked_ligand}")
                
            mols = dm.read_sdf(docked_ligand)
            if not mols or len(mols) == 0:
                raise ValueError(f"No molecules found in ligand file: {docked_ligand}")
            mol = mols[0]
        
        # Filter the protein and run minimization
        filtered_protein = filter_protein(protein_prep)
        minimizer = minimize_complex(filtered_protein, mol)
        
        pdb_string_before: str = minimizer["PDB_BEFORE"]
        pdb_string_after: str = minimizer["PDB_AFTER"]
        
        def write_temp_pdb(content: str) -> str:
            """Write PDB content to a temporary file and return the file path."""
            tmp_file = NamedTemporaryFile(delete=False, suffix=".pdb")
            with open(tmp_file.name, "w") as f:
                f.write(content)
            return tmp_file.name

        # Create temporary PDB files
        temp_pdb_path_before = write_temp_pdb(pdb_string_before)
        temp_pdb_path_after = write_temp_pdb(pdb_string_after)

        # Create a multi-state PDB with both conformations
        struct_1 = pdb.PDBFile.read(temp_pdb_path_before).get_structure(model=1)
        struct_2 = pdb.PDBFile.read(temp_pdb_path_after).get_structure(model=1)
        
        multistate = struc.stack([struct_1, struct_2])
        pdb_multistate = pdb.PDBFile()
        pdb_multistate.set_structure(multistate)
        
        # Save the complex
        pdb_multistate.write(complex_path)
        
        # Clean up temporary files
        for path in [temp_pdb_path_before, temp_pdb_path_after, filtered_protein]:
            if os.path.exists(path):
                os.unlink(path)
                
    except Exception as e:
        raise RuntimeError(f"Error during outer minimization: {e}") from e