# %%
import os
import shutil
import glob
import logging
import pathlib
from pathlib import Path
from tempfile import gettempdir, NamedTemporaryFile
from typing import List, Dict, Union

import subprocess

from openmm.app import PDBFile
from pdbfixer import PDBFixer

from pymol import cmd

from openeye import oechem
from openeye import oequacpac

# %%
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    handlers=[
                        logging.StreamHandler()
                    ])
logger = logging.getLogger(__name__)


# %%
def openeye_fixer(oemol, explicit_H=True):
    to_edit = oechem.OEGraphMol(oemol)

    oechem.OEDetermineConnectivity(to_edit)
    oechem.OEFindRingAtomsAndBonds(to_edit)
    oechem.OEPerceiveBondOrders(to_edit)
    oechem.OEAssignImplicitHydrogens(to_edit)
    oechem.OEAssignFormalCharges(to_edit)

    oechem.OEClearAromaticFlags(to_edit)

    oechem.OEFindRingAtomsAndBonds(to_edit)
    oechem.OEAssignAromaticFlags(to_edit)
    for bond in to_edit.GetBonds():
        if bond.IsAromatic():
            bond.SetIntType(5)
        elif bond.GetOrder() != 0:
            bond.SetIntType(bond.GetOrder())
        else:
            bond.SetIntType(1)

    oechem.OEKekulize(to_edit)

    if explicit_H:
        oechem.OEAddExplicitHydrogens(to_edit)
        oequacpac.OEGetReasonableProtomer(to_edit)

    return oechem.OEGraphMol(to_edit)


# %%
class Pymol_Docking:
    def __init__(self, protein_pdb: str, ligands_sdf: str):
        self.workdir: Path = Path(os.getcwd())
        self.protein_pdb: Path = Path(protein_pdb)
        self.crystal_sdf: Path = Path("Crystal.sdf")

        self.ligands_sdf: Path = Path(ligands_sdf)

    def prepare_protein(self) -> Path:
        protein_basename: str = self.protein_pdb.stem
        protein_PREP: Path = self.workdir / f"{protein_basename}_PREP.pdb"

        fixer = PDBFixer(filename=str(self.protein_pdb))
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()

        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)

        PDBFile.writeFile(fixer.topology, fixer.positions, open(protein_PREP.as_posix(), 'w'))
        return protein_PREP

    @staticmethod
    def fix_3d_mol(sdf_file: Path, TMP_basename: str) -> Path:
        TMP_fixed_sdf = Path(gettempdir()) / f"{TMP_basename}_fixed.sdf"

        ifs = oechem.oemolistream()
        ofs = oechem.oemolostream(TMP_fixed_sdf.as_posix())

        ifs.open(sdf_file.as_posix())
        for oemol in ifs.GetOEGraphMols():
            fixed_oemol = openeye_fixer(oemol)
            oechem.OEWriteMolecule(ofs, fixed_oemol)

        return TMP_fixed_sdf.resolve()

    def prepare_ligands(self):
        fixed_crystal: Path = self.fix_3d_mol(self.crystal_sdf, "crystal")
        fixed_ligands: Path = self.fix_3d_mol(self.ligands_sdf, "ligand")

        return fixed_ligands.resolve(), fixed_crystal.resolve()

    def run_docking(self, mode: str, docking_basename: str) -> Path:
        # Fix the protein
        print("Running protein preparation")
        protein_PKA: Path = self.prepare_protein()

        # Fix the ligands
        fixed_ligands, fixed_crystal = self.prepare_ligands()

        # Define the output
        smina_output: Path = self.workdir / f"{docking_basename}.sdf"
        smina_log = smina_output.with_suffix(".log")

        if mode == "Dock":
            smina_lst = [
                "smina",
                "-r", protein_PKA.as_posix(),
                "-l", fixed_ligands.as_posix(),
                "--autobox_ligand", fixed_crystal.as_posix(),
                "-o", smina_output,
                "--exhaustiveness", "32"
            ]

            print("Running docking")
            with open(smina_log, "w") as log_file:
                subprocess.run(smina_lst, check=True, stdout=log_file, stderr=log_file)

        elif mode == "Minimize":
            smina_lst = [
                "smina",
                "-r", protein_PKA.as_posix(),
                "-l", fixed_ligands.as_posix(),
                "--autobox_ligand", fixed_crystal.as_posix(),
                "-o", smina_output,
                "--minimize", "--minimize_iters", "10"
            ]
            print("Running docking")
            with open(smina_log, "w") as log_file:
                subprocess.run(smina_lst, check=True, stdout=log_file, stderr=log_file)

        return smina_output

# # Dry Example
# pymol_docking = Pymol_Docking("./LAC3.pdb", "Ligand.sdf")
# pymol_docking.run_docking("Minimize", "XXX")

# %%
@cmd.extend
def on_site_docking(protein_selection, ligand_selection, mode, outname: str):
    """
    Perform on-site docking using PyMOL selections for protein and ligand, and smina for docking or minimization.

    This function saves the selected protein and ligand as temporary files, runs the docking or minimization process
    using smina, and then loads the docked ligand structure back into PyMOL.

    Parameters:
    - protein_selection (str): The PyMOL selection string for the protein.
    - ligand_selection (str): The PyMOL selection string for the ligand.
    - mode (str): The operation mode, either "Minimize" for energy minimization or "Dock" for docking.
    - outname (str): The base name for the output file. The docked structure will be saved as "{outname}.sdf"
                     and loaded into PyMOL with the same name.

    Raises:
    - AssertionError: If the mode is not one of the expected values ("Minimize" or "Dock").

    Returns:
    None. The result of the docking or minimization is loaded into PyMOL.
    """
    # Assert the mode
    assert mode in ["Minimize", "Dock"]

    protein_name = cmd.get_object_list(protein_selection)[0]
    ligand_name = cmd.get_object_list(ligand_selection)[0]

    to_save_protein = Path(gettempdir()) / f"{protein_name}.pdb"
    to_save_ligand = Path(gettempdir()) / f"{ligand_name}.sdf"

    cmd.save(str(to_save_protein), protein_selection)
    cmd.save(str(to_save_ligand), ligand_selection)

    pymol_docking = Pymol_Docking(str(to_save_protein), str(to_save_ligand))
    docked_sdf: Path = pymol_docking.run_docking(mode, outname)

    cmd.load(str(docked_sdf))
# %%