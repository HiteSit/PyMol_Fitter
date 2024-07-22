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

    def assert_organic(selection):
        """
        Assert that the given PyMOL selection consists of only organic molecules.

        Parameters:
        - selection (str): The PyMOL selection string to check.

        Raises:
        - AssertionError: If the selection contains non-organic molecules.
        """
        # Create a temporary selection for organic molecules
        cmd.select("organic_check", f"{selection} and organic")

        # Get the total number of atoms in the original selection
        total_atoms = cmd.count_atoms(selection)

        # Get the number of atoms in the organic selection
        organic_atoms = cmd.count_atoms("organic_check")

        # Check if the number of organic atoms is equal to the total number of atoms
        assert organic_atoms == total_atoms, "Selection contains non-organic molecules."
        logger.info(f"Selection {selection} contains only organic molecules.")

        # Delete the temporary selection
        cmd.delete("organic_check")

    # Assert the mode
    assert mode in ["Minimize", "Dock"]

    protein_name = cmd.get_object_list(protein_selection)[0]
    ligand_name = cmd.get_object_list(ligand_selection)[0]
    assert_organic(ligand_name)

    to_save_protein = Path(gettempdir()) / f"{protein_name}.pdb"
    to_save_ligand = Path(gettempdir()) / f"{ligand_name}.sdf"

    cmd.save(str(to_save_protein), protein_selection)
    cmd.save(str(to_save_ligand), ligand_selection)

    pymol_docking = Pymol_Docking(str(to_save_protein), str(to_save_ligand))
    docked_sdf: Path = pymol_docking.run_docking(mode, outname)

    cmd.load(str(docked_sdf))
# %%

#################################################################################
########################### Start of pymol plugin code ##########################
#################################################################################

def _get_select_list():
    '''
    Get either a list of object names, or a list of chain selections
    '''
    loaded_objects = [name for name in cmd.get_names('all', 1) if '_cluster_' not in name]

    # if single object, try chain selections
    if len(loaded_objects) == 1:
        chains = cmd.get_chains(loaded_objects[0])
        if len(chains) > 1:
            loaded_objects = ['{} & chain {}'.format(loaded_objects[0], chain) for chain in chains]

    return loaded_objects


class Pymol_Docking_GUI(object):
    ''' Qt version of the Plugin GUI '''
    def __init__(self):
        from pymol.Qt import QtWidgets
        dialog = QtWidgets.QDialog()
        self.setupUi(dialog)
        self.populate_ligand_select_list()
        self.choose_modes()
        dialog.accepted.connect(self.accept)
        dialog.exec_()

    def accept(self):
        s1 = self.comboBox_SX.currentText()
        s2 = self.comboBox_DX.currentText()
        mode = self.mode_chooser.currentText()
        outname = self.output_chooser.toPlainText()
        on_site_docking(s1, s2, mode, outname)

    def choose_modes(self):
        possible_choices = ["Minimize", "Dock"]
        self.mode_chooser.addItems(possible_choices)

    def populate_ligand_select_list(self):
        loaded_objects = _get_select_list()

        self.comboBox_SX.clear()
        self.comboBox_DX.clear()

        self.comboBox_SX.addItems(loaded_objects)
        self.comboBox_DX.addItems(loaded_objects)

        if len(loaded_objects) > 1:
            self.comboBox_DX.setCurrentIndex(1)

    def setupUi(self, Dialog):
        from PyQt5 import QtCore, QtGui, QtWidgets
        Dialog.setObjectName("Dialog")
        Dialog.resize(662, 229)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(410, 140, 193, 28))
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayoutWidget = QtWidgets.QWidget(Dialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(40, 20, 251, 191))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.comboBox_SX = QtWidgets.QComboBox(self.verticalLayoutWidget)
        self.comboBox_SX.setObjectName("comboBox_SX")
        self.verticalLayout.addWidget(self.comboBox_SX)
        self.comboBox_DX = QtWidgets.QComboBox(self.verticalLayoutWidget)
        self.comboBox_DX.setObjectName("comboBox_DX")
        self.verticalLayout.addWidget(self.comboBox_DX)
        self.horizontalLayoutWidget = QtWidgets.QWidget(Dialog)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(370, 20, 271, 89))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.mode_chooser = QtWidgets.QComboBox(self.horizontalLayoutWidget)
        self.mode_chooser.setObjectName("mode_chooser")
        self.horizontalLayout.addWidget(self.mode_chooser)
        self.output_chooser = QtWidgets.QPlainTextEdit(self.horizontalLayoutWidget)
        self.output_chooser.setObjectName("output_chooser")
        self.horizontalLayout.addWidget(self.output_chooser)

        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)


def __init__(self):
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('Pymol Docking', Pymol_Docking_GUI)
        return
    except Exception as e:
        print(e)