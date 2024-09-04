# %%
import logging
import os
import shutil
import subprocess
from pathlib import Path
from tempfile import gettempdir

from typing import Optional, List, Tuple, Dict, Any, Union

from openeye import oechem
from openeye import oeomega
from openeye import oequacpac
from openmm.app import PDBFile

from pdbfixer import PDBFixer
from pymol import cmd

from openbabel import pybel

# %%
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    handlers=[
                        logging.StreamHandler()
                    ])
logger = logging.getLogger(__name__)


# %%
def openeye_fixer(oemol, explicit_H=True, protonation=False):
    to_edit = oechem.OEGraphMol(oemol)

    # oechem.OEDetermineConnectivity(to_edit)
    # oechem.OEFindRingAtomsAndBonds(to_edit)
    # oechem.OEPerceiveBondOrders(to_edit)
    # oechem.OEAssignImplicitHydrogens(to_edit)
    # oechem.OEAssignFormalCharges(to_edit)
    oechem.OEClearAromaticFlags(to_edit)

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
        if protonation == True:
            oequacpac.OEGetReasonableProtomer(to_edit)

    return oechem.OEGraphMol(to_edit)

def openeye_gen3d(smile) -> Path:
    oemol = oechem.OEMol()
    oechem.OESmilesToMol(oemol, smile)

    # Generate 3D coordinates
    builder = oeomega.OEConformerBuilder()
    ret_code = builder.Build(oemol)

    tmp_save = Path(gettempdir()) / "smiles_3d.sdf"
    ofs = oechem.oemolostream(tmp_save.as_posix())
    oechem.OEWriteMolecule(ofs, oechem.OEGraphMol(oemol))

    return tmp_save

def openeye_converter(input_file: Path, output_file: Path):
    ifs = oechem.oemolistream()
    ofs = oechem.oemolostream(str(output_file))
    if ifs.open(str(input_file)):
        for oemol in ifs.GetOEGraphMols():
            oechem.OEWriteMolecule(ofs, oemol)

class Plants_Docking:
    def __init__(self, protein_pdb: Path, input_ligands: Union[Path, str]):
        self.workdir: Path = Path(os.getcwd())
        self.protein_pdb: Path = protein_pdb
        self.crystal_sdf: Path = Path("Crystal.sdf")

        if isinstance(input_ligands, str):
            self.ligands_smiles: str = input_ligands
            self.input_mode = "SMILES"
        elif isinstance(input_ligands, Path):
            self.input_ligands: Path = input_ligands
            self.input_mode = "SDF"

        self.plants_env_var()

    def plants_env_var(self):
        os.environ['PATH'] = '/home/hitesit/Software/PLANTS:' + os.environ.get('PATH', '')

    def prepare_protein(self):
        protein_basename: str = self.protein_pdb.stem

        protein_tmp: Path = Path(gettempdir()) / f"{protein_basename}_tmp.pdb"
        protein_PREP: Path = self.workdir / f"{protein_basename}_PREP.mol2"

        fixer = PDBFixer(filename=str(self.protein_pdb))
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()

        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)

        PDBFile.writeFile(fixer.topology, fixer.positions, open(str(protein_tmp), 'w'))

        # Openbabel convert pdb to mol2
        prt = next(pybel.readfile("pdb", str(protein_tmp)))
        out = pybel.Outputfile("mol2", str(protein_PREP), overwrite=True)

        out.write(prt)
        out.close()

        self.protein_PREP: Path = protein_PREP
        return protein_PREP

    @staticmethod
    def fix_3d_mol(sdf_file: Path, TMP_basename: str) -> Path:
        TMP_fixed_sdf = Path(gettempdir()) / f"{TMP_basename}_fixed.mol2"

        ifs = oechem.oemolistream()
        ofs = oechem.oemolostream(TMP_fixed_sdf.as_posix())

        ifs.open(sdf_file.as_posix())
        for oemol in ifs.GetOEGraphMols():
            fixed_oemol = openeye_fixer(oemol)
            oechem.OEWriteMolecule(ofs, fixed_oemol)

        return TMP_fixed_sdf.resolve()

    def prepare_ligands(self):
        if self.input_mode == "SDF":
            logger.info("Preparing ligands from SDF")
            fixed_crystal: Path = self.fix_3d_mol(self.crystal_sdf, "crystal")
            fixed_ligands: Path = self.fix_3d_mol(self.input_ligands, "ligand")

            return fixed_ligands.resolve(), fixed_crystal.resolve()

        elif self.input_mode == "SMILES":
            logger.info("Preparing ligands from SMILES")
            fixed_crystal: Path = self.fix_3d_mol(self.crystal_sdf, "crystal")
            smiles_3d: Path = openeye_gen3d(self.ligands_smiles)
            fixed_ligands: Path = self.fix_3d_mol(smiles_3d, "ligand")

            return fixed_ligands.resolve(), fixed_crystal.resolve()

    def _define_binding_site(self):
        # temp convert the crystal.sdf to mol2
        tmp_crystal: Path = Path(gettempdir()) / "tmp.mol2"
        openeye_converter(self.crystal_sdf, tmp_crystal)

        plants_command = f"plants.64bit --mode bind {str(tmp_crystal)} {str(self.protein_PREP)}"
        print(plants_command)

        box_res = subprocess.run(plants_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, cwd=self.workdir)
        box_infos = box_res.stdout.split("\n")

        binding_site_center = box_infos[-4]
        binding_site_radius = box_infos[-3]

        return binding_site_center, binding_site_radius

    def write_conf(self, ligand_to_dock: str, n_confs: int):
        # Define the binding site
        center, radius = self._define_binding_site()

        config_path = self.workdir / "config.txt"
        config_str = f"""### PLANTS configuration file
# scoring function and search settings
scoring_function chemplp
search_speed speed1

# input
protein_file {self.protein_PREP}
ligand_file {ligand_to_dock}

# output
output_dir plants_raw
write_multi_mol2 1

# binding site definition
{center}
{radius}

# cluster algorithm
cluster_structures {n_confs}
cluster_rmsd 1.0
"""

        with open(config_path, "w") as config_file:
            config_file.write(config_str)
        return config_path

    def _retrieve_docked_ligands(self):
        plants_raw = self.workdir / "plants_raw"
        docked_ligand: Path = plants_raw / "docked_ligands.mol2"

        shutil.copy(docked_ligand, self.workdir / "docked_ligands.mol2")

    def run_plants_docking(self, n_confs: int):
        protein_PREP: Path = self.prepare_protein()
        fixed_ligands, fixed_crystal = self.prepare_ligands()
        config_path: Path = self.write_conf(str(fixed_ligands), n_confs)

        runner = f"plants.64bit --mode screen {str(config_path)}"
        subprocess.run(runner, shell=True, check=True, cwd=self.workdir)

        self._retrieve_docked_ligands()
        return


# %%
class Pymol_Docking:
    def __init__(self, protein_pdb: str, input_ligands: str):
        self.workdir: Path = Path(os.getcwd())
        self.protein_pdb: Path = Path(protein_pdb)
        self.crystal_sdf: Path = Path("Crystal.sdf")

        if os.path.exists(input_ligands):
            self.ligands_sdf: Path = Path(input_ligands)
            self.input_mode = "SDF"
        else:
            self.ligands_smiles: str = input_ligands
            self.input_mode = "SMILES"

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
        if self.input_mode == "SDF":
            logger.info("Preparing ligands from SDF")
            fixed_crystal: Path = self.fix_3d_mol(self.crystal_sdf, "crystal")
            fixed_ligands: Path = self.fix_3d_mol(self.ligands_sdf, "ligand")

            return fixed_ligands.resolve(), fixed_crystal.resolve()

        elif self.input_mode == "SMILES":
            logger.info("Preparing ligands from SMILES")
            fixed_crystal: Path = self.fix_3d_mol(self.crystal_sdf, "crystal")
            smiles_3d: Path = openeye_gen3d(self.ligands_smiles)
            fixed_ligands: Path = self.fix_3d_mol(smiles_3d, "ligand")

            return fixed_ligands.resolve(), fixed_crystal.resolve()

    def run_smina_docking(self, mode: str, docking_basename: str) -> Path:
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
            assert self.input_mode == "SDF", "Minimization only works with SDF input"
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

# %%
@cmd.extend
def off_site_docking(protein_selection, ligand_smiles, outname:str):

    protein_name = cmd.get_object_list(protein_selection)[0]
    to_save_protein: Path = Path(gettempdir()) / f"{protein_name}.pdb"

    cmd.save(to_save_protein.as_posix(), protein_selection)

    pymol_docking = Pymol_Docking(to_save_protein.as_posix(), str(ligand_smiles))
    docked_sdf: Path = pymol_docking.run_smina_docking("Dock", outname)

    cmd.load(str(docked_sdf))

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
    assert_organic(ligand_name)

    to_save_protein = Path(gettempdir()) / f"{protein_name}.pdb"
    to_save_ligand = Path(gettempdir()) / f"{ligand_name}.sdf"

    cmd.save(str(to_save_protein), protein_selection)
    cmd.save(str(to_save_ligand), ligand_selection)

    pymol_docking = Pymol_Docking(str(to_save_protein), str(to_save_ligand))
    docked_sdf: Path = pymol_docking.run_smina_docking(mode, outname)

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
        self.choose_docking_modes()

        # Choose modality
        self.choose_setting()

        # Connect modality chooser signal to a method that clears and repopulates the selectors
        self.mode_chooser_2.currentIndexChanged.connect(self.clear_and_repopulate_selectors)

        # Connect the accepted signal to a method that checks the modality
        dialog.accepted.connect(self.on_dialog_accepted)

        dialog.exec_()

    def clear_and_repopulate_selectors(self):
        """Clear and repopulate the contents of the protein and ligand selectors."""
        self.protein_chooser_1.clear()
        self.protein_chooser_2.clear()
        self.ligand_chooser_1.clear()
        self.mode_chooser.clear()

        # Repopulate the selectors with available options
        self.populate_ligand_select_list()
        self.choose_docking_modes()

    def on_dialog_accepted(self):
        modality = self.mode_chooser_2.currentText().strip()

        if modality == "In-Site":
            self.in_site_wrapper()
            logger.info("In-Site docking selected")
        elif modality == "Off-Site":
            self.off_site_wrapper()
            logger.info("Off-Site docking selected")
        else:
            logger.warning("No valid modality selected")

    def in_site_wrapper(self):
        s1 = self.protein_chooser_1.currentText()
        s2 = self.ligand_chooser_1.currentText()
        mode = self.mode_chooser.currentText()
        outname = self.output_chooser.toPlainText()
        logger.info("Running on-site docking...")
        on_site_docking(s1, s2, mode, outname)

    def off_site_wrapper(self):
        s1 = self.protein_chooser_2.currentText()
        s2 = self.smile_chooser_2.toPlainText()
        outname = self.output_chooser_2.toPlainText()
        logger.info("Running off-site docking...")
        off_site_docking(s1, str(s2), outname)

    def choose_setting(self):
        possible_choices = ["In-Site", "Off-Site"]
        self.mode_chooser_2.addItems(possible_choices)

    def choose_docking_modes(self):
        possible_choices = ["Minimize", "Dock"]
        self.mode_chooser.addItems(possible_choices)

    def populate_ligand_select_list(self):
        loaded_objects = _get_select_list()

        self.protein_chooser_1.clear()
        self.protein_chooser_2.clear()
        self.ligand_chooser_1.clear()

        self.protein_chooser_1.addItems(loaded_objects)
        self.protein_chooser_2.addItems(loaded_objects)
        self.ligand_chooser_1.addItems(loaded_objects)

    def setupUi(self, Form):
        from PyQt5 import QtCore, QtWidgets
        Form.setObjectName("Form")
        Form.resize(662, 302)
        self.buttonBox = QtWidgets.QDialogButtonBox(Form)
        self.buttonBox.setGeometry(QtCore.QRect(440, 270, 193, 28))
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayoutWidget = QtWidgets.QWidget(Form)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(170, 20, 271, 91))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.protein_chooser_1 = QtWidgets.QComboBox(self.verticalLayoutWidget)
        self.protein_chooser_1.setObjectName("protein_chooser_1")
        self.verticalLayout.addWidget(self.protein_chooser_1)
        self.ligand_chooser_1 = QtWidgets.QComboBox(self.verticalLayoutWidget)
        self.ligand_chooser_1.setObjectName("ligand_chooser_1")
        self.verticalLayout.addWidget(self.ligand_chooser_1)
        self.mode_chooser = QtWidgets.QComboBox(self.verticalLayoutWidget)
        self.mode_chooser.setObjectName("mode_chooser")
        self.verticalLayout.addWidget(self.mode_chooser)
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(Form)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(170, 140, 271, 118))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.protein_chooser_2 = QtWidgets.QComboBox(self.verticalLayoutWidget_2)
        self.protein_chooser_2.setObjectName("protein_chooser_2")
        self.verticalLayout_2.addWidget(self.protein_chooser_2)
        self.smile_chooser_2 = QtWidgets.QPlainTextEdit(self.verticalLayoutWidget_2)
        self.smile_chooser_2.setObjectName("smile_chooser_2")
        self.verticalLayout_2.addWidget(self.smile_chooser_2)
        self.output_chooser = QtWidgets.QTextEdit(Form)
        self.output_chooser.setGeometry(QtCore.QRect(460, 20, 151, 91))
        self.output_chooser.setObjectName("output_chooser")
        self.mode_chooser_2 = QtWidgets.QComboBox(Form)
        self.mode_chooser_2.setGeometry(QtCore.QRect(40, 140, 73, 22))
        self.mode_chooser_2.setObjectName("mode_chooser_2")
        self.output_chooser_2 = QtWidgets.QTextEdit(Form)
        self.output_chooser_2.setGeometry(QtCore.QRect(460, 140, 151, 111))
        self.output_chooser_2.setObjectName("output_chooser_2")

        self.buttonBox.accepted.connect(Form.accept)
        self.buttonBox.rejected.connect(Form.reject)


def __init__(self):
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('Pymol Docking', Pymol_Docking_GUI)
        return
    except Exception as e:
        print(e)