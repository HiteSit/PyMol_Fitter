import os
from pathlib import Path
from tempfile import gettempdir
import logging
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi
from pymol import cmd

logger = logging.getLogger(__name__)

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

@cmd.extend
def off_site_docking(protein_selection, ligand_smiles, outname:str):
    from .Docking_Engine import Pymol_Docking

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
    from .Docking_Engine import Pymol_Docking
    
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

# Define the dialog class that inherits from QDialog
class PymolDockingDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(PymolDockingDialog, self).__init__(parent)
        
        # Load the UI file directly
        uifile = os.path.join(os.path.dirname(__file__), 'GUI.ui')
        loadUi(uifile, self)
        
        # Populate the UI with initial data
        self.populate_ligand_select_list()
        self.choose_docking_modes()
        self.choose_setting()
        
        # Connect signals to slots
        self.mode_chooser_2.currentIndexChanged.connect(self.clear_and_repopulate_selectors)
        self.buttonBox.accepted.connect(self.on_dialog_accepted)
        # Note: buttonBox is assumed to be in GUI.ui and connected to accept/reject in the UI file
    
    @staticmethod
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

    def populate_ligand_select_list(self):
        """Populate the protein and ligand selectors with available options."""
        loaded_objects = self._get_select_list()
        self.protein_chooser_1.clear()
        self.protein_chooser_2.clear()
        self.ligand_chooser_1.clear()
        self.protein_chooser_1.addItems(loaded_objects)
        self.protein_chooser_2.addItems(loaded_objects)
        self.ligand_chooser_1.addItems(loaded_objects)

    def choose_docking_modes(self):
        """Populate the docking mode chooser."""
        possible_choices = ["Minimize", "Dock"]
        self.mode_chooser.addItems(possible_choices)

    def choose_setting(self):
        """Populate the modality chooser."""
        possible_choices = ["In-Site", "Off-Site"]
        self.mode_chooser_2.addItems(possible_choices)

    def clear_and_repopulate_selectors(self):
        """Clear and repopulate the selectors when modality changes."""
        self.protein_chooser_1.clear()
        self.protein_chooser_2.clear()
        self.ligand_chooser_1.clear()
        self.mode_chooser.clear()
        self.populate_ligand_select_list()
        self.choose_docking_modes()

    def on_dialog_accepted(self):
        """Handle dialog acceptance based on selected modality."""
        modality = self.mode_chooser_2.currentText().strip()
        if modality == "In-Site":
            self.in_site_wrapper()
            logger.info("In-Site docking selected")  # Assuming logger is defined
        elif modality == "Off-Site":
            self.off_site_wrapper()
            logger.info("Off-Site docking selected")
        else:
            logger.warning("No valid modality selected")

    def in_site_wrapper(self):
        """Wrapper for in-site docking."""
        s1 = self.protein_chooser_1.currentText()
        s2 = self.ligand_chooser_1.currentText()
        mode = self.mode_chooser.currentText()
        outname = self.output_chooser.toPlainText()
        logger.info("Running on-site docking...")
        on_site_docking(s1, s2, mode, outname)

    def off_site_wrapper(self):
        """Wrapper for off-site docking."""
        s1 = self.protein_chooser_2.currentText()
        s2 = self.smile_chooser_2.toPlainText()
        outname = self.output_chooser_2.toPlainText()
        logger.info("Running off-site docking...")
        off_site_docking(s1, str(s2), outname)

# Plugin initialization for PyMOL
def __init_plugin__(app=None):
    """Add the plugin to PyMOL's menu."""
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Pymol Docking', run_plugin_gui)

def run_plugin_gui():
    """Create and show the dialog when the menu item is clicked."""
    dialog = PymolDockingDialog()
    dialog.exec_()  # Show the dialog modally