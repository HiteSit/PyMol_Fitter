import os
from pathlib import Path
from tempfile import gettempdir
import logging
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi
from pymol import cmd

# Import the client for communicating with the Docker server
from .client import PyMOLDockingClient

# Create the client instance
docker_client = PyMOLDockingClient("http://localhost:5000")

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
    """
    Perform off-site docking using PyMOL protein selection and a SMILES string for the ligand.
    
    This function uses the Docker server to perform the docking.
    
    Parameters:
    - protein_selection (str): The PyMOL selection string for the protein.
    - ligand_smiles (str): The SMILES string for the ligand.
    - outname (str): The base name for the output file.
    """
    # Check if the Docker server is running
    if not docker_client.check_health():
        print("Error: Docker server is not running. Please start the Docker server first.")
        return
    
    # Save the protein to a temporary file
    protein_name = cmd.get_object_list(protein_selection)[0]
    to_save_protein: Path = Path(gettempdir()) / f"{protein_name}.pdb"
    cmd.save(to_save_protein.as_posix(), protein_selection)

    try:
        # Run the docking process
        print(f"Running docking with {protein_name} and SMILES: {ligand_smiles}")
        results = docker_client.dock(
            protein_file=to_save_protein,
            ligand=ligand_smiles,
            is_smiles=True,
            dock_mode="Dock",
            output_name=outname,
            output_dir=gettempdir()
        )
        
        # Load the results into PyMOL
        cmd.load(str(results["docked_ligand"]), outname)
        cmd.load(str(results["prepared_protein"]), f"{outname}_protein")
        print(f"Docking completed successfully! Results loaded as {outname} and {outname}_protein")
        
    except Exception as e:
        print(f"Error during docking: {e}")

@cmd.extend
def on_site_docking(protein_selection, ligand_selection, mode, outname: str, minimization_flag: bool = False):
    """
    Perform on-site docking using PyMOL selections for protein and ligand, and smina for docking or minimization.

    This function saves the selected protein and ligand as temporary files, runs the docking or minimization process
    using the Docker server, and then loads the docked ligand structure back into PyMOL.

    Parameters:
    - protein_selection (str): The PyMOL selection string for the protein.
    - ligand_selection (str): The PyMOL selection string for the ligand.
    - mode (str): The operation mode, either "Minimize" for energy minimization or "Dock" for docking.
    - outname (str): The base name for the output file. The docked structure will be saved as "{outname}.sdf"
                     and loaded into PyMOL with the same name.
    - minimization_flag (bool): Whether to perform minimization after docking, default is False.

    Raises:
    - AssertionError: If the mode is not one of the expected values ("Minimize" or "Dock").

    Returns:
    None. The result of the docking or minimization is loaded into PyMOL.
    """
    # Check if the Docker server is running
    if not docker_client.check_health():
        print("Error: Docker server is not running. Please start the Docker server first.")
        return
    
    # Assert the mode
    assert mode in ["Minimize", "Dock"]

    protein_name = cmd.get_object_list(protein_selection)[0]
    ligand_name = cmd.get_object_list(ligand_selection)[0]
    assert_organic(ligand_name)

    to_save_protein = Path(gettempdir()) / f"{protein_name}.pdb"
    to_save_ligand = Path(gettempdir()) / f"{ligand_name}.sdf"

    cmd.save(str(to_save_protein), protein_selection)
    cmd.save(str(to_save_ligand), ligand_selection)

    try:
        if minimization_flag:
            # Run docking and minimization
            print(f"Running docking and minimization with {protein_name} and {ligand_name}")
            results = docker_client.dock_and_minimize(
                protein_file=to_save_protein,
                ligand=to_save_ligand,
                is_smiles=False,
                dock_mode=mode,
                output_name=outname,
                output_dir=gettempdir()
            )
            
            # Load the results into PyMOL
            cmd.load(str(results["docked_ligand"]), outname)
            cmd.load(str(results["prepared_protein"]), f"{outname}_protein")
            cmd.load(str(results["minimized_complex"]), f"{outname}_complex")
            print(f"Docking and minimization completed successfully! Results loaded as {outname}, {outname}_protein, and {outname}_complex")
        else:
            # Run docking only
            print(f"Running docking with {protein_name} and {ligand_name}")
            results = docker_client.dock(
                protein_file=to_save_protein,
                ligand=to_save_ligand,
                is_smiles=False,
                dock_mode=mode,
                output_name=outname,
                output_dir=gettempdir()
            )
            
            # Load the results into PyMOL
            cmd.load(str(results["docked_ligand"]), outname)
            cmd.load(str(results["prepared_protein"]), f"{outname}_protein")
            print(f"Docking completed successfully! Results loaded as {outname} and {outname}_protein")
            
    except Exception as e:
        print(f"Error during docking: {e}")

@cmd.extend
def docker_minimize_complex(protein_selection, ligand_selection, outname: str):
    """
    Perform minimization of a protein-ligand complex using the Docker server.
    
    Parameters:
    - protein_selection (str): The PyMOL selection string for the protein.
    - ligand_selection (str): The PyMOL selection string for the ligand.
    - outname (str): The base name for the output file.
    """
    # Check if the Docker server is running
    if not docker_client.check_health():
        print("Error: Docker server is not running. Please start the Docker server first.")
        return
    
    protein_name = cmd.get_object_list(protein_selection)[0]
    ligand_name = cmd.get_object_list(ligand_selection)[0]
    assert_organic(ligand_name)

    to_save_protein = Path(gettempdir()) / f"{protein_name}.pdb"
    to_save_ligand = Path(gettempdir()) / f"{ligand_name}.sdf"

    cmd.save(str(to_save_protein), protein_selection)
    cmd.save(str(to_save_ligand), ligand_selection)

    try:
        # Run minimization
        print(f"Running minimization with {protein_name} and {ligand_name}")
        results = docker_client.minimize_complex(
            protein_file=to_save_protein,
            ligand_file=to_save_ligand,
            output_name=outname,
            output_dir=gettempdir()
        )
        
        # Load the results into PyMOL
        cmd.load(str(results["minimized_complex"]), f"{outname}_complex")
        print(f"Minimization completed successfully! Result loaded as {outname}_complex")
        
    except Exception as e:
        print(f"Error during minimization: {e}")

# Define the dialog class that inherits from QDialog
class PymolDockingDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(PymolDockingDialog, self).__init__(parent)
        
        # Load the UI file directly
        uifile = os.path.join(os.path.dirname(__file__), 'GUI.ui')
        loadUi(uifile, self)
        
        # Check if the Docker server is running
        if not docker_client.check_health():
            self.setWindowTitle("PyMOL Docking - Docker Server Not Running")
            # Create a warning message
            warning_label = QtWidgets.QLabel("Warning: Docker server is not running.\nPlease start the Docker server first.", self)
            warning_label.setStyleSheet("color: red; font-weight: bold;")
            # Add the warning to the layout
            self.layout().addWidget(warning_label)
        
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
        # Check if the Docker server is running
        if not docker_client.check_health():
            print("Error: Docker server is not running. Please start the Docker server first.")
            return
            
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
        """Wrapper for in-site docking."""
        s1 = self.protein_chooser_1.currentText()
        s2 = self.ligand_chooser_1.currentText()
        mode = self.mode_chooser.currentText()
        outname = self.output_chooser.toPlainText()
        logger.info("Running on-site docking...")
        
        minimization_flag = self.minimizer.isChecked()
        
        on_site_docking(s1, s2, mode, outname, minimization_flag)

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