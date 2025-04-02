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

    # Define crystal file path (using "./Crystal.sdf" as specified)
    crystal_file = Path("./Crystal.sdf")

    try:
        # Run the docking process using the consolidated dock_minimize function
        print(f"Running docking with {protein_name} and SMILES: {ligand_smiles}")
        results = docker_client.dock_minimize(
            protein_file=to_save_protein,
            ligand=ligand_smiles,
            crystal_file=crystal_file,
            is_smiles=True,
            minimize=False,
            dock_mode="Dock",
            output_name=outname,
            output_dir="."
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
    - mode (str): The operation mode, either "Minimize" or "Dock".
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

    # Define crystal file path (using "./Crystal.sdf" as specified)
    crystal_file = Path("./Crystal.sdf")

    try:
        # Run docking using the consolidated dock_minimize function
        print(f"Running {'docking with minimization' if minimization_flag else 'docking'} with {protein_name} and {ligand_name}")
        results = docker_client.dock_minimize(
            protein_file=to_save_protein,
            ligand=to_save_ligand,
            crystal_file=crystal_file,
            is_smiles=False,
            minimize=minimization_flag,
            dock_mode=mode,
            output_name=outname,
            output_dir="."
        )
        
        # Load the results into PyMOL
        cmd.load(str(results["docked_ligand"]), outname)
        cmd.load(str(results["prepared_protein"]), f"{outname}_protein")
        
        # If minimization was performed, load the minimized complex
        if minimization_flag and "minimized_complex" in results:
            cmd.load(str(results["minimized_complex"]), f"{outname}_complex")
            print(f"Docking and minimization completed successfully! Results loaded as {outname}, {outname}_protein, and {outname}_complex")
        else:
            print(f"Docking completed successfully! Results loaded as {outname} and {outname}_protein")
            
    except Exception as e:
        print(f"Error during docking: {e}")

@cmd.extend
def md_minimization(protein_selection, ligand_selection, outname: str = "minimized_complex"):
    """
    Perform molecular dynamics minimization on a protein-ligand complex.
    
    This function saves the selected protein and ligand as temporary files, 
    runs the minimization process using the Docker server, and then loads
    the minimized complex structure back into PyMOL.
    
    Parameters:
    - protein_selection (str): The PyMOL selection string for the protein.
    - ligand_selection (str): The PyMOL selection string for the ligand.
    - outname (str): The name for the output complex file, default is "minimized_complex".
    
    Returns:
    None. The result of the minimization is loaded into PyMOL.
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
        # Run minimization using the inplace_minimization function
        print(f"Running MD minimization with {protein_name} and {ligand_name}")
        results = docker_client.inplace_minimization(
            protein_file=to_save_protein,
            ligand=to_save_ligand,
            output_dir="."
        )
        
        # Load the results into PyMOL
        cmd.load(str(results["minimized_complex"]), outname)
        print(f"MD minimization completed successfully! Result loaded as {outname}")
        
    except Exception as e:
        print(f"Error during MD minimization: {e}")

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
            # Add the warning to the tab's layout
            if self.tab.layout() is None:
                # If tab has no layout, create one
                layout = QtWidgets.QVBoxLayout(self.tab)
                layout.addWidget(warning_label)
            else:
                self.tab.layout().addWidget(warning_label)
        
        # Populate the UI with initial data
        self.populate_ligand_select_list()
        self.choose_docking_modes()
        self.choose_setting()
        
        # Connect signals to slots
        self.mode_chooser_2.currentIndexChanged.connect(self.clear_and_repopulate_selectors)
        self.mode_chooser.currentIndexChanged.connect(self.on_mode_changed)
        self.buttonBox.accepted.connect(self.on_dialog_accepted)
        # Connect tab widget change signal
        self.tabWidget.currentChanged.connect(self.on_tab_changed)
        # Connect buttonBox_MD for the second tab
        self.buttonBox_MD.accepted.connect(self.on_md_dialog_accepted)
        # Note: buttonBox is assumed to be in GUI.ui and connected to accept/reject in the UI file
        
        # Initially set minimizer visibility based on the default mode selection
        self.on_mode_changed(self.mode_chooser.currentIndex())
    
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

    def populate_md_select_list(self):
        """Populate the MD tab protein and ligand selectors with available options."""
        loaded_objects = self._get_select_list()
        self.protein_chooser_MD.clear()
        self.ligand_chooser_MD.clear()
        self.protein_chooser_MD.addItems(loaded_objects)
        self.ligand_chooser_MD.addItems(loaded_objects)

    def choose_docking_modes(self):
        """Populate the docking mode chooser."""
        possible_choices = ["Minimize", "Dock"]
        self.mode_chooser.addItems(possible_choices)
        
    def choose_setting(self):
        """Populate the modality chooser and set initial visibility."""
        possible_choices = ["In-Site", "Off-Site"]
        self.mode_chooser_2.addItems(possible_choices)
        
        # Set initial visibility based on default selection
        modality = self.mode_chooser_2.currentText().strip()
        
        # Default visibility setup for In-Site mode
        in_site_visible = (modality == "In-Site")
        # These are the widgets for In-Site docking
        self.protein_chooser_1.setVisible(in_site_visible)
        self.ligand_chooser_1.setVisible(in_site_visible) 
        self.mode_chooser.setVisible(in_site_visible)
        self.output_chooser.setVisible(in_site_visible)
        # Don't set minimizer visibility here - it's controlled by mode_chooser
        
        # Default visibility setup for Off-Site mode
        off_site_visible = (modality == "Off-Site")
        # These are the widgets for Off-Site docking
        self.protein_chooser_2.setVisible(off_site_visible)
        self.smile_chooser_2.setVisible(off_site_visible)
        self.output_chooser_2.setVisible(off_site_visible)
        
        # You may also want to set visibility for the layout containers
        self.verticalLayoutWidget.setVisible(in_site_visible)
        self.verticalLayoutWidget_2.setVisible(off_site_visible)

    def clear_and_repopulate_selectors(self):
        """Clear and repopulate the selectors when modality changes, and toggle visibility."""
        # Clear and repopulate as before
        self.protein_chooser_1.clear()
        self.protein_chooser_2.clear()
        self.ligand_chooser_1.clear()
        self.mode_chooser.clear()
        self.populate_ligand_select_list()
        self.choose_docking_modes()
        
        # Toggle visibility based on modality
        modality = self.mode_chooser_2.currentText().strip()
        
        # For In-Site mode
        in_site_visible = (modality == "In-Site")
        self.protein_chooser_1.setVisible(in_site_visible)
        self.ligand_chooser_1.setVisible(in_site_visible)
        self.mode_chooser.setVisible(in_site_visible)
        self.output_chooser.setVisible(in_site_visible)
        # Don't set minimizer visibility here - it's controlled by mode_chooser
        self.verticalLayoutWidget.setVisible(in_site_visible)
        
        # Update the minimizer visibility based on both the modality and the selected mode
        self.on_mode_changed(self.mode_chooser.currentIndex())
        
        # For Off-Site mode
        off_site_visible = (modality == "Off-Site")
        self.protein_chooser_2.setVisible(off_site_visible)
        self.smile_chooser_2.setVisible(off_site_visible)
        self.output_chooser_2.setVisible(off_site_visible)
        self.verticalLayoutWidget_2.setVisible(off_site_visible)
        
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
        
    def on_tab_changed(self, index):
        """Handle tab change events."""
        if index == 0:
            # First tab (original docking interface)
            logger.info("Switched to docking tab")
            # You could refresh the UI here if needed
            self.populate_ligand_select_list()
        elif index == 1:
            # Second tab (MD minimization)
            logger.info("Switched to MD minimization tab")
            # Initialize MD tab functionality
            self.populate_md_select_list()
        elif index == 2:
            # Third tab (new functionality)
            logger.info("Switched to second tab")
            # Initialize any second tab functionality here
            
    def on_mode_changed(self, index):
        """Toggle minimizer visibility based on the selected docking mode."""
        mode = self.mode_chooser.currentText()
        modality = self.mode_chooser_2.currentText().strip()
        
        # Only show minimizer checkbox when in "In-Site" modality AND "Minimize" mode
        minimizer_visible = (modality == "In-Site" and mode == "Minimize")
        self.minimizer.setVisible(minimizer_visible)
        logger.info(f"Mode changed to {mode}, minimizer visibility set to {minimizer_visible}")

    def md_minimizer_wrapper(self):
        """Wrapper for MD minimization."""
        protein = self.protein_chooser_MD.currentText()
        ligand = self.ligand_chooser_MD.currentText()
        outname = "minimized_complex"  # Default name, can be edited later
        logger.info("Running MD minimization...")
        md_minimization(protein, ligand, outname)
    
    def on_md_dialog_accepted(self):
        """Handle dialog acceptance for MD tab."""
        # Check if the Docker server is running
        if not docker_client.check_health():
            print("Error: Docker server is not running. Please start the Docker server first.")
            return
        
        # Call the MD minimizer wrapper
        self.md_minimizer_wrapper()
        logger.info("MD minimization selected")

# Plugin initialization for PyMOL
def __init_plugin__(app=None):
    """Add the plugin to PyMOL's menu."""
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Pymol Docking', run_plugin_gui)

def run_plugin_gui():
    """Create and show the dialog when the menu item is clicked."""
    dialog = PymolDockingDialog()
    dialog.exec_()  # Show the dialog modally 