import os
from pathlib import Path
from tempfile import gettempdir
import logging
import re
from pymol.Qt import QtWidgets
from pymol.Qt.utils import loadUi
from pymol import cmd
from typing import Optional, List, Dict, Union

# Import the client for communicating with the Docker server
from .client import PyMOLDockingClient

# Create the client instance with configurable URL
SERVER_URL = "http://localhost:5000"
docker_client = PyMOLDockingClient(SERVER_URL)

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def assert_organic(selection: str) -> None:
    """
    Assert that the given PyMOL selection consists of only organic molecules.

    Parameters:
        selection: The PyMOL selection string to check.

    Raises:
        AssertionError: If the selection contains non-organic molecules.
    """ 
    # Create a temporary selection for organic molecules
    cmd.select("organic_check", f"{selection} and organic")

    # Compare atom counts to verify selection is organic
    total_atoms = cmd.count_atoms(selection)
    organic_atoms = cmd.count_atoms("organic_check")

    # Check if all atoms are organic
    assert organic_atoms == total_atoms, f"Selection '{selection}' contains non-organic molecules."
    logger.info(f"Selection '{selection}' contains only organic molecules.")

    # Delete the temporary selection
    cmd.delete("organic_check")

@cmd.extend
def off_site_docking(protein_selection: str, ligand_smiles: str, outname: str) -> None:
    """
    Perform off-site docking using PyMOL protein selection and a SMILES string for the ligand.
    
    This function saves the protein selection to a temporary file, sends the data to the Docker
    server for docking, and loads the results back into PyMOL.
    
    Parameters:
        protein_selection: The PyMOL selection string for the protein.
        ligand_smiles: The SMILES string for the ligand.
        outname: The base name for the output file.
    
    Returns:
        None. The docked structure is loaded into PyMOL.
    """
    # Check server connection
    if not docker_client.check_health():
        print("Error: Docker server is not running. Please start the Docker server first.")
        return
    
    try:
        # Save the protein to a temporary file
        protein_name = cmd.get_object_list(protein_selection)[0]
        to_save_protein = Path(gettempdir()) / f"{protein_name}.pdb"
        cmd.save(to_save_protein.as_posix(), protein_selection)

        # Define crystal file path
        crystal_file = Path("./Crystal.sdf")
        if not crystal_file.exists():
            print(f"Warning: Crystal file not found at {crystal_file}. Docking may fail.")

        # Run the docking process
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
        # cmd.load(str(results["prepared_protein"]), f"{outname}_protein")
        print(f"Docking completed successfully! Results loaded as {outname} and {outname}_protein")
        
    except FileNotFoundError as e:
        print(f"Error: Required file not found - {e}")
    except AssertionError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Error during docking: {e}")
        logger.exception("Detailed docking error:")

@cmd.extend
def on_site_docking(protein_selection: str, ligand_selection: str, mode: str, 
                   outname: str, minimization_flag: bool = False) -> None:
    """
    Perform on-site docking using PyMOL selections for protein and ligand.

    This function saves the selected protein and ligand as temporary files, runs 
    the docking or minimization process using the Docker server, and loads the 
    results back into PyMOL.

    Parameters:
        protein_selection: The PyMOL selection string for the protein.
        ligand_selection: The PyMOL selection string for the ligand.
        mode: The operation mode, either "Minimize" or "Dock".
        outname: The base name for the output files.
        minimization_flag: Whether to perform minimization after docking.

    Raises:
        AssertionError: If the mode is invalid or the ligand is not organic.

    Returns:
        None. The result of the docking/minimization is loaded into PyMOL.
    """
    # Check server connection
    if not docker_client.check_health():
        print("Error: Docker server is not running. Please start the Docker server first.")
        return
    
    # Validate inputs
    if mode not in ["Minimize", "Dock"]:
        print(f"Error: Invalid mode '{mode}'. Must be 'Minimize' or 'Dock'.")
        return

    try:
        # Get object names from selections
        protein_name = cmd.get_object_list(protein_selection)[0]
        ligand_name = cmd.get_object_list(ligand_selection)[0]
        assert_organic(ligand_name)

        # Create temp files
        to_save_protein = Path(gettempdir()) / f"{protein_name}.pdb"
        to_save_ligand = Path(gettempdir()) / f"{ligand_name}.sdf"

        cmd.save(str(to_save_protein), protein_selection)
        cmd.save(str(to_save_ligand), ligand_selection)

        # Define crystal file path
        crystal_file = Path("./Crystal.sdf")
        if not crystal_file.exists():
            print(f"Warning: Crystal file not found at {crystal_file}. Docking may fail.")

        # Run docking
        operation_type = "docking with minimization" if minimization_flag else mode.lower()
        print(f"Running {operation_type} with {protein_name} and {ligand_name}")
        
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
        # cmd.load(str(results["prepared_protein"]), f"{outname}_protein")
        
        # If minimization was performed, load the minimized complex
        if minimization_flag and "minimized_complex" in results:
            cmd.load(str(results["minimized_complex"]), f"{outname}_complex")
            print(f"Operation completed successfully! Results loaded as {outname}, {outname}_protein, and {outname}_complex")
        else:
            print(f"Operation completed successfully! Results loaded as {outname} and {outname}_protein")
            
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except AssertionError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Error during {mode.lower()}: {e}")
        logger.exception(f"Detailed {mode.lower()} error:")

@cmd.extend
def md_minimization(protein_selection: str, ligand_selection: Optional[str] = None, outname: str = "minimized_complex") -> None:
    """
    Perform molecular dynamics minimization on a protein-ligand complex.
    
    This function saves the selected protein and ligand as temporary files, runs the
    minimization process using the Docker server, and loads the minimized complex 
    structure back into PyMOL.
    
    Parameters:
        protein_selection: The PyMOL selection string for the protein.
        ligand_selection: The PyMOL selection string for the ligand (optional).
        outname: The name for the output complex file.
    
    Returns:
        None. The minimized complex is loaded into PyMOL.
    """
    # Check server connection
    if not docker_client.check_health():
        print("Error: Docker server is not running. Please start the Docker server first.")
        return
    
    try:
        # Get object names from selections
        protein_name = cmd.get_object_list(protein_selection)[0]
        
        # Create temp files
        to_save_protein = Path(gettempdir()) / f"{protein_name}.pdb"
        cmd.save(str(to_save_protein), protein_selection)
        
        # Process ligand if provided
        to_save_ligand = None
        if ligand_selection:
            ligand_name = cmd.get_object_list(ligand_selection)[0]
            assert_organic(ligand_name)
            to_save_ligand = Path(gettempdir()) / f"{ligand_name}.sdf"
            cmd.save(str(to_save_ligand), ligand_selection)
            print(f"Running MD minimization with {protein_name} and {ligand_name}")
        else:
            print(f"Running MD minimization with {protein_name} only")
        
        # Run minimization using the inplace_minimization function
        results = docker_client.inplace_minimization(
            protein_file=to_save_protein,
            ligand=to_save_ligand,
            output_name=outname,
            output_dir="."
        )
        
        # Load the results into PyMOL
        cmd.load(str(results["minimized_complex"]), outname)
        print(f"MD minimization completed successfully! Result loaded as {outname}")
        
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except AssertionError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"Error during MD minimization: {e}")
        logger.exception("Detailed MD minimization error:")

# Define the dialog class that inherits from QDialog
class PymolDockingDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super(PymolDockingDialog, self).__init__(parent)
        
        # Load the UI file directly
        uifile = os.path.join(os.path.dirname(__file__), 'GUI.ui')
        loadUi(uifile, self)
        
        # Create UI elements that might not be in the UI file
        # Add protein-only checkbox if not in UI
        if not hasattr(self, 'protein_only_flag'):
            self.protein_only_flag = QtWidgets.QCheckBox("Protein-only minimization", self.tab_2)
            self.protein_only_flag.setGeometry(150, 30, 200, 20)
            
        # Add output name field if not in UI
        if not hasattr(self, 'output_name_MD'):
            self.output_name_MD = QtWidgets.QLineEdit(self.tab_2)
            self.output_name_MD.setGeometry(150, 180, 301, 25)
            self.output_name_MD.setPlaceholderText("Output name (default: minimized_complex)")
            
        # Add labels if needed
        self.output_name_label = QtWidgets.QLabel("Output Name:", self.tab_2)
        self.output_name_label.setGeometry(70, 180, 80, 25)
        
        # Ensure state choosers have at least one option
        self.protein_state_chooser.addItem("1")
        self.ligand_state_chooser.addItem("1")
        
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
        # Connect protein-only checkbox toggle event
        self.protein_only_flag.toggled.connect(self.on_protein_only_toggled)
        # Connect tab widget change signal
        self.tabWidget.currentChanged.connect(self.on_tab_changed)
        # Connect buttonBox_MD for the second tab
        self.buttonBox_MD.accepted.connect(self.on_md_dialog_accepted)
        
        # Connect state update signals for MD tab
        self.protein_chooser_MD.currentIndexChanged.connect(self.update_protein_states)
        self.ligand_chooser_MD.currentIndexChanged.connect(self.update_ligand_states)
        
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
        
        # Set the ligand_chooser visibility based on protein_only_flag
        self.ligand_chooser_MD.setVisible(not self.protein_only_flag.isChecked())
        self.ligand_state_label.setVisible(not self.protein_only_flag.isChecked())
        self.ligand_state_chooser.setVisible(not self.protein_only_flag.isChecked())
        
        # Initialize state selectors
        self.update_protein_states()
        self.update_ligand_states()
    
    def update_protein_states(self):
        """Update the protein state selector with available states."""
        # Clear existing state options
        self.protein_state_chooser.clear()
        
        # Get current protein selection
        protein_selection = self.protein_chooser_MD.currentText()
        if not protein_selection:
            return
        
        try:
            # Get the object name from the selection
            # If selection includes chain specifier (e.g., "protein & chain A"), extract the object name
            obj_name = protein_selection.split()[0] if "&" in protein_selection else protein_selection
            
            # Count states for the selected object
            state_count = cmd.count_states(obj_name)
            
            # Populate state selector with states (1-based indexing in PyMOL)
            self.protein_state_chooser.addItems([str(i) for i in range(1, state_count + 1)])
            
            # Select first state by default
            if state_count > 0:
                self.protein_state_chooser.setCurrentIndex(0)
                
        except Exception as e:
            logger.error(f"Error updating protein states: {e}")
            # Default to a single state if error occurs
            self.protein_state_chooser.addItem("1")
    
    def update_ligand_states(self):
        """Update the ligand state selector with available states."""
        # Clear existing state options
        self.ligand_state_chooser.clear()
        
        # Only update if ligand selection is visible
        if not self.ligand_chooser_MD.isVisible():
            return
            
        # Get current ligand selection
        ligand_selection = self.ligand_chooser_MD.currentText()
        if not ligand_selection:
            return
        
        try:
            # Get the object name from the selection
            obj_name = ligand_selection.split()[0] if "&" in ligand_selection else ligand_selection
            
            # Count states for the selected object
            state_count = cmd.count_states(obj_name)
            
            # Populate state selector with states
            self.ligand_state_chooser.addItems([str(i) for i in range(1, state_count + 1)])
            
            # Select first state by default
            if state_count > 0:
                self.ligand_state_chooser.setCurrentIndex(0)
                
        except Exception as e:
            logger.error(f"Error updating ligand states: {e}")
            # Default to a single state if error occurs
            self.ligand_state_chooser.addItem("1")

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
            
            # Make sure state selectors are properly initialized
            if self.protein_chooser_MD.count() > 0:
                self.update_protein_states()
            if self.ligand_chooser_MD.count() > 0 and not self.protein_only_flag.isChecked():
                self.update_ligand_states()
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
        protein_state = int(self.protein_state_chooser.currentText())
        
        # Validate output name (allow only alphanumeric characters, underscore, and hyphen)
        outname = self.output_name_MD.text() if self.output_name_MD.text() else "minimized_complex"
        if not re.match(r'^[a-zA-Z0-9_\-]+$', outname):
            print(f"Error: Output name '{outname}' contains invalid characters. Using 'minimized_complex' instead.")
            logger.warning(f"Invalid output name: {outname}. Using default.")
            outname = "minimized_complex"
        
        # Only use ligand if protein_only_flag is not checked
        ligand = None
        ligand_state = None
        if not self.protein_only_flag.isChecked():
            ligand = self.ligand_chooser_MD.currentText()
            ligand_state = int(self.ligand_state_chooser.currentText())
        
        # Provide clear feedback on what's being run
        if ligand:
            print(f"Running MD minimization on protein-ligand complex...")
            print(f"  - Protein: {protein} (State: {protein_state})")
            print(f"  - Ligand: {ligand} (State: {ligand_state})")
            print(f"  - Output name: {outname}")
            logger.info(f"Running protein-ligand minimization with output name: {outname}")
            
            # Create temporary selections with specific states
            protein_with_state = f"{protein} and state {protein_state}"
            ligand_with_state = f"{ligand} and state {ligand_state}"
            
            # Call md_minimization with state-specific selections
            md_minimization(protein_with_state, ligand_with_state, outname)
        else:
            print(f"Running MD minimization on protein only...")
            print(f"  - Protein: {protein} (State: {protein_state})")
            print(f"  - Output name: {outname}")
            logger.info(f"Running protein-only minimization with output name: {outname}")
            
            # Create temporary selection with specific state
            protein_with_state = f"{protein} and state {protein_state}"
            
            # Call md_minimization with state-specific selection
            md_minimization(protein_with_state, None, outname)
    
    def on_protein_only_toggled(self, checked):
        """Handle the protein-only checkbox state change."""
        # Show/hide ligand selector based on checkbox state
        self.ligand_chooser_MD.setVisible(not checked)
        self.ligand_state_label.setVisible(not checked)
        self.ligand_state_chooser.setVisible(not checked)
        if checked:
            logger.info(f"Protein-only minimization mode activated")
            print("Protein-only minimization mode activated")
        else:
            logger.info(f"Protein-ligand minimization mode activated")
            print("Protein-ligand minimization mode activated")

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