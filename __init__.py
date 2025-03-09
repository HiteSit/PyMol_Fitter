import os

def __init_plugin__(app=None):
    """
    Initialize the plugin by adding a menu item to PyMOL's Plugin menu.
    """
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('My Plugin', run_plugin_gui)

# Global reference to the dialog to prevent garbage collection
dialog = None

def run_plugin_gui():
    """
    Launch the GUI when the menu item is clicked.
    If the dialog doesn't exist, create it and load the UI file.
    """
    global dialog
    if dialog is None:
        from pymol.Qt import QtWidgets
        from pymol.Qt.utils import loadUi
        dialog = QtWidgets.QDialog()
        
        # Construct the path to the UI file (assumes 'my_ui.ui' is in the same directory as this script)
        uifile = os.path.join(os.path.dirname(__file__), './widget/GUI.ui')
        form = loadUi(uifile, dialog)
        
        # Connect signals (example: connect a button click to a function)
        def say_hello():
            print("Hello from plugin!")
        
        # Assuming there's a button named 'button_hello' in the UI file
        form.protein_chooser_1.currentTextChanged.connect(say_hello)
    
    # Show the dialog
    dialog.show()