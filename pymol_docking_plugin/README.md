# PyMOL Docking Plugin (Client)

This directory contains the client-side code for the PyMOL Docking Docker plugin.

## Contents

- `__init__.py`: The main plugin file with PyMOL commands and GUI integration
- `client.py`: Client code to communicate with the Docker server API
- `GUI.ui`: Qt UI file defining the plugin's graphical interface

## How It Works

This plugin communicates with a Flask API running in a Docker container to perform docking calculations. The workflow is:

1. User selects a protein and ligand in PyMOL
2. Plugin sends the data to the Docker server
3. Server performs calculations and returns results
4. Plugin loads the results back into PyMOL

## Installation

1. Install the required Python packages:
   ```bash
   pip install requests pathlib
   ```

2. Copy this directory to your PyMOL plugins directory:
   ```bash
   # For Linux/macOS
   cp -r pymol_docking_plugin ~/.pymol/startup/
   
   # For Windows (adjust path as needed)
   # Copy to C:\Users\<username>\AppData\Local\PyMOL\startup
   ```

## Usage

The plugin adds the following commands to PyMOL:

- `off_site_docking`: Dock a ligand (specified by SMILES) to a protein
- `on_site_docking`: Dock a 3D ligand to a protein
- `docker_minimize_complex`: Minimize a protein-ligand complex

You can also access the plugin via the PyMOL menu: Plugin → Pymol Docking

## Troubleshooting

If you encounter issues:

1. Make sure the Docker container is running
2. Check your network connection to localhost:5000
3. Look for error messages in the PyMOL console 