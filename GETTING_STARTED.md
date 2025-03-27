# Getting Started with PyMOL Docking Docker

This guide will help you quickly get up and running with PyMOL Docking Docker.

## Prerequisites

Before you begin, make sure you have:

- Docker Desktop installed and running ([Install Docker](https://docs.docker.com/get-docker/))
- PyMOL installed ([Download PyMOL](https://pymol.org/))
- Python 3.6+ with `requests` and `pathlib` packages

## Quick Setup

### 1. Start the Docker Container

1. Open a terminal/command prompt
2. Navigate to the `docker` directory:
   ```bash
   cd docker
   ```
3. Start the Docker container:
   ```bash
   docker-compose up -d
   ```
4. Verify the container is running:
   ```bash
   docker ps | grep pymol-docking-server
   ```

### 2. Install the PyMOL Plugin

1. Install required Python packages:
   ```bash
   pip install requests pathlib
   ```

2. Copy the plugin to PyMOL's plugin directory:
   - **Linux/macOS**: Copy `pymol_docking_plugin` to `~/.pymol/startup/`
   - **Windows**: Copy `pymol_docking_plugin` to `C:\Users\<username>\AppData\Local\PyMOL\startup\`

## Basic Usage

### Start PyMOL

Open PyMOL. The plugin should be automatically loaded and available in the Plugins menu.

### Option 1: Using the GUI

1. Load your protein structure in PyMOL
2. Click on `Plugin → Pymol Docking`
3. Select the docking mode:
   - **In-Site**: For docking with a 3D ligand already loaded in PyMOL
   - **Off-Site**: For docking with a SMILES string

4. For In-Site docking:
   - Select your protein and ligand from the dropdown menus
   - Choose "Dock" or "Minimize" mode
   - Enter an output name
   - Click "OK"

5. For Off-Site docking:
   - Select your protein from the dropdown menu
   - Enter a SMILES string for the ligand
   - Enter an output name
   - Click "OK"

### Option 2: Using PyMOL Commands

You can use the plugin's commands directly in PyMOL:

1. Off-site docking (with SMILES):
   ```python
   off_site_docking protein_selection, "SMILES_STRING", output_name
   ```

2. On-site docking (with 3D ligand):
   ```python
   on_site_docking protein_selection, ligand_selection, "Dock", output_name, True
   ```

3. Minimization only:
   ```python
   docker_minimize_complex protein_selection, ligand_selection, output_name
   ```

## Example Walkthrough

1. Load a protein:
   ```python
   load pymol_docking_server/examples/LAC3.pdb
   ```

2. Option A: Off-site docking with SMILES
   ```python
   # Aspirin SMILES string
   off_site_docking LAC3, "CC(=O)OC1=CC=CC=C1C(=O)O", "aspirin_docked"
   ```

3. Option B: On-site docking with ligand file
   ```python
   load pymol_docking_server/examples/Ligand.sdf, my_ligand
   on_site_docking LAC3, my_ligand, "Dock", "ligand_docked", True
   ```

After docking, the results will be automatically loaded into PyMOL.

## Troubleshooting

If you encounter issues:

1. Check that the Docker container is running:
   ```bash
   docker ps | grep pymol-docking-server
   ```

2. Verify the API is accessible:
   ```bash
   curl http://localhost:5000/health
   ```

3. If needed, restart the Docker container:
   ```bash
   cd docker
   docker-compose restart
   ```

For more detailed information, consult the [README.md](README.md) file. 