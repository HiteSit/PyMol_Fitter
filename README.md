# PyMOL Docking Docker

A PyMOL plugin for molecular docking that uses Docker to enable cross-platform compatibility, especially for Windows users who may encounter dependency issues with computational chemistry software.

## Features

- Protein preparation with ProtoSS
- Molecular docking with Smina
- Complex minimization with OpenMM
- Support for both SDF files and SMILES strings as ligand input
- Cross-platform compatibility via Docker
- Simple PyMOL GUI integration

## Project Structure

This project is organized into two main components:

1. **PyMOL Plugin (Client-Side)**
   - Located in `pymol_docking_plugin/`
   - Contains the PyMOL plugin code and GUI
   - Communicates with the Docker container via API calls

2. **Docker Server (Computational Backend)**
   - Located in `pymol_docking_server/` and `docker/`
   - Contains the Flask API and computational code
   - Runs in a Docker container to ensure consistent environment

## How It Works

1. The PyMOL plugin (client) sends protein and ligand data to the Docker server
2. The Docker server performs the docking calculations
3. Results are sent back to the PyMOL plugin and loaded into the PyMOL session

This architecture allows the plugin to work on any platform (including Windows) that can run Docker, regardless of dependency issues.

## Requirements

- Docker
- Docker Compose (optional, but recommended)
- PyMOL
- Python 3.6+ with `requests` and `pathlib` packages

## Installation

### 1. Set Up the Docker Container

1. Build and start the Docker container:
   ```bash
   cd docker
   docker-compose up -d
   ```

2. Verify that the container is running:
   ```bash
   docker ps | grep pymol-docking-server
   ```

### 2. Install the PyMOL Plugin

1. Install the required dependencies for the plugin:
   ```bash
   pip install requests pathlib
   ```

2. Copy the `pymol_docking_plugin` directory to PyMOL's plugin directory:
   ```bash
   # For Linux/macOS
   cp -r pymol_docking_plugin ~/.pymol/startup/
   
   # For Windows (adjust path as needed)
   # Copy to C:\Users\<username>\AppData\Local\PyMOL\startup
   ```

## Usage

### From the PyMOL GUI

1. Open PyMOL
2. The plugin will be available in the Plugin menu as "Pymol Docking"
3. Select your protein and ligand
4. Choose the docking mode and other options:
   - **In-Site**: Use a 3D ligand structure already loaded in PyMOL
   - **Off-Site**: Use a SMILES string to generate a ligand
5. Click "OK" to start the docking process

### Using PyMOL Commands

The plugin adds several commands to PyMOL:

1. **Off-site docking** (using SMILES string):
   ```python
   off_site_docking protein_selection, "SMILES_STRING", output_name
   ```

2. **On-site docking** (using 3D structure):
   ```python
   on_site_docking protein_selection, ligand_selection, mode, output_name, minimization_flag
   ```
   
   Where:
   - `mode` is either "Dock" or "Minimize"
   - `minimization_flag` is `True` or `False`

3. **Complex minimization**:
   ```python
   docker_minimize_complex protein_selection, ligand_selection, output_name
   ```

### Using the Python Client Directly

For advanced usage, you can use the client API directly:

```python
from pymol_docking_plugin.client import PyMOLDockingClient

# Initialize the client
client = PyMOLDockingClient("http://localhost:5000")

# Docking a ligand to a protein
results = client.dock(
    protein_file="path/to/protein.pdb",
    ligand="path/to/ligand.sdf",  # OR SMILES string if is_smiles=True
    is_smiles=False,
    dock_mode="Dock",  # or "Minimize"
    output_name="my_docking_result",
    output_dir="./results"
)

# The results dictionary contains paths to the output files
print(results)
```

## Troubleshooting

### Docker Server Not Running

If you see a message "Docker server is not running", check:

1. The Docker container is running:
   ```bash
   docker ps | grep pymol-docking-server
   ```

2. The container is accessible on port 5000:
   ```bash
   curl http://localhost:5000/health
   ```

3. Restart the container if needed:
   ```bash
   cd docker
   docker-compose restart
   ```

### PyMOL Plugin Errors

If you encounter errors in the PyMOL plugin:

1. Check that the required Python dependencies are installed
2. Ensure the Docker server is running
3. Look at the PyMOL log for detailed error messages

## Development

### Docker Server

- `pymol_docking_server/app.py`: Flask API for handling docking requests
- `pymol_docking_server/pymol_docking_src/`: Computational code for docking
- `docker/Dockerfile`: Docker configuration
- `docker/environment.yml`: Conda environment for the Docker container

### PyMOL Plugin

- `pymol_docking_plugin/__init__.py`: PyMOL plugin initialization and commands
- `pymol_docking_plugin/client.py`: Client for communicating with the Docker server
- `pymol_docking_plugin/GUI.ui`: Qt UI file for the plugin dialog

## Examples

The `examples` directory contains sample files for testing the plugin:
- `LAC3.pdb`: Sample protein structure
- `Ligand.sdf`: Sample ligand structure
- `Crystal.sdf`: Sample crystal structure
- `Notebook.ipynb`: Jupyter notebook with usage examples

## License

[Your license information here]

## Acknowledgements

- Original PyMOL docking plugin developers
- This project uses Smina, OpenMM, ProtoSS, and other open-source tools 