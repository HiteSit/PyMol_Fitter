# PyMOL Docking Docker

A PyMOL plugin for molecular docking that enables cross-platform compatibility through Docker, especially beneficial for Windows users who may encounter dependency issues with computational chemistry software.

## Project Structure

The project is organized into three main components:

1. **PyMOL Plugin (Client)** - `pymol_docking_plugin/`
   - GUI interface for PyMOL
   - Client code to communicate with the Docker server
   - No computational code - all calculations happen in the Docker container

2. **Server Code** - `pymol_docking_server/`
   - Flask API endpoints 
   - Computational code for docking, minimization, etc.
   - Example files for testing

3. **Docker Configuration** - `docker/`
   - Dockerfile and docker-compose.yml 
   - Environment configuration
   - Data volume for persistent storage

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

1. Install the required Python packages for the client:
   ```bash
   pip install requests pathlib
   ```

2. Copy the PyMOL plugin directory to your PyMOL plugins folder:
   ```bash
   # For Linux/macOS
   cp -r pymol_docking_plugin ~/.pymol/startup/
   
   # For Windows (adjust path as needed)
   # Copy to C:\Users\<username>\AppData\Local\PyMOL\startup
   ```

## Usage

### From PyMOL

1. Open PyMOL
2. The plugin will be available in the Plugin menu as "Pymol Docking"
3. Select your protein and ligand
4. Choose docking options:
   - **In-Site**: Use a 3D ligand structure already loaded in PyMOL
   - **Off-Site**: Use a SMILES string to generate a ligand

## Troubleshooting

If you encounter issues:

1. Make sure the Docker container is running:
   ```
   docker ps | grep pymol-docking-server
   ```

2. Test the API directly:
   ```
   curl http://localhost:5000/health
   ```

3. Check PyMOL's console for error messages

## Development

### Adding New Features

1. **Server-side changes**:
   - Add code to `pymol_docking_server/pymol_docking_src/`
   - Expose functionality via the Flask API in `pymol_docking_server/app.py`
   - Rebuild the Docker container: `docker-compose up -d --build`

2. **Client-side changes**:
   - Update `pymol_docking_plugin/client.py` to call new API endpoints
   - Add new PyMOL commands in `pymol_docking_plugin/__init__.py`

## Examples

The `pymol_docking_server/examples/` directory contains sample files for testing. 