# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Docker Environment Commands

### Building and Running Docker Containers

The project uses Docker containers for computational chemistry tasks. Build and run containers using:

```bash
# CPU version (default)
cd docker
docker-compose --profile cpu up -d

# GPU version (requires NVIDIA GPU with CUDA)
cd docker
docker-compose --profile gpu up -d

# Both versions simultaneously
docker-compose --profile cpu --profile gpu up -d
```

### Health Checks and Container Management

```bash
# Check if containers are running
docker ps | grep pymol-fitter

# Check container health
curl http://localhost:5000/health  # CPU version
curl http://localhost:5001/health  # GPU version

# View container logs
docker logs pymol-fitter      # CPU version
docker logs pymol-fitter-gpu  # GPU version

# Stop containers
docker-compose down
```

### Pre-built Images

Alternative to building locally:
```bash
# Pull pre-built images
docker pull hitesit/pymol-fitter:latest  # CPU version
docker pull hitesit/pymol-fitter:gpu     # GPU version
```

## Testing

### Test Structure
- `tests/scientific_tests/` - Core computational chemistry functionality
- `tests/application_tests/` - Flask API and integration tests (planned)

### Running Tests

```bash
# Run all tests
pytest tests/

# Run only scientific tests
pytest tests/scientific_tests/

# Run specific test file
pytest tests/scientific_tests/test_docking_engine.py

# Run with coverage
pytest --cov=pymol_fitter_server --cov=pymol_fitter_plugin tests/

# Skip slow tests (marked with @pytest.mark.slow)
pytest -k "not slow" tests/

# Run only slow tests
pytest -k "slow" tests/
```

### Test Requirements

Before running tests, ensure all dependencies are installed:
```bash
pip install pytest pytest-cov
pip install biopython biotite rdkit datamol requests flask
conda install -c conda-forge pdbfixer
```

Some tests require the Docker container to be running. Start it before running Docker-dependent tests.

## Architecture Overview

### Three-Component Architecture

1. **PyMOL Plugin (Client)** - `pymol_fitter_plugin/`
   - GUI integrated with PyMOL using Qt
   - Communicates with Docker server via REST API
   - No computational code - delegates all calculations to server
   - Main files: `__init__.py` (GUI), `client.py` (API client)

2. **Server Code** - `pymol_fitter_server/`
   - Flask API (`app.py`) with three main endpoints:
     - `/health` - Health check
     - `/dock` - Molecular docking operations
     - `/minimize` - Molecular dynamics minimization
     - `/virtual_screen` - Virtual screening with multiple ligands
   - Core computation modules in `pymol_fitter_src/`:
     - `Docking_Engine.py` - Main docking workflow using Smina
     - `Protein_Preparation.py` - Protein preprocessing
     - `Protein_Minimization.py` - OpenMM-based minimization
     - `CDPK_Utils.py` - Chemical structure utilities
     - `VS_Class.py` - Virtual screening functionality

3. **Docker Configuration** - `docker/`
   - `Dockerfile.cpu` and `Dockerfile.gpu` - Environment setup
   - `docker-compose.yml` - Service orchestration
   - `environment.yml` - Conda environment with scientific dependencies

### Key Workflows

**Docking Workflow:**
1. Protein preparation (Protoss/PDBFixer fallback)
2. Ligand preparation (SMILES→3D or SDF standardization)
3. Binding site detection from crystal reference
4. Smina docking
5. Optional OpenMM minimization with AMBER ff14SB/GAFF2

**Data Flow:**
1. PyMOL plugin saves selections to temp files
2. Client encodes files as base64 and sends to server
3. Server processes in temporary directories
4. Results returned as base64-encoded files
5. Plugin loads results back into PyMOL

### Required Environment Setup

For proper docking functionality:
- A `Crystal.sdf` file must exist in the working directory to define the binding site
- This file should contain the crystal ligand structure that defines where docking will occur

## Key Dependencies and Force Fields

- **Smina** - Docking algorithm (AutoDock Vina fork)
- **OpenMM** - Molecular dynamics minimization
  - AMBER ff14SB for proteins
  - GAFF2 for small molecules
- **PoseBusters** - Pose validation
- **CDPKit** - Chemical structure processing and protonation
- **RDKit** - Chemical informatics
- **Protoss/PDBFixer** - Protein preparation

## Development Notes

### Import Path Issues
Server modules use relative imports (`from pymol_fitter_src.X import Y`). If experiencing import errors during development:

```python
# Fix import paths by changing from:
from pymol_fitter_src.Docking_Engine import ...
# To:
from pymol_fitter_server.pymol_fitter_src.Docking_Engine import ...

# Or add to PYTHONPATH:
export PYTHONPATH=/path/to/PyMol_Fitter:$PYTHONPATH
```

### GUI Development
The PyMOL plugin uses Qt with a `.ui` file (`GUI.ui`) loaded dynamically. The main dialog class `PymolDockingDialog` handles three tabs:
- Docking (In-Site/Off-Site modes)
- MD Minimization
- Virtual Screening

### Server Communication
All client-server communication uses base64-encoded files in JSON payloads. The client (`PyMOLDockingClient`) handles encoding/decoding automatically.

### Scientific Validation
The system includes pose validation using PoseBusters to assess docking quality and identify potential issues with generated poses.

## Contacts Visualization Plugin

The repository includes `dkoes_show_contacts.py` (by David Koes) for protein-ligand contact visualization:
- Press **F3** to show contacts for any selected molecule/ligand
- Displays polar contacts, hydrogen bonds, and clashes with color coding
- Useful for analyzing docking results