# PyMOL Fitter Plugin

A powerful PyMOL plugin for molecular docking and minimization that leverages Docker containerization to provide robust cross-platform compatibility and consistent computational environments, especially beneficial for Windows users who often encounter dependency issues with computational chemistry software.

## 📋 Features

- **Cross-platform compatibility** through Docker containerization
- **Multiple docking modes**:
    - In-Site: Use 3D ligand structures already loaded in PyMOL
    - Off-Site: Generate 3D conformers from SMILES strings
- **Structure preparation and optimization**:
    - Automated protein preparation
    - Automated ligand preparation and protonation
    - OpenMM-powered minimization of protein-ligand complexes
- **Evaluation and analysis**:
    - Pose assessment with PoseBusters
    - Detailed docking logs and scores
- **User-friendly GUI** integrated into PyMOL's interface

## 🏗️ Architecture

The project is organized into three main components:

1. **PyMOL Plugin (Client)** - `pymol_fitter_plugin/`
    - GUI interface integrated with PyMOL
    - Client code to communicate with the Docker server
    - No computational code - all calculations happen in the Docker container
2. **Server Code** - `pymol_fitter_server/`
    - Flask API endpoints for docking and minimization
    - Computational code using smina, OpenMM, and other libraries
    - Functions for protein preparation, ligand preparation, and pose assessment
3. **Docker Configuration** - `docker/`
    - Dockerfile and docker-compose.yml
    - Environment configuration with all dependencies
    - Data volume for persistent storage of results

## 🔧 Requirements

### System Requirements

- **GPU**: NVIDIA GPU with CUDA support (optional, for faster calculations)
- **Disk Space**: At least 10GB free space for Docker images and data

### Software Requirements

- [Docker](https://www.docker.com/products/docker-desktop/)
- [Pymol3](https://pymol.org/)

## 📥 Installation

### 1a. Pull docker container

1. Simply pull the container
    
    ```bash
    docker pull hitesit/pymol-fitter-minimizer
    ```
    

### 1b. Set Up the Docker Container

1. Clone this repository:
    
    ```bash
    git clone https://github.com/yourusername/pymol-fitter.git
    cd pymol-fitter
    ```
    
2. Build and start the Docker container:
    
    ```bash
    cd docker
    docker-compose up -d
    ```
    

### 2. Install the PyMOL Plugin

1. Install the required Python packages for the client:
    
    ```bash
    pip install requests pathlib
    ```
    
2. Install the plugin from the Pymol-Plugin GUI using the URL:
    
    ```bash
     https://github.com/HiteSit/Pymol_Fitter/blob/master/pymol_fitter_plugin.zip
    ```
    
    ![](./docs/img1.png)
    

## 🚀 Usage Guide

### Starting the Plugin

1. Open PyMOL
2. The plugin will be available in the Plugin menu as "Pymol Fitter"
3. Click on it to open the plugin interface

### In-Site Docking (Using 3D Structures)

Use this mode when you already have both protein and ligand structures loaded in PyMOL.

1. In the plugin dialog, select "In-Site" from the modality dropdown
2. Choose your protein from the "Protein" dropdown
3. Choose your ligand from the "Ligand" dropdown
4. Select the operation mode:
    - **Dock**: For full docking (finding binding poses)
    - **Minimize**: For optimizing an existing binding pose
5. Enter an output name for the results
6. Click "OK" to start the docking process

### Off-Site Docking (Using SMILES)

Use this mode when you have a protein structure but want to generate a ligand from a SMILES string.

1. In the plugin dialog, select "Off-Site" from the modality dropdown
2. Choose your protein from the "Protein" dropdown
3. Enter the SMILES string for your ligand in the text box
4. Enter an output name for the results
5. Click "OK" to start the docking process

### MD Minimization

This mode allows you to perform molecular dynamics-based minimization on a protein or protein-ligand complex.

1. In the plugin dialog, switch to the "MD Minimization" tab
2. Choose your protein from the "Protein" dropdown
3. If minimizing a complex, choose your ligand from the "Ligand" dropdown
    - Check "Protein-only minimization" if you only want to minimize the protein
4. Enter an output name for the minimized structure
5. Click "OK" to start the minimization process

### Results

After successful docking or minimization:

- The docked ligand will be loaded into PyMOL as a new object
- The prepared protein will be loaded as a separate object
- For minimization, the minimized complex will be loaded
- The console will display information about the process

## ⚙️ Technical Details

### Docking Workflow

1. **Protein Preparation**:
    - Removal of non-standard residues and heteroatoms
    - Addition of missing atoms and hydrogens
    - Optimization of hydrogen bonds and minimization
2. **Ligand Preparation**:
    - For SMILES input: 3D structure generation
    - For SDF input: Standardization and hydrogen addition
    - Conformer generation and optimization
3. **Docking Process**:
    - Binding site detection using reference ligand
    - Docking with smina (AutoDock Vina derivative)
    - Scoring and pose ranking
4. **Minimization**:
    - Force field-based energy minimization using OpenMM
    - AMBER ff14SB for proteins
    - GAFF2 for small molecules

### Server API Endpoints

The server exposes the following RESTful API endpoints:

- **GET `/health`**: Check server status
- **POST `/dock`**: Perform docking operation
- **POST `/minimize`**: Perform minimization operation

## 🔍 Troubleshooting

### Common Issues

1. **Docker Container Not Running**:
    
    ```bash
    docker ps | grep pymol-fitter-minimizer
    ```
    
    If no results, try restarting the container:
    
    ```bash
    cd docker
    docker-compose down
    docker-compose up -d
    ```
    
2. **API Connection Errors**:
Test the API directly:
    
    ```bash
    curl http://localhost:5000/health
    ```
    
    If not responding, check Docker logs:
    
    ```bash
    docker logs pymol-fitter-minimizer
    ```
    
3. **Plugin Not Appearing in PyMOL**:
    - Ensure the plugin files are in the correct directory
    - Restart PyMOL
    - Check PyMOL's plugin manager
4. **Docking/Minimization Failures**:
    - Check console output for specific errors
    - Ensure input structures are valid
    - For SMILES inputs, verify the SMILES string is correct
    - For minimization, check for clashes in the input structure

## 🔬 For Developers

### Adding New Server-side Features

1. Add new computational code to `pymol_fitter_server/pymol_fitter_src/`
2. Expose functionality via new endpoints in `pymol_fitter_server/app.py`
3. Rebuild the Docker container:
    
    ```bash
    cd docker && docker-compose up -d --build
    ```
    

### Extending the Client

1. Update `pymol_fitter_plugin/client.py` to call new API endpoints
2. Add UI elements to `pymol_fitter_plugin/GUI.ui` if needed
3. Add new PyMOL commands in `pymol_fitter_plugin/__init__.py`
4. Test the new functionality with the Docker server running

### Development Environment

For development, you may want to mount your code directory into the Docker container:

```yaml
# In docker-compose.ymlvolumes:  - ./data:/app/data  - ../pymol_fitter_server:/app
```

This allows you to modify the server code without rebuilding the container.

## 📚 Scientific Methods

### Docking Algorithm

This plugin uses smina, a fork of AutoDock Vina with enhanced scoring function options and improved minimization algorithms. The docking process:

1. Identifies the binding site using a reference ligand
2. Generates and scores multiple conformers
3. Returns the top-ranked pose(s)

### Energy Minimization

The minimization is performed using OpenMM with:

- AMBER ff14SB force field for proteins
- GAFF2 for small molecules
- Hydrogen mass repartitioning for increased time steps
- Amber-compatible water model (TIP3P)

### Pose Evaluation

PoseBusters is used to evaluate docking poses by checking:

- Atom clashes
- Bond lengths and angles
- Hydrogen bonding networks
- Comparison to reference crystal structures (if available)

## 📋 Citation

If you use this plugin in your research, please cite:

```
@software{pymol_fitter_plugin,
  author = {Your Name},
  title = {PyMOL Fitter Plugin},
  year = {2023},
  url = {https://github.com/yourusername/pymol-fitter}
}
```

Additionally, please cite the underlying tools:

- Smina (AutoDock Vina Fork)
- OpenMM
- PoseBusters
- RDKit
- PyMOL

## 🙏 Acknowledgments

This project utilizes several open-source tools and libraries:

- [PyMOL](https://pymol.org/)
- [smina](https://sourceforge.net/projects/smina/)
- [OpenMM](http://openmm.org/)
- [RDKit](https://www.rdkit.org/)
- [PoseBusters](https://github.com/OpenFreeEnergy/poseBusters)
- [CDPK](https://github.com/molinfo-vienna/CDPKit)
- [Flask](https://flask.palletsprojects.com/)

## 📝 Documentation

This documentation was created with assistance from [Claude 3.7](https://www.anthropic.com/claude), an AI assistant by Anthropic.

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.