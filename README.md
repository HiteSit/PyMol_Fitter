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
    docker pull hitesit/pymol-fitter:latest
    ```
    
    For GPU-accelerated version:
    
    ```bash
    docker pull hitesit/pymol-fitter:gpu
    ```

### 1b. Set Up the Docker Container

1. Clone this repository:
    
    ```bash
    git clone https://github.com/HiteSit/PyMol_Fitter.git
    cd PyMol_Fitter
    ```
    
2. Build and start the Docker container:
    
    ```bash
    cd docker
    docker-compose --profile cpu up -d  # Uses the CPU version (hitesit/pymol-fitter:latest)
    ```
    
    If you have an NVIDIA GPU and want to use it for acceleration:
    
    ```bash
    cd docker
    docker-compose --profile gpu up -d  # Uses the GPU version (hitesit/pymol-fitter:gpu)
    ```
    
    You can also run both services simultaneously:
    
    ```bash
    cd docker
    docker-compose --profile cpu --profile gpu up -d  # Runs both CPU and GPU versions
    ```

### 2. Install the PyMOL Plugin

1. Install the required Python packages for the client:
   
    > A helpful tip to simplify the installation of these lightweight dependencies is to download PyMOL3 and install the dependencies using the PyMOL2 prompt. Where can you find it? If you have PyMOL3, you can simply search for it in the Start Menu.
    
    ```bash
    pip install requests pathlib
    ```
    
2. Install the plugin from the Pymol-Plugin GUI using the URL:
    
    ```bash
    https://github.com/HiteSit/PyMol_Fitter/blob/master/pymol_fitter_plugin.zip
    ```
    
    ![](./docs/img1.png)
    

## 🚀 Usage Guide

### Starting the Plugin

1. Open PyMOL
2. The plugin will be available in the Plugin menu as "Pymol Fitter"
3. Click on it to open the plugin interface

### Environment Preparation

To perform docking, the plugin requires information about the binding site. This information is retrieved from a crystal ligand. For the plugin to work correctly in your current directory (`$PWD`), you must have a file named `Crystal.sdf` that will be used to define the binding site.

Here's how to create this file:
```python
cd /path/to/working_directory

# Click on the Crystal Ligand
save Crystal.sdf, sele
```

### In-Site Docking (Using 3D Structures)

Use this mode when you already have both protein and ligand structures loaded in PyMOL.

1. In the plugin dialog, select "In-Site" from the modality dropdown
2. Choose your protein from the "Protein" dropdown
3. Choose your ligand from the "Ligand" dropdown
4. Select the operation mode:
    - **Dock**: For full docking (finding binding poses)
    - **Minimize**: For optimizing an existing binding pose.
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

This mode allows you to perform molecular dynamics-based minimization (vacuum) on a protein or protein-ligand complex.

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
- The (Docker) console will display information about the process

## ⚙️ Technical Details

### Docking Workflow

1. **Protein Preparation**:
    - Removal of non-standard residues and heteroatoms
    - Addition of missing atoms and hydrogens
    - Optimization of hydrogen bonds and minimization
    - Default on Protoss, fallback on PDBFixer
2. **Ligand Preparation**:
    - For SMILES input: 3D structure generation
    - For SDF input: Standardization and hydrogen addition
    - CDPKit-based Protonation
3. **Docking Process**:
    - Binding site detection using reference ligand
    - Docking with Smina
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
    # Check for CPU version
    docker ps | grep pymol-fitter
    
    # Check for GPU version
    docker ps | grep pymol-fitter-gpu
    ```
    
    If no results, try restarting the container:
    
    ```bash
    # For CPU version (hitesit/pymol-fitter:latest)
    cd docker
    docker-compose down
    docker-compose --profile cpu up -d
    ```
    
    For GPU version:
    ```bash
    # For GPU version (hitesit/pymol-fitter:gpu)
    cd docker
    docker-compose down
    docker-compose --profile gpu up -d
    ```

2. **API Connection Errors**:
Test the API directly:
    
    ```bash
    # For CPU version (port 5000)
    curl http://localhost:5000/health
    
    # For GPU version (port 5001)
    curl http://localhost:5001/health
    ```
    
    If not responding, check Docker logs:
    
    ```bash
    # For CPU version
    docker logs pymol-fitter
    ```
    
    For GPU version:
    ```bash
    # For GPU version
    docker logs pymol-fitter-gpu
    ```
    
3. **Plugin Not Appearing in PyMOL**:
    - Ensure the plugin files are in the correct directory
    - Restart PyMOL
    - Check PyMOL's plugin manager

4. **In-Site Docking Failures**:
    - Always check the Docker console; the error handling should be descriptive enough.
    - Check for the kekulization of the ligand and consider starting from a PDBx file instead of a PDB file. If it still fails, use the Off-Site option with the kekulized SMILES of the Crystal ligand of interest.
  
5. **Minimization Failures**:
    - Check the Docker console. You should see something like the example below. This process can take several minutes depending on the complexity of the ligand.
        ```shell
        Generating a residue template for [H][C]1=[N][N]([H])[C]([H])=[C]1[c]1[c]([H])[c]([H])[c]([H])[c]([C@]([H])([N]([H])[c]2[n][c]([C]([H])([H])[H])[n][c]3[c]([H])[c]([O][C]([H])([H])[H])[c]([O][C]([H])([H])[H])[c]([H])[c]23)[C]([H])([H])[H])[c]1[H] using gaff-2.11
        ```

## 📚 Scientific Methods

### Docking Algorithm

This plugin uses Smina.

### Energy Minimization

The minimization is performed using OpenMM with:

- AMBER ff14SB force field for proteins
- GAFF2 for small molecules
  
The minimization process is highly customizable and encapsulated into a single function (`pymol_fitter_server/pymol_fitter_src/Protein_Minimization.py`). Since it leverages `SystemGenerator`, users can easily swap protein and ligand forcefields as well as incorporate explicit water models.

### Pose Evaluation

PoseBusters is used to evaluate docking poses.

## 📋 Citation

If you use this plugin in your research, please cite:

```
@software{pymol_fitter_plugin,
  author = {Riccardo Fusco},
  title = {PyMOL Fitter Plugin},
  year = {2025},
  url = {https://github.com/hitesit/pymol-fitter}
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
- [Smina](https://sourceforge.net/projects/smina/)
- [OpenMM](http://openmm.org/)
- [RDKit](https://www.rdkit.org/)
- [PoseBusters](https://github.com/OpenFreeEnergy/poseBusters)
- [CDPK](https://github.com/molinfo-vienna/CDPKit)
- [Flask](https://flask.palletsprojects.com/)

## 📝 Documentation

This documentation was created with assistance from [Claude 3.7](https://www.anthropic.com/claude).

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.