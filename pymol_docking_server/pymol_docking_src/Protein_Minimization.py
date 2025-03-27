import numpy as np
import io
from tempfile import NamedTemporaryFile
from pdbfixer import PDBFixer
from typing import List, Dict, Any, Optional, Union, Tuple, TypeVar, IO, BinaryIO, TextIO
from pathlib import Path
import os
import tempfile

# Biotite
import biotite.structure.io.pdb as pdb
from biotite.structure import AtomArray

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem

# OpenMM Application Layer
from openmm import app
from openmm.app import (
    Modeller, Simulation, PDBFile, DCDReporter, 
    StateDataReporter, CheckpointReporter
)

# OpenMM Library Layer
from openmm import (
    Platform, LangevinIntegrator, MonteCarloBarostat,
    CustomExternalForce, State, System, Context
)

# OpenMM Units
from openmm import unit
from openmm.unit import Quantity

# OPENFF
from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator

# Type for OpenMM unit quantities
UnitQuantity = Quantity


def minimize_complex(prot_path: Union[str, Path], lig_mol: Chem.rdchem.Mol) -> Dict[str, Union[str, float]]:
    """
    Prepare and minimize a protein-ligand complex using OpenMM.
    
    This function:
    1. Processes the protein structure to fix missing atoms and add hydrogens
    2. Combines the protein with the provided ligand molecule
    3. Sets up a simulation system with appropriate forcefields
    4. Performs energy minimization on the complex
    5. Returns the minimized structure as a PDB string
    
    Parameters
    ----------
    prot_path : Union[str, Path]
        Path to the protein PDB file
    lig_mol : rdkit.Chem.rdchem.Mol
        RDKit molecule object representing the ligand
        
    Returns
    -------
    Dict[str, Union[str, float]]
        Dictionary containing:
        - "PDB_BEFORE": PDB string of the initial complex
        - "PDB_AFTER": PDB string of the minimized complex
        - "energy_before_min": Energy before minimization (kJ/mol)
        - "energy_after_min": Energy after minimization (kJ/mol)
    """
    # Fix the protein
    fixer: PDBFixer = PDBFixer(str(prot_path))
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    
    # Parse the ligand
    ligand_mol: Molecule = Molecule.from_rdkit(lig_mol)
    lig_top = ligand_mol.to_topology()
    
    # Merge the ligand into the protein
    modeller: Modeller = Modeller(fixer.topology, fixer.positions)
    modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
    
    # Create the forcefield
    forcefield_kwargs: Dict[str, Any] = { 
        'constraints': app.HBonds, 
        # 'rigidWater': True, 
        # 'removeCMMotion': False, 
        'hydrogenMass': 4*unit.amu 
    }
    
    system_generator: SystemGenerator = SystemGenerator(
        forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
        small_molecule_forcefield='gaff-2.11',
        molecules=[ligand_mol],
        forcefield_kwargs=forcefield_kwargs
    )
    
    system: System = system_generator.create_system(modeller.topology)
    integrator: LangevinIntegrator = LangevinIntegrator(
        300 * unit.kelvin,
        1 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    
    platform: Platform = Platform.getPlatformByName('CUDA')
    proprieties: Dict[str, str] = {'Precision': 'mixed', 'CudaDeviceIndex': "0"}

    simulation: Simulation = Simulation(
        modeller.topology, 
        system,
        integrator, 
        platform=platform, 
        platformProperties=proprieties
    )
    
    # Set context
    context: Context = simulation.context
    context.setPositions(modeller.positions)
    
    # Minimize
    min_state: State = simulation.context.getState(getEnergy=True, getPositions=True)
    energy_before_min: UnitQuantity = min_state.getPotentialEnergy()
    simulation.minimizeEnergy()
    min_state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy_after_min: UnitQuantity = min_state.getPotentialEnergy()
    
    # with io.StringIO() as pdb_string:
    #     PDBFile.writeHeader(simulation.topology, pdb_string)
    #     PDBFile.writeModel(modeller.topology, modeller.positions, pdb_string, modelIndex=1)
    #     PDBFile.writeModel(simulation.topology, min_state.getPositions(), pdb_string, modelIndex=2)
    #     PDBFile.writeFooter(simulation.topology, pdb_string)
    #     PDB_MIN = pdb_string.getvalue()
        
    with io.StringIO() as PDB_before:
        PDBFile.writeFile(modeller.topology, modeller.positions, PDB_before)
        PDB_BEFORE: str = PDB_before.getvalue()
    
    with io.StringIO() as PDB_after:
        PDBFile.writeFile(simulation.topology, min_state.getPositions(), PDB_after)
        PDB_AFTER: str = PDB_after.getvalue()
    
    # Get the energies
    energy_before_min: float = energy_before_min.value_in_unit(unit.kilojoule_per_mole)
    energy_after_min: float = energy_after_min.value_in_unit(unit.kilojoule_per_mole)
    
    delta_energy = energy_after_min - energy_before_min
    
    return {
        "PDB_BEFORE": PDB_BEFORE,
        "PDB_AFTER": PDB_AFTER,
        "delta_energy": round(delta_energy, 2)
    }

def assign_bond_orders_from_smiles(
    pdb_input: Union[str, Path, io.StringIO], 
    smiles: str, 
    output_sdf: Optional[Union[str, Path]] = None
) -> Optional[Chem.rdchem.Mol]:
    """
    Assign bond orders to a 3D molecule using a SMILES string.
    
    Parameters:
    -----------
    pdb_input : Union[str, Path, io.StringIO]
        The PDB file as a path, PDB content as a string, or StringIO object
    smiles : str
        SMILES string of the ligand
    output_sdf : Optional[Union[str, Path]]
        Path to save the output SDF file. If None, will return the molecule object
    
    Returns:
    --------
    Optional[rdkit.Chem.rdchem.Mol]
        RDKit molecule with assigned bond orders if output_sdf is None, otherwise None
    """
    # Handle different input types for PDB
    temp_file: Optional[tempfile.NamedTemporaryFile] = None
    
    try:
        if isinstance(pdb_input, (str, Path)):
            # Check if it's a file path or PDB content string
            if os.path.exists(str(pdb_input)):
                # It's a file path
                pdb_file: str = str(pdb_input)
            else:
                # It's a PDB content string
                temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb', mode='w')
                temp_file.write(pdb_input)
                temp_file.close()
                pdb_file = temp_file.name
        
        elif isinstance(pdb_input, io.StringIO):
            # It's a StringIO object
            temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb', mode='w')
            temp_file.write(pdb_input.getvalue())
            temp_file.close()
            pdb_file = temp_file.name
            
        else:
            raise ValueError("pdb_input must be a string, Path, or StringIO object")
        
        # Read the PDB file
        mol_3d: Chem.rdchem.Mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol_3d is None:
            raise ValueError(f"Could not read PDB content")
        
        # Create molecule from SMILES
        mol_2d: Chem.rdchem.Mol = Chem.MolFromSmiles(smiles)
        if mol_2d is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        
        # Get 3D coordinates from the PDB molecule
        conf: Chem.rdchem.Conformer = mol_3d.GetConformer()
        coords: np.ndarray = np.array([conf.GetAtomPosition(i) for i in range(mol_3d.GetNumAtoms())])
        
        # Prepare the template molecule from SMILES with hydrogen atoms
        mol_template: Chem.rdchem.Mol = Chem.AddHs(mol_2d)
        
        # Match atoms between the molecules
        num_atoms_3d: int = mol_3d.GetNumAtoms()
        num_atoms_template: int = mol_template.GetNumAtoms()
        
        if num_atoms_3d != num_atoms_template:
            print(f"Warning: Different number of atoms. PDB: {num_atoms_3d}, SMILES: {num_atoms_template}")
        
        # Create a new molecule with correct bond orders
        mol_result: Chem.rdchem.Mol = AllChem.AssignBondOrdersFromTemplate(mol_template, mol_3d)
        
        # If successful, save or return
        if mol_result is None:
            raise ValueError("Failed to assign bond orders. Molecules may not be compatible.")
        
        if output_sdf:
            # Write to SDF file
            writer: Chem.SDWriter = Chem.SDWriter(str(output_sdf))
            writer.write(mol_result)
            writer.close()
            print(f"Successfully saved molecule with bond orders to {output_sdf}")
            return None
        else:
            return mol_result
            
    finally:
        # Clean up temporary file if created
        if temp_file and os.path.exists(temp_file.name):
            os.unlink(temp_file.name)

def extract_unk_residue(pdb_str: str) -> Path:
    """
    Extract the UNK residue from a PDB string.
    
    Parameters
    ----------
    pdb_str : str
        PDB string containing the structure
        
    Returns
    -------
    biotite.structure.AtomArray
        Atom array containing only the UNK residue
    """
    pbd_io: io.StringIO = io.StringIO(pdb_str)
    pdb_struct: AtomArray = pdb.PDBFile.read(pbd_io).get_structure(model=1)
    unk_struct: AtomArray = pdb_struct[pdb_struct.res_name == "UNK"]
    
    unk_file: pdb.PDBFile = pdb.PDBFile()
    unk_file.set_structure(unk_struct)
    
    with NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_file:
        unk_file.write(tmp_file.name)
        return Path(tmp_file.name)