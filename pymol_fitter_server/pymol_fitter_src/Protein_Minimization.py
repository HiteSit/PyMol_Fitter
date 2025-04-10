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
import openmm
from openmm import (
    Platform, LangevinIntegrator, MonteCarloBarostat,
    CustomExternalForce, State, System, Context
)

# OpenMM Units
from openmm import unit
from openmm.unit import Quantity

# OPENFF
from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator, SMIRNOFFTemplateGenerator

# Type for OpenMM unit quantities
UnitQuantity = Quantity


def minimize_complex(
    prot_path: Union[str, Path], 
    lig_mol: Optional[Chem.rdchem.Mol] = None,
    use_implicit_solvent: bool = True
) -> Dict[str, Union[str, float]]:
    """
    Prepare and minimize a protein-ligand complex using OpenMM.
    
    This function:
    1. Processes the protein structure to fix missing atoms and add hydrogens
    2. Combines the protein with the provided ligand molecule
    3. Sets up a simulation system with appropriate forcefields
    4. Performs energy minimization on the complex
    5. Returns the minimized structure as a PDB string
    
    Parameters:
        prot_path: Path to the protein PDB file
        lig_mol: RDKit molecule object representing the ligand
        use_implicit_solvent: If True, use implicit solvent (GBSA) for minimization
        
    Returns:
        Dictionary containing:
        - "PDB_BEFORE": PDB string of the initial complex
        - "PDB_AFTER": PDB string of the minimized complex
        - "delta_energy": Energy difference before/after minimization (kJ/mol)
        
    Raises:
        FileNotFoundError: If the protein file doesn't exist
        ValueError: If the ligand molecule is invalid
        RuntimeError: If the minimization process fails
    """
    # Check if protein file exists
    if not os.path.exists(str(prot_path)):
        raise FileNotFoundError(f"Protein file not found: {prot_path}")
    
    if lig_mol is not None and lig_mol.GetNumAtoms() == 0:
        raise ValueError("Invalid ligand molecule: molecule is None or empty")
    
    try:
        # Fix the protein
        fixer: PDBFixer = PDBFixer(str(prot_path))
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)
        
        # Create the basic forcefield kwargs
        forcefield_kwargs: Dict[str, Any] = { 
            'constraints': app.HBonds, 
            'hydrogenMass': 4*unit.amu 
        }
        
        # Determine whether to use implicit solvent and set up appropriate parameters
        if use_implicit_solvent:
            # Need to include the implicit solvent force field XML along with the protein force field
            forcefields = ['amber/ff14SB.xml', 'amber/tip3p_standard.xml', 'implicit/obc2.xml']
            
            # Direct forcefield approach instead of passing through SystemGenerator for implicit solvent
            forcefield = app.ForceField(*forcefields)
            
            # Create modeller based on whether we have a ligand
            if lig_mol is not None:
                # Parse the ligand
                ligand_mol: Molecule = Molecule.from_rdkit(lig_mol)
                lig_top = ligand_mol.to_topology()
                
                # Merge the ligand into the protein
                modeller: Modeller = Modeller(fixer.topology, fixer.positions)
                modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
                
                # Need to parameterize the ligand separately
                template_generator = SMIRNOFFTemplateGenerator(
                    molecules=[ligand_mol], forcefield='openff-2.1.1'
                )
                forcefield.registerTemplateGenerator(template_generator.generator)
            else:
                modeller: Modeller = Modeller(fixer.topology, fixer.positions)
            
            # Ensure no periodic boundary conditions are set
            modeller.topology.setPeriodicBoxVectors(None)
            
            # Create the system directly with OBC implicit solvent
            system = forcefield.createSystem(
                modeller.topology,
                nonbondedMethod=app.CutoffNonPeriodic, 
                nonbondedCutoff=1.0*unit.nanometer,
                constraints=app.HBonds,
                hydrogenMass=4*unit.amu
            )
            # Add the implicit solvent force directly
            # The implicitSolvent parameter is now specified in the force field xml file
            
        else:
            # Non-implicit solvent approach using SystemGenerator
            nonperiodic_forcefield_kwargs: Dict[str, Any] = {
                'nonbondedMethod': app.NoCutoff,
            }
            
            forcefields = ['amber/ff14SB.xml']
                
            if lig_mol is not None:
                # Parse the ligand
                ligand_mol: Molecule = Molecule.from_rdkit(lig_mol)
                lig_top = ligand_mol.to_topology()
                
                # Merge the ligand into the protein
                modeller: Modeller = Modeller(fixer.topology, fixer.positions)
                modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
                
                # Set up the system generator with the small molecule forcefield
                system_generator: SystemGenerator = SystemGenerator(
                    forcefields=forcefields,
                    small_molecule_forcefield='openff-2.1.1',
                    molecules=[ligand_mol],
                    forcefield_kwargs=forcefield_kwargs,
                    nonperiodic_forcefield_kwargs=nonperiodic_forcefield_kwargs
                )
            else:
                modeller: Modeller = Modeller(fixer.topology, fixer.positions)
                
                # Set up the system generator WITHOUT small molecule forcefield
                system_generator: SystemGenerator = SystemGenerator(
                    forcefields=forcefields,
                    forcefield_kwargs=forcefield_kwargs,
                    nonperiodic_forcefield_kwargs=nonperiodic_forcefield_kwargs
                )
            
            # Ensure no periodic boundary conditions are set
            modeller.topology.setPeriodicBoxVectors(None)
            
            # Create the system for simulation
            system: System = system_generator.create_system(modeller.topology)
        
        # Check for and remove any unwanted periodic box vectors
        # This ensures we're truly in non-periodic mode
        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                if force.getNonbondedMethod() in [
                    openmm.NonbondedForce.PME,
                    openmm.NonbondedForce.Ewald,
                    openmm.NonbondedForce.CutoffPeriodic
                ]:
                    force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)
        
        # Set up the integrator
        integrator: LangevinIntegrator = LangevinIntegrator(
            300 * unit.kelvin,
            1 / unit.picosecond,
            0.002 * unit.picoseconds,
        )
        
        # Try to use CUDA, fall back to CPU if not available
        try:
            platform: Platform = Platform.getPlatformByName('CUDA')
            proprieties: Dict[str, str] = {'Precision': 'mixed', 'CudaDeviceIndex': "0"}
        except Exception:
            platform: Platform = Platform.getPlatformByName('CPU')
            proprieties: Dict[str, str] = {}

        # Set up the simulation environment
        simulation: Simulation = Simulation(
            modeller.topology, 
            system,
            integrator, 
            platform=platform, 
            platformProperties=proprieties
        )
        
        # Set up the context and perform minimization
        context: Context = simulation.context
        context.setPositions(modeller.positions)
        
        # Minimize and get energies
        min_state: State = simulation.context.getState(getEnergy=True, getPositions=True)
        energy_before_min: UnitQuantity = min_state.getPotentialEnergy()
        
        # Perform energy minimization
        simulation.minimizeEnergy()
        min_state = simulation.context.getState(getEnergy=True, getPositions=True)
        energy_after_min: UnitQuantity = min_state.getPotentialEnergy()
        
        # Generate PDB strings for before and after states
        with io.StringIO() as PDB_before:
            PDBFile.writeFile(modeller.topology, modeller.positions, PDB_before)
            PDB_BEFORE: str = PDB_before.getvalue()
        
        with io.StringIO() as PDB_after:
            PDBFile.writeFile(simulation.topology, min_state.getPositions(), PDB_after)
            PDB_AFTER: str = PDB_after.getvalue()
        
        # Calculate energy values and difference
        energy_before_min_val: float = energy_before_min.value_in_unit(unit.kilojoule_per_mole)
        energy_after_min_val: float = energy_after_min.value_in_unit(unit.kilojoule_per_mole)
        delta_energy = energy_after_min_val - energy_before_min_val
        
        return {
            "PDB_BEFORE": PDB_BEFORE,
            "PDB_AFTER": PDB_AFTER,
            "delta_energy": round(delta_energy, 2)
        }
        
    except Exception as e:
        raise RuntimeError(f"Error during complex minimization: {str(e)}") from e

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