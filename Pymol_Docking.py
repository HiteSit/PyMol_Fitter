# %%
import logging
import os
import shutil
import subprocess
import io
import time
import warnings
from pathlib import Path
from tempfile import gettempdir
from urllib.parse import urljoin
from typing import Optional, List, Tuple, Dict, Any, Union

from pymol import cmd
from openbabel import pybel
import datamol as dm
import rdkit.Chem as rdChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import requests
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
import CDPL.Chem as Chem
import CDPL.ConfGen as ConfGen
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

#################################################################################
########################### Protein Preparation Utilities #######################
#################################################################################

# ProteinsPlus API endpoints
PROTEINS_PLUS_URL = 'https://proteins.plus/api/v2/'
UPLOAD = urljoin(PROTEINS_PLUS_URL, 'molecule_handler/upload/')
UPLOAD_JOBS = urljoin(PROTEINS_PLUS_URL, 'molecule_handler/upload/jobs/')
PROTEINS = urljoin(PROTEINS_PLUS_URL, 'molecule_handler/proteins/')
LIGANDS = urljoin(PROTEINS_PLUS_URL, 'molecule_handler/ligands/')
PROTOSS = urljoin(PROTEINS_PLUS_URL, 'protoss/')
PROTOSS_JOBS = urljoin(PROTEINS_PLUS_URL, 'protoss/jobs/')

def poll_job(job_id, poll_url, poll_interval=1, max_polls=10):
    """
    Poll the progress of a job by continuously polling the server in regular intervals and updating the job information.

    Args:
        job_id (str): UUID of the job to poll.
        poll_url (str): URL to send the polling request to.
        poll_interval (int): Time interval between polls in seconds. Default is 1 second.
        max_polls (int): Maximum number of times to poll before exiting. Default is 10.

    Returns:
        dict: Polled job information.
    """
    job = requests.get(poll_url + job_id + '/').json()
    status = job['status']
    current_poll = 0

    while status == 'pending' or status == 'running':
        print(f'Job {job_id} is {status}')
        current_poll += 1

        if current_poll >= max_polls:
            print(f'Job {job_id} has not completed after {max_polls} polling requests and {poll_interval * max_polls} seconds')
            return job

        time.sleep(poll_interval)
        job = requests.get(poll_url + job_id + '/').json()
        status = job['status']

    print(f'Job {job_id} completed with {status}')
    return job

def prepare_protein_protoss(receptor: Path, receptor_protoss: Path) -> Path:
    """
    Prepares a protein using ProtoSS.

    Args:
        receptor (Path): Path to the input protein file in PDB format.
        receptor_protoss (Path): Path where the prepared protein file should be saved.

    Returns:
        Path: Path to the prepared protein file in PDB format.
    """
    print('Preparing protein with ProtoSS ...')

    with open(receptor) as upload_file:
        query = {'protein_file': upload_file}
        job_submission = requests.post(PROTOSS, files=query).json()

    protoss_job = poll_job(job_submission['job_id'], PROTOSS_JOBS)
    protossed_protein = requests.get(PROTEINS + protoss_job['output_protein'] + '/').json()
    protein_file = io.StringIO(protossed_protein['file_string'])
    protein_structure = PDBParser().get_structure(protossed_protein['name'], protein_file)

    with open(receptor_protoss, 'w') as output_file_handle:
        pdbio = PDBIO()
        pdbio.set_structure(protein_structure)
        pdbio.save(output_file_handle)
    
    return receptor_protoss

#################################################################################
########################### CDPK Utilities ######################################
#################################################################################

def molsToFiles(mols: list[Chem.BasicMolecule], output_file: Path) -> None:
    """
    Writes a list of molecules to the specified output file.

    Args:
        mols (list): List of CDPKit molecules to write to the output file.
        output_file (str): Path to the output file.
    """
    writer = Chem.MolecularGraphWriter(str(output_file))
    for mol in mols:
        writer.write(mol)
    writer.close()

def filesToMols(sdf_input: Path) -> List[Chem.BasicMolecule]:
    """
    Retrieves the structure data of each molecule in the provided SD file.

    Args:
        sdf_input (Path): Path to the input SD file.

    Returns:
        List[Chem.BasicMolecule]: List of parsed molecules.
    """
    reader = Chem.FileSDFMoleculeReader(str(sdf_input))
    mol = Chem.BasicMolecule()
    mols_list = []
    
    try:
        while reader.read(mol):
            try:
                if not Chem.hasStructureData(mol):
                    print('Error: no structure data available for molecule', Chem.getName(mol))
                    continue
                struct_data = Chem.getStructureData(mol)
                new_mol = Chem.BasicMolecule(mol)
                mols_list.append(new_mol)
            except Exception:
                pass
    except Exception:
        pass
    
    return mols_list

def standardize(chembl_proc: Chem.ChEMBLStandardizer,
                in_mol: Chem.Molecule, out_mol: Chem.Molecule,
                proc_excluded: bool, extract_parent: bool) -> Chem.ChEMBLStandardizer.ChangeFlags:
    """
    Performs ChEMBL molecule standardization and parent structure extraction.

    Args:
        chembl_proc (Chem.ChEMBLStandardizer): Instance of the Chem.ChEMBLStandardizer class.
        in_mol (Chem.Molecule): Input molecule to standardize.
        out_mol (Chem.Molecule): Output molecule to store the standardized molecule.
        proc_excluded (bool): If True, molecules flagged as excluded will be processed.
        extract_parent (bool): If True, the parent structure will be extracted.

    Returns:
        Chem.ChEMBLStandardizer.ChangeFlags: Flags indicating the carried out modifications.
    """
    change_flags = chembl_proc.standardize(in_mol, out_mol, proc_excluded)

    if extract_parent:
        change_flags &= ~Chem.ChEMBLStandardizer.EXCLUDED
        change_flags |= chembl_proc.getParent(out_mol)
    return change_flags

def getListOfChangesString(change_flags: Chem.ChEMBLStandardizer.ChangeFlags, verbose: bool = False) -> str:
    """Returns a string listing the carried out modifications."""
    if not verbose:
        return None

    changes = '   Carried out modifications:'
    change_list = [
        (Chem.ChEMBLStandardizer.EXPLICIT_HYDROGENS_REMOVED, 'Explicit hydrogens removed'),
        (Chem.ChEMBLStandardizer.UNKNOWN_STEREO_STANDARDIZED, 'Undefined stereocenter information standardized'),
        (Chem.ChEMBLStandardizer.BONDS_KEKULIZED, 'Kekule structure generated'),
        (Chem.ChEMBLStandardizer.STRUCTURE_NORMALIZED, 'Functional groups normalized'),
        (Chem.ChEMBLStandardizer.CHARGES_REMOVED, 'Number of charged atoms reduced'),
        (Chem.ChEMBLStandardizer.TARTRATE_STEREO_CLEARED, 'Configuration of chiral tartrate atoms set to undefined'),
        (Chem.ChEMBLStandardizer.STRUCTURE_2D_CORRECTED, '2D structure corrected'),
        (Chem.ChEMBLStandardizer.ISOTOPE_INFO_CLEARED, 'Isotope information cleared'),
        (Chem.ChEMBLStandardizer.SALT_COMPONENTS_REMOVED, 'Salt components removed'),
        (Chem.ChEMBLStandardizer.SOLVENT_COMPONENTS_REMOVED, 'Solvent components removed'),
        (Chem.ChEMBLStandardizer.DUPLICATE_COMPONENTS_REMOVED, 'Duplicate components removed')
    ]

    for flag, description in change_list:
        if change_flags & flag:
            changes += '\n    * ' + description

    return changes

def generate3dConformation(mol: Chem.Molecule, struct_gen: ConfGen.StructureGenerator) -> int:
    """
    Generates a low-energy 3D structure of the argument molecule.

    Args:
        mol (Chem.Molecule): Molecule to generate a 3D structure for.
        struct_gen (ConfGen.StructureGenerator): Instance of the ConfGen.StructureGenerator class.

    Returns:
        int: Status code indicating the success of the 3D structure generation.
    """
    ConfGen.prepareForConformerGeneration(mol)
    status = struct_gen.generate(mol)
    
    if status == ConfGen.ReturnCode.SUCCESS:
        struct_gen.setCoordinates(mol)
    
    return status

def stero_enumerator(input_sdf: Path, output_sdf: Path) -> Path:
    """Enumerates stereoisomers for molecules in an SDF file."""
    supplier = rdChem.SDMolSupplier(str(input_sdf))
    writer = rdChem.SDWriter(str(output_sdf))
    
    opts = StereoEnumerationOptions(unique=True)
    for mol in supplier:
        init_name = mol.GetProp("_Name")
        isomers: List[rdChem.Mol] = tuple(EnumerateStereoisomers(mol, options=opts))
        for i, isomer in enumerate(isomers):
            iso_name = init_name + "_Iso" + str(i)
            isomer.SetProp("_Name", iso_name)
            writer.write(isomer)
    writer.close()

    return output_sdf

class CDPK_Runner:
    """Class for running CDPK operations on molecules."""
    
    def __init__(self, standardize: bool = True, protonate: bool = True, gen3d: bool = True):
        self.chembl_proc = Chem.ChEMBLStandardizer()
        self.prot_state_gen = Chem.ProtonationStateStandardizer()
        
        self.standardize = standardize
        self.protonate = protonate
        self.gen3d = gen3d
        
    def prepare_ligands(self, input_sdf: Path, output_sdf: Path):
        """Prepares ligands with standardization, protonation, and 3D generation as configured."""
        mols: List[Chem.BasicMolecule] = filesToMols(input_sdf)
        writer = Chem.MolecularGraphWriter(str(output_sdf))
        
        for mol in mols:
            in_mol = mol
            out_mol = Chem.BasicMolecule()
            
            if self.standardize:
                change_flags = standardize(self.chembl_proc, in_mol, out_mol, True, True)
                change_flags_str = getListOfChangesString(change_flags, verbose=True)
            
            if self.protonate:
                self.prot_state_gen.standardize(out_mol, Chem.ProtonationStateStandardizer.PHYSIOLOGICAL_CONDITION_STATE)
                Chem.perceiveComponents(out_mol, True)

            if self.gen3d:
                Chem.setMultiConfExportParameter(writer, False)
                struct_gen = ConfGen.StructureGenerator()
                status = generate3dConformation(out_mol, struct_gen)
                if status == ConfGen.ReturnCode.SUCCESS:
                    Chem.setMDLDimensionality(out_mol, 3)
        
            writer.write(out_mol)
        
        writer.close()

#################################################################################
########################### Logging Setup #######################################
#################################################################################

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    handlers=[
                        logging.StreamHandler()
                    ])
logger = logging.getLogger(__name__)

#################################################################################
########################### Plants Docking ######################################
#################################################################################

class Plants_Docking:
    def __init__(self, protein_pdb: Path, input_ligands: Union[Path, str]):
        self.workdir: Path = Path(os.getcwd())
        self.protein_pdb: Path = protein_pdb
        self.crystal_sdf: Path = Path("Crystal.sdf")

        if isinstance(input_ligands, str):
            self.ligands_smiles: str = input_ligands
            self.input_mode = "SMILES"
        elif isinstance(input_ligands, Path):
            self.input_ligands: Path = input_ligands
            self.input_mode = "SDF"

        self.plants_env_var()

    def plants_env_var(self):
        os.environ['PATH'] = '/home/hitesit/Software/PLANTS:' + os.environ.get('PATH', '')

    def prepare_protein(self):
        protein_basename: str = self.protein_pdb.stem

        protein_tmp: Path = Path(gettempdir()) / f"{protein_basename}_tmp.pdb"
        protein_PREP: Path = self.workdir / f"{protein_basename}_PREP.mol2"

        fixer = PDBFixer(filename=str(self.protein_pdb))
        fixer.removeHeterogens(keepWater=False)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()

        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)

        PDBFile.writeFile(fixer.topology, fixer.positions, open(str(protein_tmp), 'w'))

        # Openbabel convert pdb to mol2
        prt = next(pybel.readfile("pdb", str(protein_tmp)))
        out = pybel.Outputfile("mol2", str(protein_PREP), overwrite=True)

        out.write(prt)
        out.close()

        self.protein_PREP: Path = protein_PREP
        return protein_PREP

    @staticmethod
    def cdpk_fixer(input_sdf: Path, TMP_basename: str) -> Path:
        TMP_fixed_sdf = Path(gettempdir()) / "TMP_fixed.sdf"
        
        supplier = rdChem.SDMolSupplier(input_sdf.as_posix())
        writer = rdChem.SDWriter(TMP_fixed_sdf.as_posix())
        for mol in supplier:
            mol.SetProp("_Name", TMP_basename)
            writer.write(mol)
        writer.close()
        
        cdpk_runner = CDPK_Runner(standardize=True, protonate=True, gen3d=True)
        cdpk_runner.prepare_ligands(TMP_fixed_sdf, TMP_fixed_sdf)
        
        return TMP_fixed_sdf.resolve()

    def _convert_to_mol2(self, input_sdf: Path) -> Path:
        """Convert prepared SDF file to MOL2 format required by PLANTS."""
        tmp_mol2 = Path(gettempdir()) / "ligand_prepared.mol2"
        
        # Read the first molecule from SDF and write as MOL2
        mol = next(pybel.readfile("sdf", str(input_sdf)))
        mol.write("mol2", str(tmp_mol2), overwrite=True)
        
        return tmp_mol2

    def prepare_ligands(self):
        if self.input_mode == "SDF":
            logger.info("Preparing ligands from SDF")
            # First prepare using CDPK (in SDF format)
            fixed_ligands_sdf = self.cdpk_fixer(self.input_ligands, "ligand")
            # Convert to MOL2 for PLANTS
            fixed_ligands_mol2 = self._convert_to_mol2(fixed_ligands_sdf)
            
            # Convert crystal to MOL2 as well
            crystal_mol2 = Path(gettempdir()) / "crystal.mol2"
            crystal = next(pybel.readfile("sdf", str(self.crystal_sdf)))
            crystal.write("mol2", str(crystal_mol2), overwrite=True)

            return fixed_ligands_mol2, crystal_mol2

        elif self.input_mode == "SMILES":
            logger.info("Preparing ligands from SMILES")
            # Convert SMILES to temporary SDF
            TMP_SMILES_SDF = Path(gettempdir()) / "TMP_SMILES.sdf"
            mol = dm.to_mol(self.ligands_smiles, sanitize=True, kekulize=True)
            mol.SetProp("_Name", "ligand")
            dm.to_sdf(mol, TMP_SMILES_SDF)
            
            # Prepare using CDPK
            fixed_ligands_sdf = self.cdpk_fixer(TMP_SMILES_SDF, "ligand")
            # Convert to MOL2 for PLANTS
            fixed_ligands_mol2 = self._convert_to_mol2(fixed_ligands_sdf)
            
            # Convert crystal to MOL2
            crystal_mol2 = Path(gettempdir()) / "crystal.mol2"
            crystal = next(pybel.readfile("sdf", str(self.crystal_sdf)))
            crystal.write("mol2", str(crystal_mol2), overwrite=True)
            
            return fixed_ligands_mol2, crystal_mol2

    def _define_binding_site(self):
        # Use the MOL2 crystal file directly
        _, crystal_mol2 = self.prepare_ligands()

        plants_command = f"plants.64bit --mode bind {str(crystal_mol2)} {str(self.protein_PREP)}"
        print(plants_command)

        box_res = subprocess.run(plants_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, cwd=self.workdir)
        box_infos = box_res.stdout.split("\n")

        binding_site_center = box_infos[-4]
        binding_site_radius = box_infos[-3]

        return binding_site_center, binding_site_radius

    def write_conf(self, ligand_to_dock: str, n_confs: int):
        # Define the binding site
        center, radius = self._define_binding_site()

        config_path = self.workdir / "config.txt"
        config_str = f"""### PLANTS configuration file
# scoring function and search settings
scoring_function chemplp
search_speed speed1

# input
protein_file {self.protein_PREP}
ligand_file {ligand_to_dock}

# output
output_dir plants_raw
write_multi_mol2 1

# binding site definition
{center}
{radius}

# cluster algorithm
cluster_structures {n_confs}
cluster_rmsd 1.0
"""

        with open(config_path, "w") as config_file:
            config_file.write(config_str)
        return config_path

    def _retrieve_docked_ligands(self):
        plants_raw = self.workdir / "plants_raw"
        docked_ligand: Path = plants_raw / "docked_ligands.mol2"

        shutil.copy(docked_ligand, self.workdir / "docked_ligands.mol2")

    def run_plants_docking(self, n_confs: int):
        protein_PREP: Path = self.prepare_protein()
        fixed_ligands, fixed_crystal = self.prepare_ligands()
        config_path: Path = self.write_conf(str(fixed_ligands), n_confs)

        runner = f"plants.64bit --mode screen {str(config_path)}"
        subprocess.run(runner, shell=True, check=True, cwd=self.workdir)

        self._retrieve_docked_ligands()
        return

#################################################################################
########################### Pymol Docking #######################################
#################################################################################

class Pymol_Docking:
    def __init__(self, protein_pdb: str, input_ligands: str):
        self.workdir: Path = Path(os.getcwd())
        self.protein_pdb: Path = Path(protein_pdb)
        self.crystal_sdf: Path = Path("Crystal.sdf")
        
        self.protein_preared = None

        if os.path.exists(input_ligands):
            self.ligands_sdf: Path = Path(input_ligands)
            self.input_mode = "SDF"
        else:
            self.ligands_smiles: str = input_ligands
            self.input_mode = "SMILES"

    def prepare_protein(self) -> Path:
        protein_basename: str = self.protein_pdb.stem
        protein_PREP: Path = self.workdir / f"{protein_basename}_PREP.pdb"
        
        prepare_protein_protoss(self.protein_pdb, protein_PREP)
        
        self.protein_preared = protein_PREP
        return protein_PREP

    @staticmethod
    def cdpk_fixer(input_sdf: Path, TMP_basename: str) -> None:
        TMP_fixed_sdf = Path(gettempdir()) / "TMP_fixed.sdf"
        
        supplier = rdChem.SDMolSupplier(input_sdf.as_posix())
        writer = rdChem.SDWriter(TMP_fixed_sdf.as_posix())
        for mol in supplier:
            mol.SetProp("_Name", TMP_basename)
            writer.write(mol)
        writer.close()
        
        cdpk_runner = CDPK_Runner(standardize=True, protonate=True, gen3d=True)
        cdpk_runner.prepare_ligands(TMP_fixed_sdf, TMP_fixed_sdf)
        
        return TMP_fixed_sdf.resolve()
    
    def prepare_ligands(self):
        if self.input_mode == "SDF":
            logger.info("Preparing ligands from SDF")
            fixed_crystal = self.crystal_sdf.resolve()
            fixed_ligands = self.cdpk_fixer(self.ligands_sdf, "ligand")

            return fixed_ligands.resolve(), fixed_crystal.resolve()
        
        elif self.input_mode == "SMILES":
            TMP_SMILES_SDF = Path(gettempdir()) / "TMP_SMILES.sdf"
            mol = dm.to_mol(self.ligands_smiles, sanitize=True, kekulize=True)
            mol.SetProp("_Name", "ligand")
            dm.to_sdf(mol, TMP_SMILES_SDF)
            
            fixed_crystal = self.crystal_sdf.resolve()
            fixed_ligands = self.cdpk_fixer(TMP_SMILES_SDF, "ligand")
            
            return fixed_ligands.resolve(), fixed_crystal.resolve()
        
    @staticmethod
    def pose_buster_processer(mol_pred: Path, mol_crystal: Path, mol_prot: Path):
        from posebusters import PoseBusters
        buster = PoseBusters()
        df = buster.bust(mol_pred, mol_crystal, mol_prot)
        df.reset_index(drop=False, inplace=True)
        bool_columns = df.select_dtypes(include=[bool]).columns
        success_rate_series = df[bool_columns].mean(axis=1) * 100
        success_rate_list = success_rate_series.tolist()
        
        # Read all molecules first
        mol_list = list(rdChem.SDMolSupplier(mol_pred.as_posix()))
        
        # Create a temporary file path
        temp_path = mol_pred.parent / f"temp_{mol_pred.name}"
        
        # Write to temporary file
        with rdChem.SDWriter(temp_path.as_posix()) as writer:
            for pose, success_rate in zip(mol_list, success_rate_list):
                pose.SetProp("Compliance", str(success_rate))  # Convert to string
                writer.write(pose)
        
        # Replace original file with temporary file
        temp_path.replace(mol_pred)
        
        return mol_pred
            
    def run_smina_docking(self, mode: str, docking_basename: str) -> Path:
        # Fix the protein
        print("Running protein preparation")
        protein_PKA: Path = self.prepare_protein()

        # Fix the ligands
        fixed_ligands, fixed_crystal = self.prepare_ligands()

        # Define the output
        smina_output: Path = self.workdir / f"{docking_basename}.sdf"
        smina_log = smina_output.with_suffix(".log")

        if mode == "Dock":
            smina_lst = [
                "smina",
                "-r", protein_PKA.as_posix(),
                "-l", fixed_ligands.as_posix(),
                "--autobox_ligand", fixed_crystal.as_posix(),
                "-o", str(smina_output),
                "--exhaustiveness", "32"
            ]

            print("Running docking")
            with open(smina_log, "w") as log_file:
                subprocess.run(smina_lst, check=True, stdout=log_file, stderr=log_file)

        elif mode == "Minimize":
            assert self.input_mode == "SDF", "Minimization only works with SDF input"
            smina_lst = [
                "smina",
                "-r", protein_PKA.as_posix(),
                "-l", fixed_ligands.as_posix(),
                "--autobox_ligand", fixed_crystal.as_posix(),
                "-o", smina_output,
                "--minimize", "--minimize_iters", "10"
            ]
            print("Running docking")
            with open(smina_log, "w") as log_file:
                subprocess.run(smina_lst, check=True, stdout=log_file, stderr=log_file)
        
        logger.info("Running pose buster")
        smina_bustered = self.pose_buster_processer(smina_output, fixed_crystal, protein_PKA)
        
        return smina_bustered

def assert_organic(selection):
    """
    Assert that the given PyMOL selection consists of only organic molecules.

    Parameters:
    - selection (str): The PyMOL selection string to check.

    Raises:
    - AssertionError: If the selection contains non-organic molecules.
    """
    # Create a temporary selection for organic molecules
    cmd.select("organic_check", f"{selection} and organic")

    # Get the total number of atoms in the original selection
    total_atoms = cmd.count_atoms(selection)

    # Get the number of atoms in the organic selection
    organic_atoms = cmd.count_atoms("organic_check")

    # Check if the number of organic atoms is equal to the total number of atoms
    assert organic_atoms == total_atoms, "Selection contains non-organic molecules."
    logger.info(f"Selection {selection} contains only organic molecules.")

    # Delete the temporary selection
    cmd.delete("organic_check")

# %%
@cmd.extend
def off_site_docking(protein_selection, ligand_smiles, outname:str):

    protein_name = cmd.get_object_list(protein_selection)[0]
    to_save_protein: Path = Path(gettempdir()) / f"{protein_name}.pdb"

    cmd.save(to_save_protein.as_posix(), protein_selection)

    pymol_docking = Pymol_Docking(to_save_protein.as_posix(), str(ligand_smiles))
    docked_sdf: Path = pymol_docking.run_smina_docking("Dock", outname)

    cmd.load(str(docked_sdf))

# %%
@cmd.extend
def on_site_docking(protein_selection, ligand_selection, mode, outname: str):
    """
    Perform on-site docking using PyMOL selections for protein and ligand, and smina for docking or minimization.

    This function saves the selected protein and ligand as temporary files, runs the docking or minimization process
    using smina, and then loads the docked ligand structure back into PyMOL.

    Parameters:
    - protein_selection (str): The PyMOL selection string for the protein.
    - ligand_selection (str): The PyMOL selection string for the ligand.
    - mode (str): The operation mode, either "Minimize" for energy minimization or "Dock" for docking.
    - outname (str): The base name for the output file. The docked structure will be saved as "{outname}.sdf"
                     and loaded into PyMOL with the same name.

    Raises:
    - AssertionError: If the mode is not one of the expected values ("Minimize" or "Dock").

    Returns:
    None. The result of the docking or minimization is loaded into PyMOL.
    """

    # Assert the mode
    assert mode in ["Minimize", "Dock"]

    protein_name = cmd.get_object_list(protein_selection)[0]
    ligand_name = cmd.get_object_list(ligand_selection)[0]
    assert_organic(ligand_name)

    to_save_protein = Path(gettempdir()) / f"{protein_name}.pdb"
    to_save_ligand = Path(gettempdir()) / f"{ligand_name}.sdf"

    cmd.save(str(to_save_protein), protein_selection)
    cmd.save(str(to_save_ligand), ligand_selection)

    pymol_docking = Pymol_Docking(str(to_save_protein), str(to_save_ligand))
    docked_sdf: Path = pymol_docking.run_smina_docking(mode, outname)

    cmd.load(str(docked_sdf))
# %%

#################################################################################
########################### Start of pymol plugin code ##########################
#################################################################################

def _get_select_list():
    '''
    Get either a list of object names, or a list of chain selections
    '''
    loaded_objects = [name for name in cmd.get_names('all', 1) if '_cluster_' not in name]

    # if single object, try chain selections
    if len(loaded_objects) == 1:
        chains = cmd.get_chains(loaded_objects[0])
        if len(chains) > 1:
            loaded_objects = ['{} & chain {}'.format(loaded_objects[0], chain) for chain in chains]

    return loaded_objects


class Pymol_Docking_GUI(object):
    ''' Qt version of the Plugin GUI '''

    def __init__(self):
        from pymol.Qt import QtWidgets
        dialog = QtWidgets.QDialog()
        self.setupUi(dialog)
        self.populate_ligand_select_list()
        self.choose_docking_modes()

        # Choose modality
        self.choose_setting()

        # Connect modality chooser signal to a method that clears and repopulates the selectors
        self.mode_chooser_2.currentIndexChanged.connect(self.clear_and_repopulate_selectors)

        # Connect the accepted signal to a method that checks the modality
        dialog.accepted.connect(self.on_dialog_accepted)

        dialog.exec_()

    def clear_and_repopulate_selectors(self):
        """Clear and repopulate the contents of the protein and ligand selectors."""
        self.protein_chooser_1.clear()
        self.protein_chooser_2.clear()
        self.ligand_chooser_1.clear()
        self.mode_chooser.clear()

        # Repopulate the selectors with available options
        self.populate_ligand_select_list()
        self.choose_docking_modes()

    def on_dialog_accepted(self):
        modality = self.mode_chooser_2.currentText().strip()

        if modality == "In-Site":
            self.in_site_wrapper()
            logger.info("In-Site docking selected")
        elif modality == "Off-Site":
            self.off_site_wrapper()
            logger.info("Off-Site docking selected")
        else:
            logger.warning("No valid modality selected")

    def in_site_wrapper(self):
        s1 = self.protein_chooser_1.currentText()
        s2 = self.ligand_chooser_1.currentText()
        mode = self.mode_chooser.currentText()
        outname = self.output_chooser.toPlainText()
        logger.info("Running on-site docking...")
        on_site_docking(s1, s2, mode, outname)

    def off_site_wrapper(self):
        s1 = self.protein_chooser_2.currentText()
        s2 = self.smile_chooser_2.toPlainText()
        outname = self.output_chooser_2.toPlainText()
        logger.info("Running off-site docking...")
        off_site_docking(s1, str(s2), outname)

    def choose_setting(self):
        possible_choices = ["In-Site", "Off-Site"]
        self.mode_chooser_2.addItems(possible_choices)

    def choose_docking_modes(self):
        possible_choices = ["Minimize", "Dock"]
        self.mode_chooser.addItems(possible_choices)

    def populate_ligand_select_list(self):
        loaded_objects = _get_select_list()

        self.protein_chooser_1.clear()
        self.protein_chooser_2.clear()
        self.ligand_chooser_1.clear()

        self.protein_chooser_1.addItems(loaded_objects)
        self.protein_chooser_2.addItems(loaded_objects)
        self.ligand_chooser_1.addItems(loaded_objects)

    def setupUi(self, Form):
        from PyQt5 import QtCore, QtWidgets
        Form.setObjectName("Form")
        Form.resize(662, 302)
        self.buttonBox = QtWidgets.QDialogButtonBox(Form)
        self.buttonBox.setGeometry(QtCore.QRect(440, 270, 193, 28))
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayoutWidget = QtWidgets.QWidget(Form)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(170, 20, 271, 91))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.protein_chooser_1 = QtWidgets.QComboBox(self.verticalLayoutWidget)
        self.protein_chooser_1.setObjectName("protein_chooser_1")
        self.verticalLayout.addWidget(self.protein_chooser_1)
        self.ligand_chooser_1 = QtWidgets.QComboBox(self.verticalLayoutWidget)
        self.ligand_chooser_1.setObjectName("ligand_chooser_1")
        self.verticalLayout.addWidget(self.ligand_chooser_1)
        self.mode_chooser = QtWidgets.QComboBox(self.verticalLayoutWidget)
        self.mode_chooser.setObjectName("mode_chooser")
        self.verticalLayout.addWidget(self.mode_chooser)
        self.verticalLayoutWidget_2 = QtWidgets.QWidget(Form)
        self.verticalLayoutWidget_2.setGeometry(QtCore.QRect(170, 140, 271, 118))
        self.verticalLayoutWidget_2.setObjectName("verticalLayoutWidget_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_2)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.protein_chooser_2 = QtWidgets.QComboBox(self.verticalLayoutWidget_2)
        self.protein_chooser_2.setObjectName("protein_chooser_2")
        self.verticalLayout_2.addWidget(self.protein_chooser_2)
        self.smile_chooser_2 = QtWidgets.QPlainTextEdit(self.verticalLayoutWidget_2)
        self.smile_chooser_2.setObjectName("smile_chooser_2")
        self.verticalLayout_2.addWidget(self.smile_chooser_2)
        self.output_chooser = QtWidgets.QTextEdit(Form)
        self.output_chooser.setGeometry(QtCore.QRect(460, 20, 151, 91))
        self.output_chooser.setObjectName("output_chooser")
        self.mode_chooser_2 = QtWidgets.QComboBox(Form)
        self.mode_chooser_2.setGeometry(QtCore.QRect(40, 140, 73, 22))
        self.mode_chooser_2.setObjectName("mode_chooser_2")
        self.output_chooser_2 = QtWidgets.QTextEdit(Form)
        self.output_chooser_2.setGeometry(QtCore.QRect(460, 140, 151, 111))
        self.output_chooser_2.setObjectName("output_chooser_2")

        self.buttonBox.accepted.connect(Form.accept)
        self.buttonBox.rejected.connect(Form.reject)


def __init__(self):
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('Pymol Docking', Pymol_Docking_GUI)
        return
    except Exception as e:
        print(e)