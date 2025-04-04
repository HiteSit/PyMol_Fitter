import io
import time
import warnings
from pathlib import Path
from urllib.parse import urljoin
import tempfile
import glob

import requests
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

class ProteinPreparation_Protoss:
    """
    A class for preparing protein structures using ProtoSS service.
    """
    
    def __init__(self):
        """Initialize the ProteinPreparation class with API endpoints."""
        self.PROTEINS_PLUS_URL = 'https://proteins.plus/api/v2/'
        self.UPLOAD = urljoin(self.PROTEINS_PLUS_URL, 'molecule_handler/upload/')
        self.UPLOAD_JOBS = urljoin(self.PROTEINS_PLUS_URL, 'molecule_handler/upload/jobs/')
        self.PROTEINS = urljoin(self.PROTEINS_PLUS_URL, 'molecule_handler/proteins/')
        self.LIGANDS = urljoin(self.PROTEINS_PLUS_URL, 'molecule_handler/ligands/')
        self.PROTOSS = urljoin(self.PROTEINS_PLUS_URL, 'protoss/')
        self.PROTOSS_JOBS = urljoin(self.PROTEINS_PLUS_URL, 'protoss/jobs/')
    
    def poll_job(self, job_id, poll_url, poll_interval=1, max_polls=10):
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
        # Get the initial job information
        job = requests.get(poll_url + job_id + '/').json()
        status = job['status']
        current_poll = 0

        # Continuously poll the job until it is completed or maximum polls reached
        while status == 'pending' or status == 'running':
            print(f'Job {job_id} is {status}')
            current_poll += 1

            # Check if maximum polls reached
            if current_poll >= max_polls:
                print(f'Job {job_id} has not completed after {max_polls} polling requests and {poll_interval * max_polls} seconds')
                return job

            # Wait for the specified interval before polling again
            time.sleep(poll_interval)

            # Poll the job again to get updated status
            job = requests.get(poll_url + job_id + '/').json()
            status = job['status']

        print(f'Job {job_id} completed with {status}')
        return job

    def prepare_protein_protoss(self, input_pdb: Path, output_pdb: Path) -> Path:
        """
        Prepares a protein using ProtoSS.

        Args:
            input_pdb (Path): Path to the input protein file in PDB format.
            output_pdb (Path): Path to save the prepared protein file.

        Returns:
            Path: Path to the prepared protein file in PDB format.
        """
        # Print log message
        print('Preparing protein with ProtoSS ...')

        # Convert CIF to PDB if needed
        temp_pdb = None
        if input_pdb.suffix.lower() == '.cif':
            print('Converting CIF to PDB format...')
            temp_pdb = Path(tempfile.mktemp(suffix='.pdb'))
            
            # Read the CIF file
            cif_file = pdbx.CIFFile.read(str(input_pdb))
            
            # Get the structure from the CIF file
            structure = pdbx.get_structure(cif_file, model=1)  # Get the first model
            
            # Create a PDB file object and set the structure
            pdb_file = pdb.PDBFile()
            pdb_file.set_structure(structure)
            
            # Write the PDB file
            pdb_file.write(str(temp_pdb))
            
            # Use the temporary PDB file for processing
            input_pdb_for_processing = temp_pdb
            print(f'Converted CIF to PDB: {temp_pdb}')
        else:
            input_pdb_for_processing = input_pdb

        # Open the receptor protein file
        with open(input_pdb_for_processing) as upload_file:
            # Create the query with the protein file
            query = {'protein_file': upload_file}
            # Submit the job to ProtoSS and get the job submission response
            job_submission = requests.post(self.PROTOSS, files=query).json()

        # Poll the job status until it is completed
        protoss_job = self.poll_job(job_submission.get('job_id'), self.PROTOSS_JOBS)

        # Get the output protein information from the job
        protossed_protein = requests.get(self.PROTEINS + protoss_job['output_protein'] + '/').json()

        # Create a StringIO object with the protein file string
        protein_file = io.StringIO(protossed_protein['file_string'])

        # Parse the protein structure from the StringIO object
        protein_structure = PDBParser().get_structure(protossed_protein['name'], protein_file)
        
        # Ensure the output directory exists
        output_pdb.parent.mkdir(parents=True, exist_ok=True)

        # Open the output file in write mode
        with output_pdb.open('w') as output_file_handle:
            # Create a PDBIO object
            pdbio = PDBIO()
            # Set the protein structure for saving
            pdbio.set_structure(protein_structure)
            # Save the protein structure to the output file
            pdbio.save(output_file_handle)
        
        # Clean up temporary file if created
        if temp_pdb and temp_pdb.exists():
            temp_pdb.unlink()
        
        # Return the path to the prepared protein file
        return output_pdb
    
    def __call__(self, input_pdb: Path, output_pdb: Path) -> Path:
        """
        Call method that wraps prepare_protein_protoss for easier usage.
        
        Args:
            input_pdb (Path): Path to the input protein file in PDB or CIF format.
            output_pdb (Path): Path to save the prepared protein file.
            
        Returns:
            Path: Path to the prepared protein file in PDB format.
        """
        return self.prepare_protein_protoss(input_pdb, output_pdb)
    
    
class ProteinPreparation_PDBFixer:
    """
    A class for preparing protein structures using PDBFixer.
    """
    
    def __call__(self, input_pdb: Path, output_pdb: Path) -> Path:
        """
        Call method that wraps prepare_protein_pdb_fixer for easier usage.
        """
        
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        fixer = PDBFixer(filename=str(input_pdb))
        fixer.removeHeterogens(True)

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()

        fixer.addMissingHydrogens(7.0)

        PDBFile.writeFile(fixer.topology, fixer.positions, open(str(output_pdb), 'w'), keepIds=True)
        return output_pdb