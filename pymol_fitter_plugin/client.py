import requests
import base64
import json
from pathlib import Path
from typing import Dict, Union, Any, Optional, Tuple, List

class PyMOLDockingClient:
    """
    Client for interacting with the PyMOL Docking Docker API.
    
    This client provides methods to communicate with the PyMOL Docking server
    for molecular docking and minimization operations.
    """
    
    def __init__(self, base_url: str = "http://localhost:5000"):
        """
        Initialize the client with the base URL of the API.
        
        Parameters:
            base_url: The base URL of the API, default is http://localhost:5000
        """
        self.base_url = base_url
        
    def check_health(self) -> bool:
        """
        Check if the API is healthy and accessible.
        
        Returns:
            True if the API is healthy, False otherwise
        """
        try:
            response = requests.get(f"{self.base_url}/health", timeout=5)
            return response.status_code == 200 and response.json().get("status") == "healthy"
        except requests.RequestException:
            return False
            
    def _encode_file(self, file_path: Union[str, Path]) -> str:
        """
        Encode a file as base64.
        
        Parameters:
            file_path: Path to the file to encode
            
        Returns:
            Base64 encoded string of the file
            
        Raises:
            FileNotFoundError: If the file does not exist
            IOError: If there's an error reading the file
        """
        file_path = Path(file_path)
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
            
        try:
            with open(file_path, "rb") as f:
                return base64.b64encode(f.read()).decode("utf-8")
        except IOError as e:
            raise IOError(f"Error reading file {file_path}: {e}") from e
            
    def _decode_and_save(self, base64_data: str, output_path: Union[str, Path]) -> Path:
        """
        Decode a base64 string and save it to a file.
        
        Parameters:
            base64_data: Base64 encoded string
            output_path: Path to save the decoded file
            
        Returns:
            Path to the saved file
            
        Raises:
            ValueError: If the base64 data is invalid
            IOError: If there's an error writing to the file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            with open(output_path, "wb") as f:
                f.write(base64.b64decode(base64_data))
            return output_path
        except base64.binascii.Error as e:
            raise ValueError(f"Invalid base64 data: {e}") from e
        except IOError as e:
            raise IOError(f"Error writing to file {output_path}: {e}") from e

    def dock_minimize(self,
                     protein_file: Union[str, Path],
                     ligand: Union[str, Path, str],
                     crystal_file: Union[str, Path],
                     is_smiles: bool = False,
                     minimize: bool = False,
                     dock_mode: str = "Dock",
                     output_name: str = "docked",
                     output_dir: Union[str, Path] = ".") -> Dict[str, Path]:
        """
        Perform docking and optional minimization of a ligand against a protein.
        
        Parameters:
            protein_file: Path to the protein PDB file
            ligand: Path to the ligand SDF file or SMILES string
            crystal_file: Path to the crystal ligand SDF file for binding site
            is_smiles: Whether the ligand is a SMILES string
            minimize: Whether to perform minimization after docking
            dock_mode: Docking mode, either "Dock" or "Minimize"
            output_name: Base name for the output files
            output_dir: Directory to save the output files
            
        Returns:
            Dictionary with paths to the output files
            
        Raises:
            ConnectionError: If the server is not accessible
            ValueError: If input validation fails
            RuntimeError: If the server returns an error response
        """
        if not self.check_health():
            raise ConnectionError("API is not healthy. Make sure the server is running and accessible.")
        
        # Validate inputs
        if dock_mode not in ["Dock", "Minimize"]:
            raise ValueError("dock_mode must be either 'Dock' or 'Minimize'")
        
        # Prepare the request data
        try:
            if is_smiles:
                # If ligand is a SMILES string, use it directly
                ligand_data = str(ligand)
            else:
                # If ligand is a file path, encode it as base64
                ligand_data = self._encode_file(ligand)
            
            protein_data = self._encode_file(protein_file)
            crystal_data = self._encode_file(crystal_file)
            
            # Prepare the payload
            payload = {
                "protein": protein_data,
                "ligand": ligand_data,
                "crystal": crystal_data,
                "is_smiles": is_smiles,
                "minimize": minimize,
                "dock_mode": dock_mode,
                "output_name": output_name
            }
            
            # Make the API call
            response = requests.post(f"{self.base_url}/dock", json=payload, timeout=3600)
            
            if response.status_code != 200:
                error_message = response.json().get('message', 'Unknown error')
                raise RuntimeError(f"API call failed: {error_message}")
            
            response_data = response.json()
            
            # Save the returned files
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            result = {}
            
            # Save the log
            result["log"] = self._decode_and_save(
                response_data["log"],
                output_dir / f"{output_name}.log"
            )
            
            # Save docked ligand
            result["docked_ligand"] = self._decode_and_save(
                response_data["docked_ligand"],
                output_dir / f"{output_name}.sdf"
            )
            
            # Save prepared protein
            result["prepared_protein"] = self._decode_and_save(
                response_data["prepared_protein"],
                output_dir / f"{output_name}_prepared_protein.pdb"
            )
            
            # Save minimized complex if available
            if minimize and "minimized_complex" in response_data:
                result["minimized_complex"] = self._decode_and_save(
                    response_data["minimized_complex"],
                    output_dir / f"{output_name}_complex.pdb"
                )
            
            return result
            
        except FileNotFoundError as e:
            raise ValueError(f"File not found: {e}") from e
        except requests.RequestException as e:
            raise ConnectionError(f"Network error: {e}") from e
    
    def inplace_minimization(
        self,
        protein_file: Union[str, Path],
        ligand: Optional[Union[str, Path]] = None,
        output_name: str = "minimized_complex",
        output_dir: Union[str, Path] = "."
    ) -> Dict[str, Path]:
        """
        Perform in-place minimization of a protein-ligand complex.
        
        Parameters:
            protein_file: Path to the protein PDB file
            ligand: Path to the ligand SDF file (optional)
            output_name: Name for the output complex file
            output_dir: Directory to save the output files
            
        Returns:
            Dictionary with paths to the output files
            
        Raises:
            ConnectionError: If the server is not accessible
            ValueError: If input validation fails
            RuntimeError: If the server returns an error response
        """
        if not self.check_health():
            raise ConnectionError("API is not healthy. Make sure the server is running and accessible.")
        
        try:
            protein_data = self._encode_file(protein_file)
            
            payload = {
                "protein": protein_data,
                "output_name": output_name
            }
            
            # Add ligand data if provided
            if ligand is not None:
                ligand_data = self._encode_file(ligand)
                payload["ligand"] = ligand_data
            
            response = requests.post(f"{self.base_url}/minimize", json=payload, timeout=3600)
            
            if response.status_code != 200:
                error_message = response.json().get('message', 'Unknown error')
                raise RuntimeError(f"API call failed: {error_message}")
            
            response_data = response.json()
            
            # Save the files
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            result = {}
            
            result["minimized_complex"] = self._decode_and_save(
                response_data["minimized_complex"],
                output_dir / f"{output_name}.pdb"
            )
            
            return result
            
        except FileNotFoundError as e:
            raise ValueError(f"File not found: {e}") from e
        except requests.RequestException as e:
            raise ConnectionError(f"Network error: {e}") from e

    # ───────────────────────── virtual screening ──────────────────────────

    def virtual_screen(
        self,
        protein_file: Union[str, Path],
        ligand: Union[str, Path],
        crystal_file: Union[str, Path],
        is_smiles: bool = True,
        output_name: str = "vs_run",
        output_dir: Union[str, Path] = ".",
    ) -> Dict[str, Path]:
        """Run multi‑ligand virtual screening.

        `ligand` may be a dot‑separated SMILES string (set `is_smiles=True`) or
        the *filename* of a multi‑record SDF residing in the server working
        directory (set `is_smiles=False`).
        """
        if not self.check_health():
            raise ConnectionError("API is not healthy. Make sure the server is running and accessible.")

        # Validate
        if not is_smiles and not isinstance(ligand, (str, Path)):
            raise ValueError("When is_smiles=False the ligand parameter must be a file path or string filename")

        try:
            protein_data = self._encode_file(protein_file)
            crystal_data = self._encode_file(crystal_file)

            payload: Dict[str, Any] = {
                "protein": protein_data,
                "crystal": crystal_data,
                "is_smiles": is_smiles,
                "output_name": output_name,
            }

            if is_smiles:
                payload["ligand"] = str(ligand)
            else:
                # send just the file name; server expects file already present
                payload["ligand"] = str(ligand)

            resp = requests.post(f"{self.base_url}/virtual_screen", json=payload, timeout=7200)
            if resp.status_code != 200:
                err = resp.json().get("message", "Unknown error")
                raise RuntimeError(f"API call failed: {err}")

            data = resp.json()

            out_dir = Path(output_dir)
            out_dir.mkdir(parents=True, exist_ok=True)

            results: Dict[str, Path] = {}

            # Save summary CSV
            summary_path = out_dir / f"{output_name}_summary.csv"
            self._decode_and_save(data["summary_csv"], summary_path)
            results["summary_csv"] = summary_path

            # Save prepared protein
            prep_path = out_dir / f"{output_name}_prepared_protein.pdb"
            self._decode_and_save(data["prepared_protein"], prep_path)
            results["prepared_protein"] = prep_path

            # Save each docked ligand + log (iterate with index)
            lig_paths: List[Path] = []
            for idx, enc in enumerate(data["docked_ligands"], start=1):
                p = out_dir / f"{output_name}_{idx}.sdf"
                self._decode_and_save(enc, p)
                lig_paths.append(p)
            results["docked_ligands"] = lig_paths

            log_paths: List[Path] = []
            for idx, enc in enumerate(data["logs"], start=1):
                p = out_dir / f"{output_name}_{idx}.log"
                self._decode_and_save(enc, p)
                log_paths.append(p)
            results["logs"] = log_paths

            return results

        except FileNotFoundError as e:
            raise ValueError(f"File not found: {e}") from e
        except requests.RequestException as e:
            raise ConnectionError(f"Network error: {e}") from e