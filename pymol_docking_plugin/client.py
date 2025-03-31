import requests
import base64
import json
from pathlib import Path
from typing import Dict, Union, Any, Optional, Tuple

class PyMOLDockingClient:
    """
    Client for interacting with the PyMOL Docking Docker API.
    """
    
    def __init__(self, base_url: str = "http://localhost:5000"):
        """
        Initialize the client with the base URL of the API.
        
        Parameters:
        -----------
        base_url : str
            The base URL of the API, default is http://localhost:5000
        """
        self.base_url = base_url
        
    def check_health(self) -> bool:
        """
        Check if the API is healthy.
        
        Returns:
        --------
        bool
            True if the API is healthy, False otherwise
        """
        try:
            response = requests.get(f"{self.base_url}/health")
            return response.status_code == 200 and response.json().get("status") == "healthy"
        except Exception:
            return False
            
    def _encode_file(self, file_path: Union[str, Path]) -> str:
        """
        Encode a file as base64.
        
        Parameters:
        -----------
        file_path : Union[str, Path]
            Path to the file to encode
            
        Returns:
        --------
        str
            Base64 encoded string of the file
        """
        with open(file_path, "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")
            
    def _decode_and_save(self, base64_data: str, output_path: Union[str, Path]) -> Path:
        """
        Decode a base64 string and save it to a file.
        
        Parameters:
        -----------
        base64_data : str
            Base64 encoded string
        output_path : Union[str, Path]
            Path to save the decoded file
            
        Returns:
        --------
        Path
            Path to the saved file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, "wb") as f:
            f.write(base64.b64decode(base64_data))
            
        return output_path

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
        -----------
        protein_file : Union[str, Path]
            Path to the protein PDB file
        ligand : Union[str, Path, str]
            Path to the ligand SDF file or SMILES string
        crystal_file : Union[str, Path]
            Path to the crystal ligand SDF file for defining the binding site
        is_smiles : bool
            Whether the ligand is a SMILES string, default is False
        minimize : bool
            Whether to perform minimization after docking, default is False
        dock_mode : str
            Docking mode, either "Dock" or "Minimize", default is "Dock"
        output_name : str
            Base name for the output files, default is "docked"
        output_dir : Union[str, Path]
            Directory to save the output files, default is current directory
            
        Returns:
        --------
        Dict[str, Path]
            Dictionary with paths to the output files:
            - "docked_ligand": Path to the docked ligand SDF file
            - "prepared_protein": Path to the prepared protein PDB file
            - "minimized_complex": Path to the minimized complex PDB file (if minimize=True)
        """
        if not self.check_health():
            raise ConnectionError("API is not healthy. Make sure the server is running and accessible.")
        
        # Prepare the request data
        if is_smiles:
            # If ligand is a SMILES string, use it directly
            ligand_data = ligand
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
        response = requests.post(f"{self.base_url}/dock", json=payload)
        
        if response.status_code != 200:
            raise Exception(f"API call failed: {response.json().get('message', 'Unknown error')}")
        
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
    
    # def dock(self, 
    #          protein_file: Union[str, Path], 
    #          ligand: Union[str, Path, str], 
    #          is_smiles: bool = False,
    #          dock_mode: str = "Dock",
    #          output_name: str = "docked",
    #          output_dir: Union[str, Path] = ".",
    #          crystal_file: Optional[Union[str, Path]] = None) -> Dict[str, Path]:
    #     """
    #     Perform docking of a ligand against a protein.
        
    #     Parameters:
    #     -----------
    #     protein_file : Union[str, Path]
    #         Path to the protein PDB file
    #     ligand : Union[str, Path, str]
    #         Path to the ligand SDF file or SMILES string
    #     is_smiles : bool
    #         Whether the ligand is a SMILES string, default is False
    #     dock_mode : str
    #         Docking mode, either "Dock" or "Minimize", default is "Dock"
    #     output_name : str
    #         Base name for the output files, default is "docked"
    #     output_dir : Union[str, Path]
    #         Directory to save the output files, default is current directory
    #     crystal_file : Optional[Union[str, Path]]
    #         Path to the crystal ligand SDF file for defining the binding site, default is None
            
    #     Returns:
    #     --------
    #     Dict[str, Path]
    #         Dictionary with paths to the output files:
    #         - "docked_ligand": Path to the docked ligand SDF file
    #         - "prepared_protein": Path to the prepared protein PDB file
    #     """
    #     if not self.check_health():
    #         raise ConnectionError("API is not healthy. Make sure the server is running and accessible.")
        
    #     # Prepare the request data
    #     if is_smiles:
    #         # If ligand is a SMILES string, use it directly
    #         ligand_data = ligand
    #     else:
    #         # If ligand is a file path, encode it as base64
    #         ligand_data = self._encode_file(ligand)
        
    #     protein_data = self._encode_file(protein_file)
        
    #     # Prepare the payload
    #     payload = {
    #         "protein": protein_data,
    #         "ligand": ligand_data,
    #         "is_smiles": is_smiles,
    #         "dock_mode": dock_mode,
    #         "output_name": output_name
    #     }
        
    #     # Add crystal data if provided
    #     if crystal_file:
    #         payload["crystal_sdf"] = self._encode_file(crystal_file)
        
    #     # Make the API call
    #     response = requests.post(f"{self.base_url}/dock", json=payload)
        
    #     if response.status_code != 200:
    #         raise Exception(f"API call failed: {response.json().get('message', 'Unknown error')}")
        
    #     response_data = response.json()
        
    #     # Save the returned files
    #     output_dir = Path(output_dir)
    #     output_dir.mkdir(parents=True, exist_ok=True)
        
    #     docked_ligand_path = self._decode_and_save(
    #         response_data["docked_ligand"], 
    #         output_dir / f"{output_name}.sdf"
    #     )
        
    #     prepared_protein_path = self._decode_and_save(
    #         response_data["prepared_protein"], 
    #         output_dir / f"{output_name}_prepared_protein.pdb"
    #     )
        
    #     return {
    #         "docked_ligand": docked_ligand_path,
    #         "prepared_protein": prepared_protein_path
    #     }
    
    # def minimize_complex(self, 
    #                     protein_file: Union[str, Path], 
    #                     ligand_file: Union[str, Path],
    #                     output_name: str = "minimized",
    #                     output_dir: Union[str, Path] = ".",
    #                     crystal_file: Optional[Union[str, Path]] = None) -> Dict[str, Path]:
    #     """
    #     Perform minimization of a protein-ligand complex.
        
    #     Parameters:
    #     -----------
    #     protein_file : Union[str, Path]
    #         Path to the protein PDB file
    #     ligand_file : Union[str, Path]
    #         Path to the ligand SDF file
    #     output_name : str
    #         Base name for the output files, default is "minimized"
    #     output_dir : Union[str, Path]
    #         Directory to save the output files, default is current directory
    #     crystal_file : Optional[Union[str, Path]]
    #         Path to the crystal ligand SDF file for reference, default is None
            
    #     Returns:
    #     --------
    #     Dict[str, Path]
    #         Dictionary with paths to the output files:
    #         - "minimized_complex": Path to the minimized complex PDB file
    #     """
    #     if not self.check_health():
    #         raise ConnectionError("API is not healthy. Make sure the server is running and accessible.")
        
    #     # Prepare the request data
    #     protein_data = self._encode_file(protein_file)
    #     ligand_data = self._encode_file(ligand_file)
        
    #     # Prepare the payload
    #     payload = {
    #         "protein": protein_data,
    #         "ligand": ligand_data,
    #         "output_name": output_name
    #     }
        
    #     # Add crystal data if provided
    #     if crystal_file:
    #         payload["crystal_sdf"] = self._encode_file(crystal_file)
        
    #     # Make the API call
    #     response = requests.post(f"{self.base_url}/minimize_complex", json=payload)
        
    #     if response.status_code != 200:
    #         raise Exception(f"API call failed: {response.json().get('message', 'Unknown error')}")
        
    #     response_data = response.json()
        
    #     # Save the returned files
    #     output_dir = Path(output_dir)
    #     output_dir.mkdir(parents=True, exist_ok=True)
        
    #     minimized_complex_path = self._decode_and_save(
    #         response_data["minimized_complex"], 
    #         output_dir / f"{output_name}_complex.pdb"
    #     )
        
    #     return {
    #         "minimized_complex": minimized_complex_path
    #     }
    
    # def dock_and_minimize(self, 
    #                      protein_file: Union[str, Path], 
    #                      ligand: Union[str, Path, str], 
    #                      is_smiles: bool = False,
    #                      dock_mode: str = "Dock",
    #                      output_name: str = "docked",
    #                      output_dir: Union[str, Path] = ".",
    #                      crystal_file: Optional[Union[str, Path]] = None) -> Dict[str, Path]:
    #     """
    #     Perform docking of a ligand against a protein followed by minimization of the complex.
        
    #     Parameters:
    #     -----------
    #     protein_file : Union[str, Path]
    #         Path to the protein PDB file
    #     ligand : Union[str, Path, str]
    #         Path to the ligand SDF file or SMILES string
    #     is_smiles : bool
    #         Whether the ligand is a SMILES string, default is False
    #     dock_mode : str
    #         Docking mode, either "Dock" or "Minimize", default is "Dock"
    #     output_name : str
    #         Base name for the output files, default is "docked"
    #     output_dir : Union[str, Path]
    #         Directory to save the output files, default is current directory
    #     crystal_file : Optional[Union[str, Path]]
    #         Path to the crystal ligand SDF file for defining the binding site, default is None
            
    #     Returns:
    #     --------
    #     Dict[str, Path]
    #         Dictionary with paths to the output files:
    #         - "docked_ligand": Path to the docked ligand SDF file
    #         - "prepared_protein": Path to the prepared protein PDB file
    #         - "minimized_complex": Path to the minimized complex PDB file
    #     """
    #     if not self.check_health():
    #         raise ConnectionError("API is not healthy. Make sure the server is running and accessible.")
        
    #     # Prepare the request data
    #     if is_smiles:
    #         # If ligand is a SMILES string, use it directly
    #         ligand_data = ligand
    #     else:
    #         # If ligand is a file path, encode it as base64
    #         ligand_data = self._encode_file(ligand)
        
    #     protein_data = self._encode_file(protein_file)
        
    #     # Prepare the payload
    #     payload = {
    #         "protein": protein_data,
    #         "ligand": ligand_data,
    #         "is_smiles": is_smiles,
    #         "dock_mode": dock_mode,
    #         "output_name": output_name
    #     }
        
    #     # Add crystal data if provided
    #     if crystal_file:
    #         payload["crystal_sdf"] = self._encode_file(crystal_file)
        
    #     # Make the API call
    #     response = requests.post(f"{self.base_url}/dock_and_minimize", json=payload)
        
    #     if response.status_code != 200:
    #         raise Exception(f"API call failed: {response.json().get('message', 'Unknown error')}")
        
    #     response_data = response.json()
        
    #     # Save the returned files
    #     output_dir = Path(output_dir)
    #     output_dir.mkdir(parents=True, exist_ok=True)
        
    #     docked_ligand_path = self._decode_and_save(
    #         response_data["docked_ligand"], 
    #         output_dir / f"{output_name}.sdf"
    #     )
        
    #     prepared_protein_path = self._decode_and_save(
    #         response_data["prepared_protein"], 
    #         output_dir / f"{output_name}_prepared_protein.pdb"
    #     )
        
    #     minimized_complex_path = self._decode_and_save(
    #         response_data["minimized_complex"], 
    #         output_dir / f"{output_name}_complex.pdb"
    #     )
        
    #     return {
    #         "docked_ligand": docked_ligand_path,
    #         "prepared_protein": prepared_protein_path,
    #         "minimized_complex": minimized_complex_path
    #     } 