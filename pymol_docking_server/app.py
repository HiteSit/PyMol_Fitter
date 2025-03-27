from flask import Flask, request, jsonify
import os
import tempfile
from pathlib import Path
import logging
import base64
import json
from typing import Dict, Union, Any

# Import the PyMOL docking module
from pymol_docking_src.Docking_Engine import Pymol_Docking

app = Flask(__name__)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@app.route('/health', methods=['GET'])
def health_check():
    """Simple health check endpoint."""
    return jsonify({"status": "healthy"})

@app.route('/dock', methods=['POST'])
def dock():
    """
    Endpoint to perform docking.
    
    Expected JSON payload:
    {
        "protein": "base64_encoded_pdb_file",
        "ligand": "either_base64_encoded_sdf_or_smiles_string",
        "is_smiles": true/false,
        "dock_mode": "Dock" or "Minimize",
        "output_name": "output_prefix",
        "crystal_sdf": "optional_base64_encoded_sdf_file"
    }
    
    Returns:
    {
        "docked_ligand": "base64_encoded_sdf_file",
        "prepared_protein": "base64_encoded_pdb_file",
        "success": true/false,
        "message": "status message"
    }
    """
    try:
        # Get JSON data from request
        data = request.get_json()
        
        if not data:
            return jsonify({"success": False, "message": "No data provided"}), 400
        
        # Extract parameters
        protein_data = data.get("protein")
        ligand_data = data.get("ligand")
        is_smiles = data.get("is_smiles", False)
        dock_mode = data.get("dock_mode", "Dock")
        output_name = data.get("output_name", "docked")
        crystal_sdf_data = data.get("crystal_sdf")
        
        if not protein_data:
            return jsonify({"success": False, "message": "No protein data provided"}), 400
        if not ligand_data:
            return jsonify({"success": False, "message": "No ligand data provided"}), 400
        
        # Create temporary directory for files
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)
            
            # Save protein to file
            protein_path = temp_dir_path / "protein.pdb"
            if protein_data.startswith("data:"):
                # Handle data URL format
                protein_data = protein_data.split(",")[1]
            with open(protein_path, "wb") as f:
                f.write(base64.b64decode(protein_data))
            
            # Save ligand to file or use SMILES
            if is_smiles:
                ligand_input = ligand_data  # Use SMILES string directly
            else:
                ligand_path = temp_dir_path / "ligand.sdf"
                if ligand_data.startswith("data:"):
                    ligand_data = ligand_data.split(",")[1]
                with open(ligand_path, "wb") as f:
                    f.write(base64.b64decode(ligand_data))
                ligand_input = str(ligand_path)
            
            # Create output paths
            output_path = temp_dir_path / f"{output_name}.sdf"
            
            # Handle crystal SDF if provided
            crystal_sdf_path = None
            if crystal_sdf_data:
                crystal_sdf_path = temp_dir_path / "crystal.sdf"
                if crystal_sdf_data.startswith("data:"):
                    crystal_sdf_data = crystal_sdf_data.split(",")[1]
                with open(crystal_sdf_path, "wb") as f:
                    f.write(base64.b64decode(crystal_sdf_data))
            
            # Initialize Pymol_Docking
            docking = Pymol_Docking(str(protein_path), ligand_input, 
                                  str(crystal_sdf_path) if crystal_sdf_path else None)
            
            # Run docking
            docked_path, protein_prep_path = docking.run_smina_docking(dock_mode, output_name)
            
            # Read result files
            with open(docked_path, "rb") as f:
                docked_ligand_data = base64.b64encode(f.read()).decode("utf-8")
            
            with open(protein_prep_path, "rb") as f:
                prepared_protein_data = base64.b64encode(f.read()).decode("utf-8")
            
            return jsonify({
                "success": True,
                "message": f"Docking completed successfully with mode: {dock_mode}",
                "docked_ligand": docked_ligand_data,
                "prepared_protein": prepared_protein_data
            })
            
    except Exception as e:
        logger.error(f"Error during docking: {e}", exc_info=True)
        return jsonify({"success": False, "message": f"Error: {str(e)}"}), 500

@app.route('/minimize_complex', methods=['POST'])
def minimize_complex():
    """
    Endpoint to perform minimization of a protein-ligand complex.
    
    Expected JSON payload:
    {
        "protein": "base64_encoded_pdb_file",
        "ligand": "base64_encoded_sdf_file",
        "output_name": "output_prefix",
        "crystal_sdf": "optional_base64_encoded_sdf_file"
    }
    
    Returns:
    {
        "minimized_complex": "base64_encoded_pdb_file",
        "success": true/false,
        "message": "status message"
    }
    """
    try:
        # Get JSON data from request
        data = request.get_json()
        
        if not data:
            return jsonify({"success": False, "message": "No data provided"}), 400
        
        # Extract parameters
        protein_data = data.get("protein")
        ligand_data = data.get("ligand")
        output_name = data.get("output_name", "minimized")
        crystal_sdf_data = data.get("crystal_sdf")
        
        if not protein_data:
            return jsonify({"success": False, "message": "No protein data provided"}), 400
        if not ligand_data:
            return jsonify({"success": False, "message": "No ligand data provided"}), 400
        
        # Create temporary directory for files
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)
            
            # Save protein to file
            protein_path = temp_dir_path / "protein.pdb"
            if protein_data.startswith("data:"):
                protein_data = protein_data.split(",")[1]
            with open(protein_path, "wb") as f:
                f.write(base64.b64decode(protein_data))
            
            # Save ligand to file
            ligand_path = temp_dir_path / "ligand.sdf"
            if ligand_data.startswith("data:"):
                ligand_data = ligand_data.split(",")[1]
            with open(ligand_path, "wb") as f:
                f.write(base64.b64decode(ligand_data))
            
            # Handle crystal SDF if provided
            crystal_sdf_path = None
            if crystal_sdf_data:
                crystal_sdf_path = temp_dir_path / "crystal.sdf"
                if crystal_sdf_data.startswith("data:"):
                    crystal_sdf_data = crystal_sdf_data.split(",")[1]
                with open(crystal_sdf_path, "wb") as f:
                    f.write(base64.b64decode(crystal_sdf_data))
            
            # Initialize Pymol_Docking
            docking = Pymol_Docking(str(protein_path), str(ligand_path),
                                   str(crystal_sdf_path) if crystal_sdf_path else None)
            
            # Run complex minimization
            minimized_complex_path = docking.run_complex_minimization(protein_path, ligand_path)
            
            # Read result file
            with open(minimized_complex_path, "rb") as f:
                minimized_complex_data = base64.b64encode(f.read()).decode("utf-8")
            
            return jsonify({
                "success": True,
                "message": "Complex minimization completed successfully",
                "minimized_complex": minimized_complex_data
            })
            
    except Exception as e:
        logger.error(f"Error during complex minimization: {e}", exc_info=True)
        return jsonify({"success": False, "message": f"Error: {str(e)}"}), 500

@app.route('/dock_and_minimize', methods=['POST'])
def dock_and_minimize():
    """
    Endpoint to perform docking followed by minimization in one call.
    
    Expected JSON payload:
    {
        "protein": "base64_encoded_pdb_file",
        "ligand": "either_base64_encoded_sdf_or_smiles_string",
        "is_smiles": true/false,
        "dock_mode": "Dock" or "Minimize",
        "output_name": "output_prefix",
        "crystal_sdf": "optional_base64_encoded_sdf_file"
    }
    
    Returns:
    {
        "docked_ligand": "base64_encoded_sdf_file",
        "prepared_protein": "base64_encoded_pdb_file",
        "minimized_complex": "base64_encoded_pdb_file",
        "success": true/false,
        "message": "status message"
    }
    """
    try:
        # Get JSON data from request
        data = request.get_json()
        
        if not data:
            return jsonify({"success": False, "message": "No data provided"}), 400
        
        # Extract parameters
        protein_data = data.get("protein")
        ligand_data = data.get("ligand")
        is_smiles = data.get("is_smiles", False)
        dock_mode = data.get("dock_mode", "Dock")
        output_name = data.get("output_name", "docked")
        crystal_sdf_data = data.get("crystal_sdf")
        
        if not protein_data:
            return jsonify({"success": False, "message": "No protein data provided"}), 400
        if not ligand_data:
            return jsonify({"success": False, "message": "No ligand data provided"}), 400
        
        # Create temporary directory for files
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)
            
            # Save protein to file
            protein_path = temp_dir_path / "protein.pdb"
            if protein_data.startswith("data:"):
                protein_data = protein_data.split(",")[1]
            with open(protein_path, "wb") as f:
                f.write(base64.b64decode(protein_data))
            
            # Save ligand to file or use SMILES
            if is_smiles:
                ligand_input = ligand_data  # Use SMILES string directly
            else:
                ligand_path = temp_dir_path / "ligand.sdf"
                if ligand_data.startswith("data:"):
                    ligand_data = ligand_data.split(",")[1]
                with open(ligand_path, "wb") as f:
                    f.write(base64.b64decode(ligand_data))
                ligand_input = str(ligand_path)
            
            # Handle crystal SDF if provided
            crystal_sdf_path = None
            if crystal_sdf_data:
                crystal_sdf_path = temp_dir_path / "crystal.sdf"
                if crystal_sdf_data.startswith("data:"):
                    crystal_sdf_data = crystal_sdf_data.split(",")[1]
                with open(crystal_sdf_path, "wb") as f:
                    f.write(base64.b64decode(crystal_sdf_data))
            
            # Initialize Pymol_Docking
            docking = Pymol_Docking(str(protein_path), ligand_input,
                                   str(crystal_sdf_path) if crystal_sdf_path else None)
            
            # Run docking
            docked_path, protein_prep_path = docking.run_smina_docking(dock_mode, output_name)
            
            # Run minimization
            minimized_complex_path = docking.run_complex_minimization(protein_prep_path, docked_path)
            
            # Read result files
            with open(docked_path, "rb") as f:
                docked_ligand_data = base64.b64encode(f.read()).decode("utf-8")
            
            with open(protein_prep_path, "rb") as f:
                prepared_protein_data = base64.b64encode(f.read()).decode("utf-8")
                
            with open(minimized_complex_path, "rb") as f:
                minimized_complex_data = base64.b64encode(f.read()).decode("utf-8")
            
            return jsonify({
                "success": True,
                "message": f"Docking and minimization completed successfully",
                "docked_ligand": docked_ligand_data,
                "prepared_protein": prepared_protein_data,
                "minimized_complex": minimized_complex_data
            })
            
    except Exception as e:
        logger.error(f"Error during docking and minimization: {e}", exc_info=True)
        return jsonify({"success": False, "message": f"Error: {str(e)}"}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True) 