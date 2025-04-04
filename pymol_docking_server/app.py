from flask import Flask, request, jsonify
import os
import tempfile
from pathlib import Path
import logging
import base64
import json
from typing import Dict, Union, Any

# Import the PyMOL docking module
from pymol_docking_src.Docking_Engine import Pymol_Docking, outer_minimization

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

@app.route('/dock', methods=["POST"])
def dock_minimize():
    try:
        data = request.get_json()
            
        protein_data = data.get("protein")
        ligand_data = data.get("ligand")
        crystal_data = data.get("crystal")
        is_smiles = data.get("is_smiles", False)
        minimization_bool = data.get("minimize", False)
        dock_mode = data.get("dock_mode", "Dock")
        output_name = data.get("output_name", "docked")
        
        if not protein_data:
            return jsonify({"success": False, "message": "No protein data provided"}), 400
        if not ligand_data:
            return jsonify({"success": False, "message": "No ligand data provided"}), 400
        if not crystal_data:
            return jsonify({"success": False, "message": "No crystal data provided"}), 400
        
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
            
            # Save the Crystal
            crystal_sdf_path = temp_dir_path / "crystal.sdf"
            if crystal_data.startswith("data:"):
                crystal_data = crystal_data.split(",")[1]
            with open(crystal_sdf_path, "wb") as f:
                f.write(base64.b64decode(crystal_data))
                
            # Initialize Pymol_Docking
            docking_class = Pymol_Docking(str(protein_path), ligand_input, str(crystal_sdf_path), is_smiles)
            docked_path, protein_prep_path, docked_log = docking_class.run_smina_docking(dock_mode, output_name)
            
            if docked_path is None or protein_prep_path is None:
                return jsonify({
                    "success": False,
                    "message": "Docking failed: most probably there is a problem in the kekulization of the ligand. Try again from the SMILES.",
                    "log": "Critical: No docking results found, smina output file was empty"
                }), 500
            
            # Encode results
            with open(docked_path, "rb") as f:
                docked_ligand_data = base64.b64encode(f.read()).decode("utf-8")
            
            with open(protein_prep_path, "rb") as f:
                prepared_protein_data = base64.b64encode(f.read()).decode("utf-8")
            
            with open(docked_log, "rb") as f:
                docked_log_content = base64.b64encode(f.read()).decode("utf-8")
            
            if minimization_bool:
                protein_complex = docking_class.run_complex_minimization(protein_path, docked_path)
                with open(protein_complex, "rb") as f:
                    protein_complex_data = base64.b64encode(f.read()).decode("utf-8")
                    
                return jsonify({
                    "success": True,
                    "message": f"Docking completed successfully with mode: {dock_mode}",
                    "log": docked_log_content,
                    "docked_ligand": docked_ligand_data,
                    "prepared_protein": prepared_protein_data,
                    "minimized_complex": protein_complex_data
                })
                
            else:
                return jsonify({
                    "success": True,
                    "message": f"Docking completed successfully with mode: {dock_mode}",
                    "log": docked_log_content,
                    "docked_ligand": docked_ligand_data,
                    "prepared_protein": prepared_protein_data,
                })
                
    except Exception as e:
        logger.error(f"Error during docking: {e}", exc_info=True)
        return jsonify({"success": False, "message": f"Error: {str(e)}"}), 500

@app.route('/minimize', methods=['POST'])
def minimize():
    try:
        data = request.get_json()
        
        protein_data = data.get("protein")
        ligand_data = data.get("ligand")
        output_name = data.get("output_name", "minimized_complex")
        
        if not protein_data:
            return jsonify({"success": False, "message": "No protein data provided"}), 400
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)
            
            # Protein
            protein_path = temp_dir_path / "protein.pdb"
            if protein_data.startswith("data:"):
                protein_data = protein_data.split(",")[1]
            with open(protein_path, "wb") as f:
                f.write(base64.b64decode(protein_data))
            
            # Ligand (optional)
            ligand_path = None
            if ligand_data:
                ligand_path = temp_dir_path / "ligand.sdf"
                if ligand_data.startswith("data:"):
                    ligand_data = ligand_data.split(",")[1]
                with open(ligand_path, "wb") as f:
                    f.write(base64.b64decode(ligand_data))
            
            complex_path = temp_dir_path / f"{output_name}.pdb"
            outer_minimization(protein_path, ligand_path, complex_path)
            
            with open(complex_path, "rb") as f:
                complex_data = base64.b64encode(f.read()).decode("utf-8")
                
            return jsonify({
                "success": True,
                "message": "Minimization completed successfully",
                "minimized_complex": complex_data
            })
            
    except Exception as e:
        logger.error(f"Error during minimization: {e}", exc_info=True)
        return jsonify({"success": False, "message": f"Error: {str(e)}"}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True) 