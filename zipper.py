import os
import zipfile
import shutil

def zip_directory(source_dir, output_zip):
    """
    Zip a directory while preserving the directory structure inside the zip file.
    
    Args:
        source_dir (str): Path to the directory to be zipped
        output_zip (str): Path to the output zip file
    """
    # Check if source directory exists
    if not os.path.exists(source_dir):
        print(f"Error: Directory '{source_dir}' does not exist.")
        return False
    
    # Create a temporary directory to work with
    temp_dir = "temp_zip_dir"
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)
    
    # Copy the source directory into the temporary directory
    dir_name = os.path.basename(source_dir)
    shutil.copytree(source_dir, os.path.join(temp_dir, dir_name))
    
    # Create the zip file
    with zipfile.ZipFile(output_zip, 'w', zipfile.ZIP_DEFLATED) as zipf:
        # Walk through the temporary directory
        for root, dirs, files in os.walk(temp_dir):
            for file in files:
                file_path = os.path.join(root, file)
                # Calculate the archive name (path within the zip file)
                arcname = os.path.relpath(file_path, temp_dir)
                zipf.write(file_path, arcname)
    
    # Clean up the temporary directory
    shutil.rmtree(temp_dir)
    print(f"Successfully created zip file: {output_zip}")
    return True

if __name__ == "__main__":
    source_directory = "pymol_fitter_plugin"
    output_zip_file = "pymol_fitter_plugin.zip"
    
    zip_directory(source_directory, output_zip_file)
