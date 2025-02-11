from pathlib import Path
from typing import List, Optional, Union, Dict

import CDPL.Chem as Chem
import CDPL.ConfGen as ConfGen

import datamol as dm
import rdkit.Chem as rdChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

def molsToFiles(mols: list[Chem.BasicMolecule], output_file: Path) -> None:
    """
    Writes a list of molecules to the specified output file.

    Parameters:
    - mols (list): List of CDPKit molecules to write to the output file.
    - output_file (str): Path to the output file.
    """
    # Create a writer for the output molecules (format specified by file extension)
    writer = Chem.MolecularGraphWriter(str(output_file))

    for mol in mols:
        writer.write(mol)
    
    writer.close()
    
def filesToMols(sdf_inut: Path) -> List[Chem.BasicMolecule]:
    """
    Retrieves the structure data of each molecule in the provided SD file and outputs it to the console.

    Parameters:
    - sdf_input (Path): Path to the input SD file.
    """
    # Create reader for MDL SD-files
    reader = Chem.FileSDFMoleculeReader(str(sdf_inut))

    # create an instance of the default implementation of the Chem.Molecule interface
    mol = Chem.BasicMolecule()
    
    # Store the BasicMolecule object in a list
    mols_list = []
    # Iterate over each molecule in the file and retrieve structure data
    try:
        while reader.read(mol):
            try:
                if not Chem.hasStructureData(mol):
                    print('Error: no structure data available for molecule', Chem.getName(mol))
                    continue

                struct_data = Chem.getStructureData(mol)
                new_mol = Chem.BasicMolecule(mol)
                mols_list.append(new_mol)
            except Exception as e:
                pass
            
    except Exception as e:
        pass
    
    return mols_list

def standardize(chembl_proc: Chem.ChEMBLStandardizer,
                in_mol: Chem.Molecule, out_mol: Chem.Molecule,
                proc_excluded: bool, extract_parent: bool) -> Chem.ChEMBLStandardizer.ChangeFlags:
    """
    Performs ChEMBL molecule standardization and parent structure extraction (optional) for a given input molecule using a provided Chem.ChEMBLStandardizer instance.

    Parameters:
    - chembl_proc (Chem.ChEMBLStandardizer): Instance of the Chem.ChEMBLStandardizer class.
    - in_mol (Chem.Molecule): Input molecule to standardize.
    - out_mol (Chem.Molecule): Output molecule to store the standardized molecule.
    - proc_excluded (bool): If True, molecules flagged as excluded will be processed.
    - extract_parent (bool): If True, the parent structure will be extracted.

    Returns:
    - Chem.ChEMBLStandardizer.ChangeFlags: Flags indicating the carried out modifications.
    """
    # here, the standardization is carried out on a copy of the read input molecule
    # (if only one molecule instance gets provided as argument, modifications will be made in-place)
    change_flags = chembl_proc.standardize(in_mol, out_mol, proc_excluded)

    if extract_parent: # perform parent structure extraction (optional)
        change_flags &= ~Chem.ChEMBLStandardizer.EXCLUDED  # clear excluded flag possibly set by the standardization
                                                       # procedure (might change after salt stripping)
        change_flags |= chembl_proc.getParent(out_mol)     # extract parent structure (in-place) and add information
                                                       # about the carried out modifcations
    return change_flags

def getListOfChangesString(change_flags: Chem.ChEMBLStandardizer.ChangeFlags, verbose: bool = False) -> str:
    """
    Returns a string listing the carried out modifications.

    Parameters:
    - change_flags (Chem.ChEMBLStandardizer.ChangeFlags): Flags indicating the carried out modifications.
    - verbose (bool): If True, the string will contain a detailed list of the carried out modifications.

    Returns:
    - str: String listing the carried out modifications.
    """
    if not verbose:
        return None

    changes = '   Carried out modifications:'

    # List of possible changes
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

def getLogMessage(change_flags: Chem.ChEMBLStandardizer.ChangeFlags,
                  proc_excluded: bool, extract_parent: bool, mol_id: str, verbose: bool = False) -> str:
    """
    Returns a log message describing the carried out modifications.

    Parameters:
    - change_flags (Chem.ChEMBLStandardizer.ChangeFlags): Flags indicating the carried out modifications.
    - proc_excluded (bool): If True, molecules flagged as excluded will be processed.
    - extract_parent (bool): If True, the parent structure will be extracted.
    - mol_id (str): Identifier of the molecule.

    Returns:
    - str: Log message describing the carried out modifications.
    """
    if (change_flags & Chem.ChEMBLStandardizer.EXCLUDED) and proc_excluded:
        return f'Molecule {mol_id}: discarded (flagged as excluded)'

    if not proc_excluded and (change_flags & Chem.ChEMBLStandardizer.EXCLUDED):
        return f'Molecule {mol_id}: forwarded unchanged (flagged as excluded)'

    if change_flags:
        return f'Molecule {mol_id}: modified\n{getListOfChangesString(change_flags, verbose)}'

    return f'Molecule {mol_id}: forwarded unchanged'

def generate3dConformation(mol: Chem.Molecule, struct_gen: ConfGen.StructureGenerator) -> int:
    """
    Generates a low-energy 3D structure of the argument molecule using the provided initialized ConfGen.StructureGenerator instance.

    Parameters:
    - mol (Chem.Molecule): Molecule to generate a 3D structure for.
    - struct_gen (ConfGen.StructureGenerator): Instance of the ConfGen.StructureGenerator class.

    Returns:
    - int: Status code indicating the success of the 3D structure generation.
    """
    # prepare the molecule for 3D structure generation
    ConfGen.prepareForConformerGeneration(mol)

    # generate the 3D structure
    status = struct_gen.generate(mol)

    # if sucessful, store the generated conformer ensemble as
    # per atom 3D coordinates arrays (= the way conformers are represented in CDPKit)
    if status == ConfGen.ReturnCode.SUCCESS:
        struct_gen.setCoordinates(mol)

    # return status code
    return status

def stero_enumerator(input_sdf: Path, output_sdf: Path) -> Path:
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
    def __init__(self,
                 standardize: bool = True,
                 protonate: bool = True,
                 gen3d: bool = True,
                 ):
        
        self.chembl_proc = Chem.ChEMBLStandardizer()
        self.prot_state_gen = Chem.ProtonationStateStandardizer()
        
        self.standardize = standardize
        self.protonate = protonate
        self.gen3d = gen3d
        
    def prepare_ligands(self, input_sdf: Path, output_sdf: Path):
        
        mols: List[Chem.BasicMolecule] = filesToMols(input_sdf)
        writer = Chem.MolecularGraphWriter(str(output_sdf))
        
        for mol in mols:
            in_mol = mol
            out_mol = Chem.BasicMolecule()
            
            # Standardization
            if self.standardize:
                change_flags = standardize(self.chembl_proc, in_mol, out_mol, True, True)
                change_flags_str = getListOfChangesString(change_flags, verbose=True)
            
            # Protonation
            if self.protonate:
                self.prot_state_gen.standardize(out_mol, Chem.ProtonationStateStandardizer.PHYSIOLOGICAL_CONDITION_STATE)
                Chem.perceiveComponents(out_mol, True)

            if self.gen3d:
                Chem.setMultiConfExportParameter(writer, False)
                struct_gen = ConfGen.StructureGenerator()
                status_to_str = { ConfGen.ReturnCode.UNINITIALIZED                  : 'uninitialized',
                                ConfGen.ReturnCode.TIMEOUT                                          : 'max. processing time exceeded',
                                ConfGen.ReturnCode.ABORTED                                          : 'aborted',
                                ConfGen.ReturnCode.FORCEFIELD_SETUP_FAILED                          : 'force field setup failed',
                                ConfGen.ReturnCode.FORCEFIELD_MINIMIZATION_FAILED                   : 'force field structure refinement failed',
                                ConfGen.ReturnCode.FRAGMENT_LIBRARY_NOT_SET                         : 'fragment library not available',
                                ConfGen.ReturnCode.FRAGMENT_CONF_GEN_FAILED                         : 'fragment conformer generation failed',
                                ConfGen.ReturnCode.FRAGMENT_CONF_GEN_TIMEOUT                        : 'fragment conformer generation timeout',
                                ConfGen.ReturnCode.FRAGMENT_ALREADY_PROCESSED                       : 'fragment already processed',
                                ConfGen.ReturnCode.TORSION_DRIVING_FAILED                           : 'torsion driving failed',
                                ConfGen.ReturnCode.CONF_GEN_FAILED                                  : 'conformer generation failed' }
                
                status = generate3dConformation(out_mol, struct_gen)
                if status == ConfGen.ReturnCode.SUCCESS:
                    Chem.setMDLDimensionality(out_mol, 3)
        
            writer.write(out_mol)
        
        writer.close()