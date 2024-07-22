import shutil
from pathlib import Path
from tempfile import gettempdir
import os

import unittest

from Pymol_Docking import Pymol_Docking

class TestPymolDocking(unittest.TestCase):
    def setUp(self):
        self.protein_pdb: Path =  Path("sample_data/LAC3.pdb")
        self.crystal_sdf: Path = Path("sample_data/Crystal.sdf")
        self.ligand_sdf: Path = Path("sample_data/Ligand.sdf")

        self.protein_basename: str = self.protein_pdb.stem
        self.protein_PREP: Path = Path(f"{self.protein_basename}_PREP.pdb")

        self.pymol_docking_class = Pymol_Docking(self.protein_pdb.as_posix(), self.ligand_sdf.as_posix())

    def tearDown(self):
        if self.protein_PREP.exists():
            self.protein_PREP.unlink()

        if Path("XXX.log").exists():
            Path("XXX.log").unlink()

        if Path("XXX.sdf").exists():
            Path("XXX.sdf").unlink()

        if Path("Crystal.sdf").exists():
            Path("Crystal.sdf").unlink()

    def test_prepare_protein(self):
        protein_PREP: Path = self.pymol_docking_class.prepare_protein()
        self.assertTrue(protein_PREP.exists())

    def test_fix3d_mol(self):
        self.pymol_docking_class.fix_3d_mol(self.ligand_sdf, "XXX")

        tmp_path_prep = Path(gettempdir()) / "XXX_fixed.sdf"
        self.assertTrue(tmp_path_prep.exists())

    def test_prepare_ligands(self):
        fixed_ligands, fixed_crystal = self.pymol_docking_class.prepare_ligands()
        self.assertTrue(fixed_ligands.exists())
        self.assertTrue(fixed_crystal.exists())

    def test_run_docking(self):
        shutil.copy(self.crystal_sdf, "Crystal.sdf")

        self.pymol_docking_class.run_docking("Dock", "XXX")
        self.assertTrue(Path("XXX.sdf").exists())
        self.assertTrue(Path("XXX.log").exists())