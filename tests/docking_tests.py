import shutil
from pathlib import Path
from tempfile import gettempdir
import os

import unittest

from Pymol_Docking import Pymol_Docking
from Pymol_Docking import openeye_gen3d, openeye_fixer

class TestPymolDocking_SDF(unittest.TestCase):
    def setUp(self):
        self.protein_pdb: Path =  Path("sample_data/LAC3.pdb")
        self.crystal_sdf: Path = Path("sample_data/Crystal.sdf")
        shutil.copy(self.crystal_sdf, "Crystal.sdf")

        self.input_ligands: Path = Path("sample_data/Ligand.sdf")

        self.protein_basename: str = self.protein_pdb.stem
        self.protein_PREP: Path = Path(f"{self.protein_basename}_PREP.pdb")

        self.pymol_docking_class = Pymol_Docking(self.protein_pdb.as_posix(), self.input_ligands.as_posix())

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
        self.pymol_docking_class.fix_3d_mol(self.input_ligands, "XXX")

        tmp_path_prep = Path(gettempdir()) / "XXX_fixed.sdf"
        self.assertTrue(tmp_path_prep.exists())

    def test_prepare_ligands(self):
        fixed_ligands, fixed_crystal = self.pymol_docking_class.prepare_ligands()
        self.assertTrue(fixed_ligands.exists())
        self.assertTrue(fixed_crystal.exists())

    def test_run_docking(self):
        self.pymol_docking_class.run_smina_docking("Dock", "XXX")
        self.assertTrue(Path("XXX.sdf").exists())
        self.assertTrue(Path("XXX.log").exists())

class TestPymolDocking_SMILE(unittest.TestCase):
    def setUp(self):
        self.protein_pdb: Path =  Path("sample_data/LAC3.pdb")
        self.crystal_sdf: Path = Path("sample_data/Crystal.sdf")
        shutil.copy(self.crystal_sdf, "Crystal.sdf")

        self.input_ligands: str = "CC1=C2NC(=O)C(CC(=O)NC(CC3=CC=CC=C3)C3=CC=CC=C3)OC2=CC=C1"

        self.protein_basename: str = self.protein_pdb.stem
        self.protein_PREP: Path = Path(f"{self.protein_basename}_PREP.pdb")

        self.pymol_docking_class = Pymol_Docking(self.protein_pdb.as_posix(), self.input_ligands)

    def tearDown(self):
        if self.protein_PREP.exists():
            self.protein_PREP.unlink()

        if Path("XXX.log").exists():
            Path("XXX.log").unlink()

        if Path("XXX.sdf").exists():
            Path("XXX.sdf").unlink()

        if Path("Crystal.sdf").exists():
            Path("Crystal.sdf").unlink()

    def test_openeye_gen3d(self):
        mypath = openeye_gen3d(self.input_ligands)
        self.assertTrue(mypath.exists())

    def test_prepare_ligands(self):
        fixed_ligands, fixed_crystal = self.pymol_docking_class.prepare_ligands()
        self.assertTrue(fixed_ligands.exists())
        self.assertTrue(fixed_crystal.exists())

    def test_run_docking(self):
        self.pymol_docking_class.run_smina_docking("Dock", "XXX")
        self.assertTrue(Path("XXX.sdf").exists())
        size = os.path.getsize("XXX.sdf")
        self.assertGreater(size, 0)

        self.assertTrue(Path("XXX.log").exists())

    def test_off_site_docking(self):
        from Pymol_Docking import off_site_docking
        from pymol import cmd
        cmd.reinitialize()

        protein_PREP: Path = self.pymol_docking_class.prepare_protein()
        cmd.load(self.protein_PREP.as_posix(), "protein_")
        off_site_docking("protein_", "[H]N([C@H](C)C1=CC=CC=C1C)C1=CN=C2N=C(Cl)C3=C(N([H])C=C3)N12", "XXX")

        size = os.path.getsize("XXX.sdf")
        self.assertTrue(Path("XXX.sdf").exists())
        self.assertGreater(size, 0)