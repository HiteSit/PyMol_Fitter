import shutil
from pathlib import Path
from tempfile import gettempdir
import os

import unittest

from Pymol_Docking import Plants_Docking
from Pymol_Docking import openeye_gen3d, openeye_fixer

class TestPymolDocking_SDF(unittest.TestCase):
    def setUp(self):
        self.protein_pdb: Path =  Path("sample_data/LAC3.pdb")
        self.crystal_sdf: Path = Path("sample_data/Crystal.sdf")
        shutil.copy(self.crystal_sdf, "Crystal.sdf")

        self.input_ligands: Path = Path("sample_data/Ligand.sdf")

        self.plants_class_sdf = Plants_Docking(self.protein_pdb, self.input_ligands)
        self.plants_class_smile = Plants_Docking(self.protein_pdb, "C1=CC=CC=C1")

    def tearDown(self):
        plants_raw = Path("plants_raw")
        if plants_raw.exists():
            shutil.rmtree(plants_raw)

    def test_prepare_protein(self):
        protein_PREP: Path = self.plants_class_sdf.prepare_protein()
        binding_site = self.plants_class_sdf._define_binding_site()

        self.assertTrue(protein_PREP.exists())

    def test_prepare_ligands(self):
        fixed_ligand, fixed_crystal = self.plants_class_smile.prepare_ligands()
        self.assertTrue(fixed_ligand.exists())
        self.assertTrue(fixed_crystal.exists())

        fixed_ligand, fixed_crystal = self.plants_class_sdf.prepare_ligands()
        self.assertTrue(fixed_ligand.exists())
        self.assertTrue(fixed_crystal.exists())

    def test_write_conf(self):
        protein_PREP: Path = self.plants_class_smile.prepare_protein()
        fixed_ligand, fixed_crystal = self.plants_class_smile.prepare_ligands()
        config_str = self.plants_class_smile.write_conf(fixed_ligand, 10)

    def test_run_plants_docking(self):
        plants_raw = Path("plants_raw")
        docked_ligand:Path = plants_raw / "docked_ligands.mol2"

        self.plants_class_sdf.run_plants_docking(10)
        self.assertTrue(docked_ligand.exists())
