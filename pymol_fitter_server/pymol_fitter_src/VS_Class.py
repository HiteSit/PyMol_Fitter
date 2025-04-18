"""Re‑engineered `Pymol_Docking` class
===================================
Provides **single‑ligand docking/minimisation** (legacy behaviour) *and*
a dedicated **virtual‑screening** workflow triggered at construction
(`virtual_screening=True`).  The VS method enforces docking‑only (no
post‑dynamics minimisation), executes independent smina runs in
parallel, and writes a CSV summary keyed by the authoritative
`smina` score tag **`minimizedAffinity`**.

Copy this file over the old class in
`pymol_fitter_server/pymol_fitter_src/Docking_Engine.py`.
"""

from __future__ import annotations

# ─────────────────────────────── imports ────────────────────────────────
import csv
import logging
import multiprocessing as mp
import os
import subprocess
from functools import partial
from pathlib import Path
from tempfile import NamedTemporaryFile, gettempdir
from typing import Dict, List, Literal, Optional, Tuple

import datamol as dm
from rdkit import Chem

from .Protein_Preparation import ProteinPreparation_Protoss, ProteinPreparation_PDBFixer
from .CDPK_Utils import CDPK_Runner
from .Protein_Minimization import minimize_complex

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# authoritative SDF tag used by smina for affinity
_SMINA_TAG = "minimizedAffinity"  # cf. docs & community threads


class Pymol_Docking:
    """Protein–ligand docking utility with optional **virtual screening**."""

    # ───────────────────────── initialisation ──────────────────────────

    def __init__(
        self,
        protein_pdb: str,
        input_ligands: str,
        crystal_sdf: Optional[str] = None,
        *,
        is_smiles: bool = False,
        virtual_screening: bool = False,
    ) -> None:
        self.workdir = Path(os.getcwd())
        self.protein_pdb = Path(protein_pdb)
        self.crystal_sdf = Path(crystal_sdf) if crystal_sdf else Path("Crystal.sdf")
        self.virtual_screening = virtual_screening

        if is_smiles:
            self.input_mode: Literal["SMILES", "SDF"] = "SMILES"
            self.ligands_smiles = input_ligands
        else:
            self.input_mode = "SDF"
            self.ligands_sdf = Path(input_ligands)

        self._prepared_protein: Optional[Path] = None

    # ──────────────────────── protein preparation ──────────────────────

    @staticmethod
    def _filter_protein_only_aa(pdb_file: Path) -> Path:
        """Strip heteroatoms; keep amino‑acids only."""
        import biotite.structure as struc
        from biotite.structure.io import pdb

        struct = pdb.PDBFile.read(pdb_file).get_structure(model=1)
        struct = struct[struc.filter_amino_acids(struct)]
        out = NamedTemporaryFile(suffix=".pdb", delete=False)
        pdb.PDBFile.set_structure(pdb.PDBFile(), struct)
        pdb.PDBFile().write(out.name)
        return Path(out.name)

    def prepare_protein(self) -> Path:
        """Prepare receptor once; cache result.

        * **Protoss** is first tried on the *original* PDB/CIF.  If the
          cloud service is unreachable or returns an invalid job (the
          `output_protein` field can be `null`), we fall back to the
          local **PDBFixer** pathway.
        * A minimal "filter to amino‑acids only" is applied **only** to
          the PDBFixer branch (Protoss expects the full structure).
        """
        if self._prepared_protein and self._prepared_protein.exists():
            return self._prepared_protein

        src = self.protein_pdb.resolve()
        out_path = self.workdir / f"{src.stem}_Prep.pdb"

        # 1) try remote Protoss on the untouched file
        try:
            pp = ProteinPreparation_Protoss()
            self._prepared_protein = pp(src, out_path)
            logger.info("Protein prepared via Protoss")
            return self._prepared_protein
        except Exception as e:
            logger.warning(f"Protoss unavailable or failed ({e}); switching to PDBFixer")

        # 2) local PDBFixer after stripping heteroatoms (safer)
        try:
            filtered = self._filter_protein_only_aa(src)
            self._prepared_protein = ProteinPreparation_PDBFixer()(filtered, out_path)
            logger.info("Protein prepared via PDBFixer")
            return self._prepared_protein
        except Exception as e2:
            raise RuntimeError(f"Protein preparation failed using both Protoss and PDBFixer: {e2}")


        filtered = self._filter_protein_only_aa(self.protein_pdb)
        out_path = self.workdir / f"{self.protein_pdb.stem}_Prep.pdb"
        try:
            self._prepared_protein = ProteinPreparation_Protoss()(filtered, out_path)
        except Exception as e:
            logger.warning(f"Protoss failed ({e}); using PDBFixer")
            self._prepared_protein = ProteinPreparation_PDBFixer()(filtered, out_path)
        return self._prepared_protein

    # ─────────────────────── ligand preparation (single) ───────────────

    @staticmethod
    def _add_hs(path: Path) -> Path:
        tmp = NamedTemporaryFile(suffix=".sdf", delete=False)
        dm.to_sdf([dm.add_hs(m, add_coords=True) for m in dm.read_sdf(path)], tmp.name)
        return Path(tmp.name)

    @staticmethod
    def _cdpk_fix(path: Path, mode: Literal["Dock", "Minimize"]) -> Path:
        tmp = NamedTemporaryFile(suffix=".sdf", delete=False)
        CDPK_Runner(gen3d=(mode == "Dock")).prepare_ligands(path, Path(tmp.name))
        return Path(tmp.name)

    def _prepare_single_ligand(self, mode: Literal["Dock", "Minimize"]) -> Tuple[Path, Path]:
        if self.input_mode == "SDF":
            fixed = self._cdpk_fix(self.ligands_sdf, mode)
        else:
            mol = dm.to_mol(self.ligands_smiles, sanitize=True, kekulize=True)
            mol.SetProp("_Name", "ligand")
            tmp = Path(gettempdir()) / "tmp_smiles.sdf"; dm.to_sdf(mol, tmp)
            fixed = self._cdpk_fix(tmp, mode)
        return self._add_hs(fixed), self.crystal_sdf.resolve()

    # ────────────────── explode to per‑ligand SDFs (VS) ─────────────────

    def _explode_multi_ligands(self, tmpdir: Path) -> List[Tuple[Path, str]]:
        """Return list of (sdf_path, ligand_name) for VS."""
        ligs: List[Tuple[Path, str]] = []
        if self.input_mode == "SDF":
            mols = dm.read_sdf(self.ligands_sdf, sanitize=False)
            assert len(mols) > 1, "virtual_screening=True but SDF has ≤1 record"
            for i, m in enumerate(mols, 1):
                name = m.GetProp("_Name") or f"lig_{i}"
                p = tmpdir / f"{name}.sdf"; dm.to_sdf(m, p); ligs.append((p, name))
        else:
            smi_list = [s for s in self.ligands_smiles.split(".") if s]
            assert len(smi_list) > 1, "virtual_screening=True but only one SMILES"
            for i, smi in enumerate(smi_list, 1):
                mol = dm.to_mol(smi, sanitize=True, kekulize=True)
                mol.SetProp("_Name", f"lig_{i}")
                p = tmpdir / f"lig_{i}.sdf"; dm.to_sdf(mol, p); ligs.append((p, f"lig_{i}"))
        return ligs

    # ────────────────────────── smina wrapper ──────────────────────────

    @staticmethod
    def _smina_single(
        protein: Path,
        crystal: Path,
        lig_tuple: Tuple[Path, str],
        exhaust: int,
        outdir: Path,
    ) -> Dict[str, object]:
        lig_path, lig_name = lig_tuple
        out_base = outdir / lig_name
        out_sdf, out_log = out_base.with_suffix(".sdf"), out_base.with_suffix(".log")

        cmd = [
            "smina", "-r", protein, "-l", lig_path,
            "--autobox_ligand", crystal,
            "--exhaustiveness", str(exhaust), "--cpu", "1",
            "-o", out_sdf,
        ]
        with open(out_log, "w") as lf:
            subprocess.run(cmd, check=True, stdout=lf, stderr=lf)

        # pull affinity
        aff = None
        try:
            mol = dm.read_sdf(out_sdf, sanitize=False)[0]
            if mol.HasProp(_SMINA_TAG):
                aff = float(mol.GetProp(_SMINA_TAG))
        except Exception:
            pass

        return {"name": lig_name, "sdf": out_sdf, "log": out_log, "aff": aff}

    # ─────────────────────────── public API ────────────────────────────

    def run_single(self, mode: Literal["Dock", "Minimize"], basename: str) -> Dict[str, object]:
        """Legacy single‑ligand path (kept for backward compatibility)."""
        fixed_lig, crystal = self._prepare_single_ligand(mode)
        protein = self.prepare_protein()

        result = self._smina_single(protein, crystal, (fixed_lig, basename), 32, self.workdir)

        payload = {
            "success": True,
            "docked_ligand": result["sdf"],
            "log": result["log"],
            "prepared_protein": protein,
        }
        if mode == "Minimize":
            payload["minimized_complex"] = minimize_complex(protein, dm.read_sdf(result["sdf"])[0])
        return payload

    def run_virtual_screen(self, basename: str) -> Dict[str, object]:
        """Execute docking‑only screen across **multiple** ligands."""
        assert self.virtual_screening, "Instantiate with virtual_screening=True first"

        protein = self.prepare_protein()
        crystal = self.crystal_sdf.resolve()

        tmpdir = Path(gettempdir()) / f"vs_{os.getpid()}"; tmpdir.mkdir(exist_ok=True)
        lig_list = self._explode_multi_ligands(tmpdir)

        cpu = max(1, int(os.cpu_count() * 0.9))
        worker = partial(self._smina_single, protein, crystal, exhaust=8, outdir=self.workdir)

        with mp.Pool(cpu) as pool:
            res = pool.map(worker, lig_list)

        # summary CSV
        csv_path = self.workdir / f"{basename}_summary.csv"
        with open(csv_path, "w", newline="") as fh:
            wr = csv.writer(fh); wr.writerow(["ligand", _SMINA_TAG])
            for r in res: wr.writerow([r["name"], r["aff"]])

        return {
            "success": True,
            "docked_ligands": [r["sdf"] for r in res],
            "logs": [r["log"] for r in res],
            "summary": csv_path,
        }

    # convenience alias --------------------------------------------------------
    def run(self, mode: Literal["Dock", "Minimize"], basename: str) -> Dict[str, object]:
        if self.virtual_screening:
            assert mode == "Dock", "Minimization not supported in VS mode"
            return self.run_virtual_screen(basename)
        return self.run_single(mode, basename)
