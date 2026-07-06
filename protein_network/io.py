from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from biopandas.mmcif import PandasMmcif
from biopandas.pdb import PandasPdb


def _to_int(value: object) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def load_atoms(structure_path: str) -> pd.DataFrame:
    """Load a PDB/mmCIF structure and normalize atom columns."""
    path = Path(structure_path)
    ext = path.suffix.lower()

    if ext in {".cif", ".mmcif"}:
        mmcif = PandasMmcif().read_mmcif(str(path))
        df = mmcif.df["ATOM"].copy()
        chain_col = "auth_asym_id" if "auth_asym_id" in df.columns else "label_asym_id"
        normalized = pd.DataFrame(
            {
                "atom_name": df["label_atom_id"],
                "residue_name": df["label_comp_id"],
                "residue_id": df["auth_seq_id"].map(_to_int),
                "chain_id": df[chain_col],
                "x": pd.to_numeric(df["Cartn_x"], errors="coerce"),
                "y": pd.to_numeric(df["Cartn_y"], errors="coerce"),
                "z": pd.to_numeric(df["Cartn_z"], errors="coerce"),
            }
        )
        return normalized

    if ext == ".pdb":
        pdb = PandasPdb().read_pdb(str(path))
        df = pdb.df["ATOM"].copy()
        normalized = pd.DataFrame(
            {
                "atom_name": df["atom_name"],
                "residue_name": df["residue_name"],
                "residue_id": df["residue_number"].map(_to_int),
                "chain_id": df["chain_id"],
                "x": pd.to_numeric(df["x_coord"], errors="coerce"),
                "y": pd.to_numeric(df["y_coord"], errors="coerce"),
                "z": pd.to_numeric(df["z_coord"], errors="coerce"),
            }
        )
        return normalized

    raise ValueError(f"Unsupported structure format: {path}")


def build_residue_table(atoms_df: pd.DataFrame, atom_name: str = "CA") -> pd.DataFrame:
    """Build residue-level table from atoms using one representative atom per residue."""
    residues = atoms_df.loc[atoms_df["atom_name"] == atom_name].copy()
    residues = residues.dropna(subset=["residue_id", "chain_id", "x", "y", "z"])
    residues = residues.sort_values(["chain_id", "residue_id"]).drop_duplicates(
        subset=["chain_id", "residue_id"], keep="first"
    )
    residues = residues.reset_index(drop=True)
    residues["node_id"] = residues.index
    return residues


def parse_mmcif_annotations(structure_path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Extract chain->entity and entity->description maps from mmCIF annotations."""
    path = Path(structure_path)
    if path.suffix.lower() not in {".cif", ".mmcif"}:
        return {}, {}

    raw = MMCIF2Dict(str(path))

    entity_ids = raw.get("_entity.id", [])
    descriptions = raw.get("_entity.pdbx_description", [])
    chains = raw.get("_struct_asym.id", [])
    chain_entities = raw.get("_struct_asym.entity_id", [])

    if isinstance(entity_ids, str):
        entity_ids = [entity_ids]
    if isinstance(descriptions, str):
        descriptions = [descriptions]
    if isinstance(chains, str):
        chains = [chains]
    if isinstance(chain_entities, str):
        chain_entities = [chain_entities]

    entity_to_description = {
        str(entity_id): str(description)
        for entity_id, description in zip(entity_ids, descriptions, strict=False)
    }
    chain_to_entity = {
        str(chain): str(entity_id)
        for chain, entity_id in zip(chains, chain_entities, strict=False)
    }
    return chain_to_entity, entity_to_description

