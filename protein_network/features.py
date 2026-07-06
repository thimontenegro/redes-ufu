from __future__ import annotations

RESIDUE_PROPERTIES = {
    "ALA": {"hydrophobic": True, "charge": 0, "polar": False, "donor": False, "acceptor": False},
    "VAL": {"hydrophobic": True, "charge": 0, "polar": False, "donor": False, "acceptor": False},
    "LEU": {"hydrophobic": True, "charge": 0, "polar": False, "donor": False, "acceptor": False},
    "ILE": {"hydrophobic": True, "charge": 0, "polar": False, "donor": False, "acceptor": False},
    "MET": {"hydrophobic": True, "charge": 0, "polar": False, "donor": False, "acceptor": True},
    "PHE": {"hydrophobic": True, "charge": 0, "polar": False, "donor": False, "acceptor": False},
    "TRP": {"hydrophobic": True, "charge": 0, "polar": False, "donor": True, "acceptor": False},
    "PRO": {"hydrophobic": True, "charge": 0, "polar": False, "donor": False, "acceptor": False},
    "SER": {"hydrophobic": False, "charge": 0, "polar": True, "donor": True, "acceptor": True},
    "THR": {"hydrophobic": False, "charge": 0, "polar": True, "donor": True, "acceptor": True},
    "ASN": {"hydrophobic": False, "charge": 0, "polar": True, "donor": True, "acceptor": True},
    "GLN": {"hydrophobic": False, "charge": 0, "polar": True, "donor": True, "acceptor": True},
    "TYR": {"hydrophobic": False, "charge": 0, "polar": True, "donor": True, "acceptor": True},
    "CYS": {"hydrophobic": False, "charge": 0, "polar": True, "donor": True, "acceptor": True},
    "GLY": {"hydrophobic": False, "charge": 0, "polar": False, "donor": False, "acceptor": False},
    "LYS": {"hydrophobic": False, "charge": 1, "polar": True, "donor": True, "acceptor": False},
    "ARG": {"hydrophobic": False, "charge": 1, "polar": True, "donor": True, "acceptor": False},
    "HIS": {"hydrophobic": False, "charge": 1, "polar": True, "donor": True, "acceptor": True},
    "ASP": {"hydrophobic": False, "charge": -1, "polar": True, "donor": False, "acceptor": True},
    "GLU": {"hydrophobic": False, "charge": -1, "polar": True, "donor": False, "acceptor": True},
}


def affinity_score(res1: str, res2: str) -> int:
    """Simple physicochemical affinity score between two residue types."""
    p1 = RESIDUE_PROPERTIES.get(res1, RESIDUE_PROPERTIES["GLY"])
    p2 = RESIDUE_PROPERTIES.get(res2, RESIDUE_PROPERTIES["GLY"])

    score = 0
    if p1["hydrophobic"] and p2["hydrophobic"]:
        score += 2
    if p1["charge"] * p2["charge"] == -1:
        score += 3
    if p1["polar"] and p2["polar"]:
        score += 1
    if (p1["donor"] and p2["acceptor"]) or (p2["donor"] and p1["acceptor"]):
        score += 1
    return score

