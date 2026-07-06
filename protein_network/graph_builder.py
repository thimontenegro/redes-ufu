from __future__ import annotations

import numpy as np
import pandas as pd
import networkx as nx
from sklearn.neighbors import NearestNeighbors

from .features import affinity_score


def pairwise_distances(coords: np.ndarray) -> np.ndarray:
    diff = coords[:, None, :] - coords[None, :, :]
    return np.sqrt((diff * diff).sum(axis=2))


def build_residue_graph(
    residues: pd.DataFrame,
    cutoff: float = 8.0,
    weight_mode: str = "inverse_distance",
    chain_mode: str = "any",
    min_affinity: int = 0,
    connect_backbone: bool = True,
) -> nx.Graph:
    """Build an undirected residue graph from CA coordinates."""
    required_columns = {"node_id", "residue_name", "residue_id", "chain_id", "x", "y", "z"}
    missing = required_columns.difference(residues.columns)
    if missing:
        raise ValueError(f"Residue table missing columns: {sorted(missing)}")

    graph = nx.Graph()
    for _, row in residues.iterrows():
        graph.add_node(
            int(row["node_id"]),
            residue_name=str(row["residue_name"]),
            residue_id=int(row["residue_id"]),
            chain_id=str(row["chain_id"]),
        )

    coords = residues[["x", "y", "z"]].to_numpy(dtype=float)
    neighbors = NearestNeighbors(radius=cutoff, algorithm="kd_tree")
    neighbors.fit(coords)
    radius_dists, radius_idx = neighbors.radius_neighbors(coords, return_distance=True)

    seen_pairs: set[tuple[int, int]] = set()
    for i in range(len(residues)):
        row_i = residues.iloc[i]
        for dist, j in zip(radius_dists[i], radius_idx[i], strict=False):
            if j <= i:
                continue
            if (i, int(j)) in seen_pairs:
                continue
            seen_pairs.add((i, int(j)))

            row_j = residues.iloc[j]
            same_chain = str(row_i["chain_id"]) == str(row_j["chain_id"])

            if chain_mode == "same" and not same_chain:
                continue
            if chain_mode == "different" and same_chain:
                continue

            dist = float(dist)
            affinity = affinity_score(str(row_i["residue_name"]), str(row_j["residue_name"]))
            is_backbone_neighbor = same_chain and abs(int(row_i["residue_id"]) - int(row_j["residue_id"])) == 1

            if affinity < min_affinity:
                continue

            if weight_mode == "unit":
                weight = 1.0
            elif weight_mode == "hybrid_affinity":
                weight = (1.0 + affinity) / (dist + 1e-6)
            else:
                weight = 1.0 / (dist + 1e-6)

            graph.add_edge(
                int(row_i["node_id"]),
                int(row_j["node_id"]),
                distance=dist,
                affinity=affinity,
                weight=float(weight),
                same_chain=same_chain,
                backbone_neighbor=bool(is_backbone_neighbor),
            )

    if connect_backbone:
        rows = residues.sort_values(["chain_id", "residue_id"])
        prev_by_chain: dict[str, pd.Series] = {}
        for row in rows.itertuples(index=False):
            chain = str(row.chain_id)
            prev = prev_by_chain.get(chain)
            if prev is not None and int(row.residue_id) - int(prev.residue_id) == 1:
                u = int(prev.node_id)
                v = int(row.node_id)
                if not graph.has_edge(u, v):
                    coord_u = np.array([prev.x, prev.y, prev.z], dtype=float)
                    coord_v = np.array([row.x, row.y, row.z], dtype=float)
                    dist = float(np.linalg.norm(coord_u - coord_v))
                    graph.add_edge(
                        u,
                        v,
                        distance=dist,
                        affinity=affinity_score(str(prev.residue_name), str(row.residue_name)),
                        weight=1.0 / (dist + 1e-6),
                        same_chain=True,
                        backbone_neighbor=True,
                    )
            prev_by_chain[chain] = row

    return graph

