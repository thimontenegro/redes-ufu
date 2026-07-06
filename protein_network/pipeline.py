from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from .communities import (
    detect_kernighan_lin,
    detect_louvain,
    kernighan_lin_hexon_penton_counts,
    partition_stats,
)
from .graph_builder import build_residue_graph
from .io import build_residue_table, load_atoms, parse_mmcif_annotations
from .metrics import compute_centralities, degree_distribution, top_k
from .validation import validate_penton_hexon


@dataclass
class PipelineConfig:
    structure_path: str
    output_dir: str = "redes/output"
    cutoff: float = 8.0
    chain_mode: str = "any"
    weight_mode: str = "hybrid_affinity"
    min_affinity: int = 0
    connect_backbone: bool = True
    weighted_metrics: bool = False


def _safe_name(path: str) -> str:
    return Path(path).stem.replace(" ", "_")


def _plot_degree_distribution(distribution: dict[int, int], output_path: Path) -> None:
    fig = plt.figure()
    fig.set_size_inches(8.0, 4.0)
    xs = list(distribution.keys())
    ys = list(distribution.values())
    plt.bar(xs, ys, color="#4677c9")
    plt.title("Distribuicao de graus")
    plt.xlabel("Grau")
    plt.ylabel("Numero de nos")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def _plot_communities_2d(graph: nx.Graph, residues, partition: dict[int, int], output_path: Path) -> None:
    if graph.number_of_nodes() == 0:
        return

    pos = {
        int(row.node_id): np.array([float(row.x), float(row.y)])
        for row in residues.itertuples(index=False)
    }

    colors = [partition.get(node, -1) for node in graph.nodes()]
    sizes = [20 + 3 * graph.degree(node) for node in graph.nodes()]

    fig = plt.figure()
    fig.set_size_inches(8.0, 8.0)
    nx.draw_networkx_edges(graph, pos=pos, alpha=0.05, width=0.25)
    nodes = nx.draw_networkx_nodes(
        graph,
        pos=pos,
        node_color=colors,
        node_size=sizes,
        cmap=plt.cm.tab20,
    )
    plt.colorbar(nodes, label="Comunidade")
    plt.title("Comunidades sobre coordenadas reais (XY)")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def _top_report(graph: nx.Graph, centralities: dict, k: int = 10) -> dict:
    report = {}
    for metric_name, scores in centralities.items():
        report[metric_name] = []
        for node, score in top_k(scores, k=k):
            attrs = graph.nodes[node]
            report[metric_name].append(
                {
                    "node": int(node),
                    "score": float(score),
                    "residue": int(attrs.get("residue_id", -1)),
                    "chain": str(attrs.get("chain_id", "")),
                    "residue_name": str(attrs.get("residue_name", "")),
                }
            )
    return report


def run_pipeline(config: PipelineConfig) -> dict:
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    atoms = load_atoms(config.structure_path)
    residues = build_residue_table(atoms, atom_name="CA")

    graph = build_residue_graph(
        residues,
        cutoff=config.cutoff,
        weight_mode=config.weight_mode,
        chain_mode=config.chain_mode,
        min_affinity=config.min_affinity,
        connect_backbone=config.connect_backbone,
    )

    degree_dist = degree_distribution(graph)
    centralities = compute_centralities(graph, weighted=config.weighted_metrics)

    louvain_partition = detect_louvain(graph)
    kl_partition = detect_kernighan_lin(graph)

    chain_to_entity, entity_to_desc = parse_mmcif_annotations(config.structure_path)
    louvain_validation = validate_penton_hexon(
        residues, louvain_partition, chain_to_entity, entity_to_desc
    )
    kl_counts = kernighan_lin_hexon_penton_counts(
        graph,
        kl_partition,
        chain_to_entity,
        entity_to_desc,
    )

    prefix = _safe_name(config.structure_path)
    _plot_degree_distribution(degree_dist, output_dir / f"{prefix}_degree_distribution.png")
    _plot_communities_2d(graph, residues, louvain_partition, output_dir / f"{prefix}_louvain_xy.png")

    result = {
        "config": asdict(config),
        "n_atoms": int(len(atoms)),
        "n_residues": int(len(residues)),
        "n_nodes": int(graph.number_of_nodes()),
        "n_edges": int(graph.number_of_edges()),
        "density": float(nx.density(graph)) if graph.number_of_nodes() > 1 else 0.0,
        "degree_distribution": {str(k): int(v) for k, v in degree_dist.items()},
        "louvain": {
            "stats": partition_stats(louvain_partition),
            "validation": louvain_validation,
        },
        "kernighan_lin": {
            "stats": partition_stats(kl_partition),
            "community_label_counts": kl_counts,
        },
        "top_central_nodes": _top_report(graph, centralities, k=10),
        "entity_descriptions": entity_to_desc,
    }

    result_path = output_dir / f"{prefix}_summary.json"
    result_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


