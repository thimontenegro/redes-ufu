from __future__ import annotations

from collections import Counter
from typing import Dict

import networkx as nx


def degree_distribution(graph: nx.Graph) -> Dict[int, int]:
    counts = Counter(graph.degree(node) for node in graph.nodes())
    return dict(sorted(counts.items()))


def compute_centralities(graph: nx.Graph, weighted: bool = False) -> dict:
    if graph.number_of_nodes() == 0:
        return {
            "degree": {},
            "betweenness": {},
            "closeness": {},
            "eigenvector": {},
        }

    weight_key = "weight" if weighted else None
    degree = nx.degree_centrality(graph)

    n_nodes = graph.number_of_nodes()
    if n_nodes > 3000:
        # Aproximacao para grafos muito grandes (ex.: 6B1T) para manter tempo viavel.
        k = min(500, n_nodes)
        betweenness = nx.betweenness_centrality(graph, k=k, weight=weight_key, seed=42)

        ranked_nodes = sorted(graph.degree, key=lambda item: item[1], reverse=True)
        sampled_nodes = [node for node, _ in ranked_nodes[: min(1000, n_nodes)]]
        closeness = {node: nx.closeness_centrality(graph, u=node) for node in sampled_nodes}
    else:
        betweenness = nx.betweenness_centrality(graph, weight=weight_key)
        closeness = nx.closeness_centrality(graph)

    try:
        # Tenta convergência via power-iteration com mais iterações
        eigenvector = nx.eigenvector_centrality(graph, weight=weight_key, max_iter=1000, tol=1e-5)
    except nx.PowerIterationFailedConvergence:
        # Para grafos grandes (ex.: 6B1T) a power-iteration pode não convergir;
        # usa solver numpy/ARPACK que é sempre estável.
        try:
            eigenvector = nx.eigenvector_centrality_numpy(graph, weight=weight_key)
        except nx.NetworkXException:
            eigenvector = {}
    except nx.NetworkXException:
        eigenvector = {}

    return {
        "degree": degree,
        "betweenness": betweenness,
        "closeness": closeness,
        "eigenvector": eigenvector,
    }


def top_k(scores: dict, k: int = 10) -> list[tuple[int, float]]:
    return sorted(scores.items(), key=lambda item: item[1], reverse=True)[:k]

