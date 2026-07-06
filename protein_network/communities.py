from __future__ import annotations

from collections import defaultdict

import networkx as nx
from community import community_louvain
from networkx.algorithms import community as nx_community


def detect_louvain(graph: nx.Graph) -> dict[int, int]:
    if graph.number_of_nodes() == 0:
        return {}
    return community_louvain.best_partition(graph, weight="weight")


def detect_kernighan_lin(graph: nx.Graph) -> dict[int, int]:
    if graph.number_of_nodes() == 0:
        return {}

    comp = graph
    if not nx.is_connected(graph):
        largest = max(nx.connected_components(graph), key=len)
        comp = graph.subgraph(largest).copy()

    left, right = nx_community.kernighan_lin_bisection(comp, weight="weight")
    partition = {node: 0 for node in left}
    partition.update({node: 1 for node in right})
    return partition


def kernighan_lin_hexon_penton_counts(
    graph: nx.Graph,
    partition: dict[int, int],
    chain_to_entity: dict[str, str],
    entity_to_description: dict[str, str],
) -> dict[str, dict[str, int]]:
    """Count hexon/penton/other nodes per community for a given partition."""
    if not partition:
        return {}

    counts: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    for node, community in partition.items():
        attrs = graph.nodes.get(node, {})
        chain = str(attrs.get("chain_id", ""))
        entity_id = chain_to_entity.get(chain)
        description = entity_to_description.get(entity_id, "").lower()

        if "hexon" in description:
            label = "hexon"
        elif "penton" in description:
            label = "penton"
        else:
            label = "other"

        counts[int(community)][label] += 1

    return {
        str(comm): {
            "hexon": labels.get("hexon", 0),
            "penton": labels.get("penton", 0),
            "other": labels.get("other", 0),
        }
        for comm, labels in sorted(counts.items())
    }


def partition_stats(partition: dict[int, int]) -> dict:
    if not partition:
        return {"n_communities": 0, "sizes": {}}

    sizes: dict[int, int] = {}
    for comm in partition.values():
        sizes[comm] = sizes.get(comm, 0) + 1

    return {
        "n_communities": len(sizes),
        "sizes": dict(sorted(sizes.items(), key=lambda kv: kv[1], reverse=True)),
    }

