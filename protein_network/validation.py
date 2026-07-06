from __future__ import annotations

from collections import defaultdict

import pandas as pd


EMPTY_RESULT = {
    "evaluated_nodes": 0,
    "labeled_nodes": 0,
    "purity": None,
    "community_label_counts": {},
}


def _binary_label(description: str) -> str:
    desc = description.lower()
    if "hexon" in desc:
        return "hexon"
    if "penton" in desc:
        return "penton"
    return "other"


def _compute_purity(community_counts: dict[int, dict[str, int]]) -> float | None:
    if not community_counts:
        return None

    total = 0
    dominant = 0
    for labels in community_counts.values():
        community_total = sum(labels.values())
        total += community_total
        dominant += max(labels.values(), default=0)

    if total == 0:
        return None
    return dominant / total


def validate_penton_hexon(
    residues: pd.DataFrame,
    partition: dict[int, int],
    chain_to_entity: dict[str, str],
    entity_to_description: dict[str, str],
) -> dict:
    """Evaluate how well communities separate penton and hexon annotations."""
    if not partition or not chain_to_entity or not entity_to_description:
        return EMPTY_RESULT.copy()

    community_counts: dict[int, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    labeled_nodes = 0

    for _, row in residues.iterrows():
        node_id = int(row["node_id"])
        if node_id not in partition:
            continue

        chain = str(row["chain_id"])
        entity = chain_to_entity.get(chain)
        description = entity_to_description.get(entity, "")
        label = _binary_label(description)

        community = int(partition[node_id])
        community_counts[community][label] += 1
        if label in {"hexon", "penton"}:
            labeled_nodes += 1

    purity = _compute_purity(community_counts) if labeled_nodes > 0 else None

    evaluated_nodes = sum(sum(v.values()) for v in community_counts.values())

    serializable_counts = {
        str(comm): dict(sorted(labels.items())) for comm, labels in sorted(community_counts.items())
    }

    return {
        "evaluated_nodes": int(evaluated_nodes),
        "labeled_nodes": labeled_nodes,
        "purity": purity,
        "community_label_counts": serializable_counts,
    }

