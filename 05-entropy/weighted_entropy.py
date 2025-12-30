#!/usr/bin/env python3

"""
weighted_entropy.py

Compute per-read taxonomy entropy H(t_true, t_pred) using NCBI taxonomy.

Entropy definition (for t_true != t_pred):

    H(t_true, t_pred) = R(rank(L)) * log2(1 + k_L) * (alpha_up * u + alpha_down * d)

where:
    - L        = lowest common ancestor (LCA) of true and predicted taxids
    - u        = edges UP from true taxid to L
    - d        = edges DOWN from L to predicted taxid
    - k_L      = number of children of L
    - R(rank)  = rank weight (higher ranks get larger weight)
"""

import argparse
import math
from collections import defaultdict
from typing import Dict, Tuple, Optional

import pandas as pd


# ----------------------------
# Default hyperparameters
# ----------------------------

DEFAULT_ALPHA_UP = 0.3
DEFAULT_ALPHA_DOWN = 1.0

DEFAULT_RANK_WEIGHTS: Dict[str, float] = {
    # below- / at-species: small weights
    "subspecies": 0.15,
    "forma": 0.15,
    "varietas": 0.15,
    "species": 0.20,

    # typical core ranks
    "genus": 0.30,
    "family": 0.40,
    "order": 0.50,
    "class": 0.60,
    "phylum": 0.80,
    "division": 0.80,        # plants/fungi use "division"
    "kingdom": 1.00,
    "superkingdom": 1.00,
}

DEFAULT_FALLBACK_RANK_WEIGHT = 0.50  # for unranked/odd ranks


# ----------------------------
# Taxonomy loading & helpers
# ----------------------------

def load_taxonomy(nodes_path: str) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, int], Dict[str, int]]:
    """
    Load NCBI taxonomy from nodes.dmp.

    Returns:
        parent: taxid -> parent taxid
        rank:   taxid -> rank string
        depth:  taxid -> depth (edges from root)
        branch_size: taxid -> number of children
    """
    parent: Dict[str, str] = {}
    rank: Dict[str, str] = {}
    children: Dict[str, list] = defaultdict(list)

    with open(nodes_path, "r") as f:
        for line in f:
            # nodes.dmp is pipe-separated, like: tax_id \t|\t parent_tax_id \t|\t rank \t| ...
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 3:
                continue
            tax_id = parts[0]
            parent_id = parts[1]
            rank_name = parts[2]

            parent[tax_id] = parent_id
            rank[tax_id] = rank_name
            children[parent_id].append(tax_id)

    # Compute depths from root via DFS with memoization
    depth: Dict[str, int] = {}

    def get_depth(t: str) -> int:
        if t in depth:
            return depth[t]
        p = parent.get(t)
        if p is None or p == t:
            depth[t] = 0
        else:
            # guard against missing parent loops
            if p not in parent:
                depth[t] = 0
            else:
                depth[t] = get_depth(p) + 1
        return depth[t]

    for t in parent.keys():
        get_depth(t)

    branch_size: Dict[str, int] = {t: len(children.get(t, [])) for t in parent.keys()}

    return parent, rank, depth, branch_size


def lineage_to_root(t: str, parent: Dict[str, str]) -> list:
    """Return [t, parent(t), ..., root]."""
    path = []
    seen = set()
    while t is not None and t not in seen:
        path.append(t)
        seen.add(t)
        p = parent.get(t)
        if p is None or p == t:
            break
        t = p
    return path


def lca_features(t1: str,
                 t2: str,
                 parent: Dict[str, str],
                 depth: Dict[str, int],
                 branch_size: Dict[str, int],
                 rank: Dict[str, str]) -> Tuple[str, int, int, int, str]:
    """
    Compute LCA and geometric features for (t1, t2).

    Returns:
        L:          LCA taxid
        up:         edges up from t1 to L
        down:       edges down from L to t2
        k_children: number of children of L
        rank_L:     rank name of L
    """
    if t1 == t2:
        L = t1
        up = 0
        down = 0
    else:
        path1 = lineage_to_root(t1, parent)
        path2 = lineage_to_root(t2, parent)
        set1 = set(path1)
        L = None
        for node in path2:
            if node in set1:
                L = node
                break
        if L is None:
            # In a well-formed NCBI taxonomy this shouldn't happen.
            # If it does, fall back to treating t1 as LCA.
            L = t1

        up = max(depth.get(t1, 0) - depth.get(L, 0), 0)
        down = max(depth.get(t2, 0) - depth.get(L, 0), 0)

    k = branch_size.get(L, 0)
    rank_L = rank.get(L, "no_rank")

    return L, up, down, k, rank_L


# ----------------------------
# Entropy components
# ----------------------------

def rank_weight(rank_name: str,
                rank_weights: Dict[str, float],
                fallback: float = DEFAULT_FALLBACK_RANK_WEIGHT) -> float:
    """Map a rank string to a scalar weight."""
    return rank_weights.get(rank_name, fallback)


def local_branch_entropy(k_children: int) -> float:
    """
    Structural 'branch entropy' from number of children.
    Use log2(1 + k) to grow slower than k and be 0 when k=0.
    """
    if k_children <= 0:
        return 0.0
    return math.log2(1.0 + k_children)


def entropy_for_pair(true_taxid: str,
                     pred_taxid: str,
                     parent: Dict[str, str],
                     rank: Dict[str, str],
                     depth: Dict[str, int],
                     branch_size: Dict[str, int],
                     alpha_up: float = DEFAULT_ALPHA_UP,
                     alpha_down: float = DEFAULT_ALPHA_DOWN,
                     rank_weights: Optional[Dict[str, float]] = None,
                     fallback_rank_weight: float = DEFAULT_FALLBACK_RANK_WEIGHT,
                     unclassified_sentinels=None,
                     unclassified_entropy: Optional[float] = None) -> Optional[float]:
    """
    Compute entropy H(true_taxid, pred_taxid) for a single read.

    Returns:
        H (float) or None if cannot be computed.
    """
    if rank_weights is None:
        rank_weights = DEFAULT_RANK_WEIGHTS

    if unclassified_sentinels is None:
        unclassified_sentinels = {"0", "", None}

    # Normalize input
    if isinstance(true_taxid, float) and math.isnan(true_taxid):
        return None
    if isinstance(pred_taxid, float) and math.isnan(pred_taxid):
        return None

    t1 = str(true_taxid)
    t2 = str(pred_taxid)

    # Unclassified / missing prediction: assign fixed entropy (or None)
    if t2 in unclassified_sentinels:
        return unclassified_entropy

    # Require both taxids to exist in taxonomy
    if t1 not in parent or t2 not in parent:
        return None

    # Exact match
    if t1 == t2:
        return 0.0

    L, u, d, k, rank_L = lca_features(t1, t2, parent, depth, branch_size, rank)

    branch_H = local_branch_entropy(k)
    if branch_H == 0.0:
        # If no children, use path length alone (treat branch_H as 1.0)
        branch_H = 1.0

    R = rank_weight(rank_L, rank_weights, fallback_rank_weight)
    path_penalty = alpha_up * u + alpha_down * d

    H = R * branch_H * path_penalty
    return H


# ----------------------------
# I/O and main driver
# ----------------------------

READID_CANDIDATES = ["read_id", "readID", "SequenceID", "seq_id", "read"]
TRUE_TAXID_CANDIDATES = ["true_taxid", "truth_taxid", "true", "gold_taxid"]
PRED_TAXID_CANDIDATES = ["pred_taxid", "predicted_taxid", "call_taxid", "pred"]


def pick_column(df: pd.DataFrame, candidates) -> str:
    for name in candidates:
        if name in df.columns:
            return name
    raise ValueError(f"None of the expected columns {candidates} found in input file. "
                     f"Available columns: {list(df.columns)}")


def main():
    parser = argparse.ArgumentParser(description="Compute taxonomy entropy for true/predicted taxid pairs.")
    parser.add_argument("--nodes-dmp", required=True,
                        help="Path to NCBI nodes.dmp file.")
    parser.add_argument("--input-tsv", required=True,
                        help="Tab-separated file with read IDs, true_taxid, pred_taxid.")
    parser.add_argument("--output-tsv", required=True,
                        help="Where to write per-read entropy results (TSV).")
    parser.add_argument("--alpha-up", type=float, default=DEFAULT_ALPHA_UP,
                        help=f"Penalty per 'up' edge (default {DEFAULT_ALPHA_UP}).")
    parser.add_argument("--alpha-down", type=float, default=DEFAULT_ALPHA_DOWN,
                        help=f"Penalty per 'down' edge (default {DEFAULT_ALPHA_DOWN}).")
    parser.add_argument("--unclassified-entropy", type=float, default=None,
                        help="Entropy value to assign when prediction is unclassified (taxid 0/blank). "
                             "Default: None (leave entropy as NA).")

    args = parser.parse_args()

    print(f"Loading taxonomy from: {args.nodes_dmp}")
    parent, rank, depth, branch_size = load_taxonomy(args.nodes_dmp)
    print(f"Loaded {len(parent)} taxonomy nodes.")

    print(f"Reading input TSV: {args.input_tsv}")
    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str)

    read_col = pick_column(df, READID_CANDIDATES)
    true_col = pick_column(df, TRUE_TAXID_CANDIDATES)
    pred_col = pick_column(df, PRED_TAXID_CANDIDATES)

    print(f"Using columns: read_id={read_col}, true_taxid={true_col}, pred_taxid={pred_col}")

    entropies = []
    lca_taxids = []
    ups = []
    downs = []
    lca_ranks = []
    branch_sizes = []

    for idx, row in df.iterrows():
        t1 = row[true_col]
        t2 = row[pred_col]

        # Compute entropy
        H = entropy_for_pair(
            true_taxid=t1,
            pred_taxid=t2,
            parent=parent,
            rank=rank,
            depth=depth,
            branch_size=branch_size,
            alpha_up=args.alpha_up,
            alpha_down=args.alpha_down,
            unclassified_entropy=args.unclassified_entropy,
        )
        entropies.append(H)

        # Compute LCA-related diagnostics only when both taxids are in taxonomy
        L = None
        u = None
        d = None
        k = None
        rank_L = None

        if isinstance(t1, str) and isinstance(t2, str) and (t1 in parent) and (t2 in parent):
            if t1 == t2:
                L = t1
                u = 0
                d = 0
                k = branch_size.get(t1, 0)
                rank_L = rank.get(t1, "no_rank")
            else:
                L, u, d, k, rank_L = lca_features(t1, t2, parent, depth, branch_size, rank)

        lca_taxids.append(L)
        ups.append(u)
        downs.append(d)
        lca_ranks.append(rank_L)
        branch_sizes.append(k)

    df_out = df.copy()
    df_out["entropy"] = entropies
    df_out["lca_taxid"] = lca_taxids
    df_out["lca_rank"] = lca_ranks
    df_out["up_from_true"] = ups
    df_out["down_to_pred"] = downs
    df_out["branch_size_LCA"] = branch_sizes

    print(f"Writing per-read entropy to: {args.output_tsv}")
    df_out.to_csv(args.output_tsv, sep="\t", index=False)

    # Simple summary
    try:
        valid = df_out["entropy"].dropna().astype(float)
    except Exception:
        valid = pd.Series([], dtype=float)

    if len(valid) > 0:
        print("Entropy summary (valid pairs only):")
        print(f"  n = {len(valid)}")
        print(f"  mean = {valid.mean():.4f}")
        print(f"  median = {valid.median():.4f}")
        print(f"  min = {valid.min():.4f}")
        print(f"  max = {valid.max():.4f}")
    else:
        print("No valid entropy values computed (check taxids and input columns).")


if __name__ == "__main__":
    main()
