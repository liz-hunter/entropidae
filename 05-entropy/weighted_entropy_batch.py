#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import multiprocessing as mp
import os
import re
from pathlib import Path
from typing import Dict, Tuple, Optional, List

# Import your existing code (file must be named weighted_entropy.py)
from weighted_entropy import (
    load_taxonomy,
    lineage_to_root,
    local_branch_entropy,
    rank_weight,
    DEFAULT_RANK_WEIGHTS,
    DEFAULT_FALLBACK_RANK_WEIGHT,
)

# ----------------------------
# Globals for worker processes
# ----------------------------
PARENT = None
RANK = None
DEPTH = None
BRANCH_SIZE = None

FNAME_RE = re.compile(r"^(?P<dataset>[^_]+)_(?P<filename>.+)_db(?P<db>\d+)\.out(?:\.gz)?$")
TAXID_RE = re.compile(r"taxid\s+(\d+)")


def open_text(path: str, mode: str, compresslevel: int = 3):
    if path.endswith(".gz"):
        return gzip.open(path, mode, compresslevel=compresslevel)
    return open(path, mode)


def init_worker(nodes_dmp: str):
    """Load taxonomy once per worker (fork shares memory on Linux; spawn loads per worker)."""
    global PARENT, RANK, DEPTH, BRANCH_SIZE
    if PARENT is None:
        PARENT, RANK, DEPTH, BRANCH_SIZE = load_taxonomy(nodes_dmp)


def parse_pred_taxid(field3: str, status: str) -> str:
    """Extract predicted taxid from Kraken field 3 (which may be numeric or 'Name (taxid N)')."""
    if status == "U":
        return "0"
    s = (field3 or "").strip()
    if not s:
        return "0"
    if s.isdigit():
        return s
    m = TAXID_RE.search(s)
    if m:
        return m.group(1)
    # If you truly used --use-names without taxid in the field, you won't be able to recover taxid here.
    return "0"


def compute_cached_for_pred(
    true_taxid: str,
    pred_taxid: str,
    cache: Dict[str, Tuple[Optional[float], str, str, Optional[int], Optional[int], Optional[int]]],
    true_lineage_set: set,
    true_depth: int,
    alpha_up: float,
    alpha_down: float,
    unclassified_entropy: Optional[float],
    unclassified_sentinels: set,
):
    """
    Cache results by pred_taxid (true_taxid constant for the file).

    Returns:
      (entropy, lca_taxid, lca_rank, up_from_true, down_to_pred, branch_size_LCA)
    """
    global PARENT, RANK, DEPTH, BRANCH_SIZE

    if pred_taxid in cache:
        return cache[pred_taxid]

    if pred_taxid in unclassified_sentinels:
        res = (unclassified_entropy, "", "", None, None, None)
        cache[pred_taxid] = res
        return res

    if (true_taxid not in PARENT) or (pred_taxid not in PARENT):
        res = (None, "", "", None, None, None)
        cache[pred_taxid] = res
        return res

    if true_taxid == pred_taxid:
        k = BRANCH_SIZE.get(true_taxid, 0)
        rank_L = RANK.get(true_taxid, "no_rank")
        res = (0.0, true_taxid, rank_L, 0, 0, k)
        cache[pred_taxid] = res
        return res

    # LCA: walk pred lineage until it hits true lineage set
    L = ""
    for node in lineage_to_root(pred_taxid, PARENT):
        if node in true_lineage_set:
            L = node
            break
    if not L:
        L = true_taxid

    depth_L = DEPTH.get(L, 0)
    depth_pred = DEPTH.get(pred_taxid, 0)

    u = max(true_depth - depth_L, 0)
    d = max(depth_pred - depth_L, 0)

    k = BRANCH_SIZE.get(L, 0)
    rank_L = RANK.get(L, "no_rank")

    branch_H = local_branch_entropy(k)
    if branch_H == 0.0:
        branch_H = 1.0

    Rw = rank_weight(rank_L, DEFAULT_RANK_WEIGHTS, DEFAULT_FALLBACK_RANK_WEIGHT)
    path_penalty = alpha_up * u + alpha_down * d
    H = Rw * branch_H * path_penalty

    res = (H, L, rank_L, u, d, k)
    cache[pred_taxid] = res
    return res


def load_key(key_path: str) -> Dict[Tuple[str, str], str]:
    """
    Read taxid2filename.txt (TSV with header).
    Expected columns: filename, basename, dataset, taxid
    """
    mapping: Dict[Tuple[str, str], str] = {}
    with open(key_path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name.strip(): i for i, name in enumerate(header)}

        for required in ("filename", "dataset", "taxid"):
            if required not in col:
                raise ValueError(f"Key file missing required column '{required}'. Found: {header}")

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            filename = parts[col["filename"]].strip()
            dataset = parts[col["dataset"]].strip()
            taxid = parts[col["taxid"]].strip()
            if filename and dataset and taxid:
                mapping[(dataset, filename)] = taxid
    return mapping


def discover_jobs(kraken_dir: str, key_map: Dict[Tuple[str, str], str], recursive: bool):
    """
    Scan kraken_dir for *_dbN.out(.gz), parse dataset+filename, look up true taxid in key_map.
    Returns list of dict jobs and list of missing_key filenames.
    """
    p = Path(kraken_dir)
    it = p.rglob("*") if recursive else p.iterdir()

    jobs = []
    missing_key = []

    for fp in it:
        if not fp.is_file():
            continue
        name = fp.name
        if not (name.endswith(".out") or name.endswith(".out.gz")):
            continue

        m = FNAME_RE.match(name)
        if not m:
            continue

        dataset = m.group("dataset")
        filename = m.group("filename")
        db = m.group("db")

        taxid = key_map.get((dataset, filename))
        if taxid is None:
            missing_key.append(str(fp))
            continue

        jobs.append(
            {
                "path": str(fp),
                "dataset": dataset,
                "filename": filename,
                "db": db,
                "true_taxid": taxid,
            }
        )

    return jobs, missing_key


def read_jobs_tsv(path: str) -> List[Dict[str, str]]:
    """
    Read a 2-column TSV: kraken_out_path<TAB>true_taxid
    Returns list of job dicts compatible with process_one().
    """
    jobs: List[Dict[str, str]] = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue

            kraken_path = parts[0].strip()
            true_taxid = parts[1].strip()

            name = Path(kraken_path).name
            m = FNAME_RE.match(name)
            if m:
                jobs.append(
                    {
                        "path": kraken_path,
                        "dataset": m.group("dataset"),
                        "filename": m.group("filename"),
                        "db": m.group("db"),
                        "true_taxid": true_taxid,
                    }
                )
            else:
                jobs.append(
                    {
                        "path": kraken_path,
                        "dataset": "",
                        "filename": name,
                        "db": "",
                        "true_taxid": true_taxid,
                    }
                )
    return jobs


def process_one(job):
    """
    job contains: path, dataset, filename, db, true_taxid, out_path, params...
    """
    global PARENT, RANK, DEPTH, BRANCH_SIZE

    path = job["path"]
    true_taxid = str(job["true_taxid"]).strip()
    out_path = job["out_path"]

    alpha_up = job["alpha_up"]
    alpha_down = job["alpha_down"]
    unclassified_entropy = job["unclassified_entropy"]
    write_diag = job["write_diag"]
    compresslevel = job["compresslevel"]
    skip_existing = job["skip_existing"]

    if skip_existing and os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return {**job, "status": "skipped_existing", "n_reads": 0, "n_valid": 0, "mean_entropy": ""}

    # Precompute true lineage
    if true_taxid in PARENT:
        true_lineage_set = set(lineage_to_root(true_taxid, PARENT))
        true_depth = DEPTH.get(true_taxid, 0)
    else:
        true_lineage_set = set()
        true_depth = 0

    unclassified_sentinels = {"0", "", "NA", "None", None}

    cache: Dict[str, Tuple[Optional[float], str, str, Optional[int], Optional[int], Optional[int]]] = {}

    n_reads = 0
    n_valid = 0
    sum_entropy = 0.0

    try:
        with open_text(path, "rt") as fin, open_text(out_path, "wt", compresslevel=compresslevel) as fout:
            if write_diag:
                fout.write(
                    "read_id\ttrue_taxid\tpred_taxid\tentropy\tlca_taxid\tlca_rank\tup_from_true\tdown_to_pred\tbranch_size_LCA\n"
                )
            else:
                fout.write("read_id\ttrue_taxid\tpred_taxid\tentropy\n")

            for line in fin:
                if not line.strip():
                    continue
                # Split only first 3 fields; last field is huge and not needed
                parts = line.rstrip("\n").split("\t", 3)
                if len(parts) < 3:
                    continue

                status = parts[0]
                read_id = parts[1]
                pred_taxid = parse_pred_taxid(parts[2], status=status)

                n_reads += 1

                H, L, rank_L, u, d, k = compute_cached_for_pred(
                    true_taxid=true_taxid,
                    pred_taxid=pred_taxid,
                    cache=cache,
                    true_lineage_set=true_lineage_set,
                    true_depth=true_depth,
                    alpha_up=alpha_up,
                    alpha_down=alpha_down,
                    unclassified_entropy=unclassified_entropy,
                    unclassified_sentinels=unclassified_sentinels,
                )

                ent_str = "" if H is None else str(H)
                if H is not None:
                    try:
                        sum_entropy += float(H)
                        n_valid += 1
                    except Exception:
                        pass

                if write_diag:
                    fout.write(
                        f"{read_id}\t{true_taxid}\t{pred_taxid}\t{ent_str}\t{L}\t{rank_L}\t"
                        f"{'' if u is None else u}\t{'' if d is None else d}\t{'' if k is None else k}\n"
                    )
                else:
                    fout.write(f"{read_id}\t{true_taxid}\t{pred_taxid}\t{ent_str}\n")

        mean_entropy = "" if n_valid == 0 else str(sum_entropy / n_valid)
        return {**job, "status": "ok", "n_reads": n_reads, "n_valid": n_valid, "mean_entropy": mean_entropy}

    except Exception as e:
        return {**job, "status": f"error: {e}", "n_reads": n_reads, "n_valid": n_valid, "mean_entropy": ""}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes-dmp", required=True)

    ap.add_argument(
        "--jobs-tsv",
        default=None,
        help="Optional: 2-column TSV (kraken_out_path<TAB>true_taxid). If set, --key-file/--kraken-dir are ignored.",
    )

    ap.add_argument("--key-file", required=False, help="taxid2filename.txt (TSV with filename/dataset/taxid columns)")
    ap.add_argument("--kraken-dir", required=False, help="Directory containing *_dbN.out files")

    ap.add_argument("--outdir", required=True)
    ap.add_argument("--summary-tsv", required=True)

    ap.add_argument("--recursive", action="store_true")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 1))

    ap.add_argument("--alpha-up", type=float, default=0.3)
    ap.add_argument("--alpha-down", type=float, default=1.0)
    ap.add_argument("--unclassified-entropy", type=float, default=None)

    ap.add_argument("--gzip", action="store_true", help="Write outputs as .gz")
    ap.add_argument("--compresslevel", type=int, default=3)
    ap.add_argument("--no-diagnostics", action="store_true")
    ap.add_argument("--skip-existing", action="store_true")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    missing_key: List[str] = []

    # Choose job source
    if args.jobs_tsv:
        jobs = read_jobs_tsv(args.jobs_tsv)
        if not jobs:
            raise SystemExit(f"No jobs found in {args.jobs_tsv}")
    else:
        if not args.key_file or not args.kraken_dir:
            raise SystemExit("Either provide --jobs-tsv OR provide both --key-file and --kraken-dir.")
        key_map = load_key(args.key_file)
        jobs, missing_key = discover_jobs(args.kraken_dir, key_map, recursive=args.recursive)
        if not jobs:
            raise SystemExit("No Kraken .out files discovered that match *_dbN.out and exist in key map.")

    # Load taxonomy once in parent (helps fork); workers will load if needed
    global PARENT, RANK, DEPTH, BRANCH_SIZE
    PARENT, RANK, DEPTH, BRANCH_SIZE = load_taxonomy(args.nodes_dmp)

    # Build task list
    tasks = []
    for j in jobs:
        stem = Path(j["path"]).name
        suffix = ".entropy.tsv.gz" if args.gzip else ".entropy.tsv"
        out_path = str(outdir / (stem + suffix))

        tasks.append(
            {
                **j,
                "out_path": out_path,
                "alpha_up": args.alpha_up,
                "alpha_down": args.alpha_down,
                "unclassified_entropy": args.unclassified_entropy,
                "write_diag": (not args.no_diagnostics),
                "compresslevel": args.compresslevel,
                "skip_existing": args.skip_existing,
            }
        )

    # Multiprocessing
    n_workers = max(1, int(args.jobs))
    start_methods = mp.get_all_start_methods()
    ctx = mp.get_context("fork") if "fork" in start_methods else mp.get_context("spawn")

    if n_workers == 1:
        results = [process_one(t) for t in tasks]
    else:
        with ctx.Pool(processes=n_workers, initializer=init_worker, initargs=(args.nodes_dmp,)) as pool:
            results = pool.map(process_one, tasks)

    # Write summary TSV
    with open(args.summary_tsv, "w") as f:
        f.write("path\tdataset\tfilename\tdb\ttrue_taxid\toutput\tstatus\tn_reads\tn_valid_entropy\tmean_entropy\n")
        for r in results:
            f.write(
                f"{r['path']}\t{r.get('dataset','')}\t{r.get('filename','')}\t{r.get('db','')}\t{r['true_taxid']}\t"
                f"{r['out_path']}\t{r['status']}\t{r.get('n_reads','')}\t{r.get('n_valid','')}\t{r.get('mean_entropy','')}\n"
            )

    # Optional: dump missing-key list (only applies in key/discover mode)
    if missing_key:
        miss_path = str(outdir / "missing_key_for_outputs.txt")
        with open(miss_path, "w") as f:
            for p in missing_key:
                f.write(p + "\n")


if __name__ == "__main__":
    main()
