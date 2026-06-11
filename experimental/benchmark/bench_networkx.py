#!/usr/bin/env python3
"""
bench_networkx.py — benchmark nx.number_of_walks on a .mtx graph

Usage:
    python3 bench_networkx.py graph.mtx [k]

Notes:
    nx.number_of_walks always computes ALL pairs — there is no single-source
    mode in NetworkX. The all-pairs time is reported; single-source is shown
    as the all-pairs time (NetworkX offers no advantage for one source).

Output (stdout): TAB-separated lines prefixed with NXPY for shell parsing.
Output (stderr): human-readable progress.
"""

import sys
import time
import scipy.io
import scipy.sparse as sp
import networkx as nx

def now():
    return time.perf_counter()

def human(t):
    return f"{t:.4f} sec"

def main():
    if len(sys.argv) < 2:
        print("Usage: bench_networkx.py graph.mtx [k]", file=sys.stderr)
        sys.exit(1)

    mtx_file   = sys.argv[1]
    k          = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    graph_name = mtx_file.split("/")[-1].replace(".mtx", "")

    # ── load graph ────────────────────────────────────────────────────────────
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"NetworkX benchmark: {graph_name}  k={k}", file=sys.stderr)

    t0  = now()
    raw = scipy.io.mmread(mtx_file)
    mat = sp.csr_matrix(raw)
    try:
        G = nx.from_scipy_sparse_array(mat)
    except AttributeError:
        G = nx.from_scipy_sparse_matrix(mat)   # older NetworkX
    load_time = now() - t0

    n   = G.number_of_nodes()
    e   = G.number_of_edges()
    print(f"nodes : {n:,}   edges: {e:,}", file=sys.stderr)
    print(f"load  : {human(load_time)}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)

    # ── warmup ────────────────────────────────────────────────────────────────
    print("warmup ...", file=sys.stderr, end=" ", flush=True)
    t0 = now()
    nx.number_of_walks(G, k)
    print(f"{human(now() - t0)}", file=sys.stderr)

    # ── timed run: all-pairs ──────────────────────────────────────────────────
    print(f"timing nx.number_of_walks (k={k}) ...", file=sys.stderr,
          end=" ", flush=True)
    t0 = now()
    walks = nx.number_of_walks(G, k)
    t_all = now() - t0
    print(human(t_all), file=sys.stderr)

    # ── single-source note ────────────────────────────────────────────────────
    # NetworkX has no single-source API — number_of_walks always does all pairs.
    # We report the same all-pairs time for single-source to show the cost.
    print("\nNote: NetworkX has no single-source mode.", file=sys.stderr)
    print("      number_of_walks always computes all pairs.", file=sys.stderr)

    # ── stdout: machine-readable summary lines ────────────────────────────────
    print(f"NXPY\tall-pairs\t{graph_name}\tk={k}\t{t_all:.6f}")
    print(f"NXPY\tsingle-src\t{graph_name}\tk={k}\t{t_all:.6f}\tNOTE:same-as-all-pairs")

    print(f"\nSummary  graph={graph_name}  k={k}", file=sys.stderr)
    print(f"  nx.number_of_walks all-pairs : {human(t_all)}", file=sys.stderr)
    print(f"  single-source (no NX API)    : {human(t_all)}  (all-pairs cost)",
          file=sys.stderr)

if __name__ == "__main__":
    main()
