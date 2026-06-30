#!/usr/bin/env python3
"""Reproduce and verify A006672 a(11) = r(C_4, K_{1,11}) = 16.

r(C_4, K_{1,n}) is the least m such that every red/blue edge-coloring of K_m
contains a red C_4 or a blue star K_{1,n}. Equivalently it is the least m for
which there is NO C_4-free graph on m vertices with minimum degree >= m - n
(take the red graph; blue-degree = m-1 - red-degree, so avoiding K_{1,n} means
blue-degree <= n-1, i.e. red-degree >= m - n).

This script certifies a(11) = 16 with both bounds, no solver required:

  LOWER  a(11) >= 16: exhibit a C_4-free graph on 15 vertices with min degree 4
         (so blue-degree <= 10 < 11). A good coloring of K_15 exists => a(11) > 15.
  UPPER  a(11) <= 16: Reiman/KST pair counting. A C_4-free graph on m vertices
         has every pair of vertices sharing <= 1 common neighbor, so
         sum_v C(deg v, 2) <= C(m, 2). Min degree d => m*C(d,2) <= C(m,2) =>
         d(d-1) <= m-1. For m=16 a min degree of 5 needs 20 <= 15: impossible,
         so every graph on 16 vertices has a vertex of degree <= 4, i.e.
         blue-degree >= 11 somewhere -> a blue K_{1,11} (or a red C_4). a(11) <= 16.
"""
from itertools import combinations

# 4-regular C_4-free witness on 15 vertices (labels 0..14), 30 edges.
WITNESS_EDGES = [
    (0,1),(0,2),(0,3),(0,4),(1,2),(1,5),(1,11),(2,6),(2,8),(3,4),
    (3,7),(3,12),(4,9),(4,10),(5,7),(5,9),(5,11),(6,7),(6,8),(6,10),
    (7,12),(8,9),(8,13),(9,13),(10,11),(10,14),(11,14),(12,13),(12,14),(13,14),
]


def check_witness(edges, n_vertices, n_star):
    adj = {v: set() for v in range(n_vertices)}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    # C_4-free  <=>  every pair of vertices has at most one common neighbor
    c4_free = all(len(adj[u] & adj[v]) <= 1 for u, v in combinations(range(n_vertices), 2))
    min_deg = min(len(adj[v]) for v in range(n_vertices))
    max_blue = (n_vertices - 1) - min_deg          # complement (blue) max degree
    return c4_free, min_deg, max_blue, max_blue < n_star


def reiman_forbids(m, d):
    """True iff a C_4-free graph on m vertices CANNOT have min degree d."""
    return d * (d - 1) > m - 1


if __name__ == "__main__":
    n = 11
    c4_free, min_deg, max_blue, ok_lower = check_witness(WITNESS_EDGES, 15, n)
    print("LOWER BOUND (witness on 15 vertices):")
    print(f"  C_4-free: {c4_free}   min red-degree: {min_deg}   "
          f"max blue-degree: {max_blue} (< {n}? {ok_lower})")
    lower_ok = c4_free and ok_lower
    print(f"  => a(11) >= 16: {lower_ok}")

    print("UPPER BOUND (Reiman/KST):")
    print(f"  C_4-free graph on 16 vertices with min degree 5 impossible? "
          f"{reiman_forbids(16, 5)}  (d(d-1)=20 > m-1=15)")
    upper_ok = reiman_forbids(16, 5)
    print(f"  => a(11) <= 16: {upper_ok}")

    print(f"\na(11) = 16 certified: {lower_ok and upper_ok}")
