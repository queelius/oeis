#!/usr/bin/env python3
"""Reproduce and verify A316632: Sprague-Grundy value of Node-Kayles on P_3 X P_n.

Node-Kayles on a graph G: players alternately pick a vertex not adjacent to any
previously picked vertex; the last player to move wins (normal play). A move at
vertex v deletes its closed neighborhood N[v] = {v} U neighbors(v) from the
remaining graph. The Sprague-Grundy value is mex over all moves; the value of a
disconnected position is the XOR of its components' values.

A316632(n) is this Grundy value for the 3 X n grid graph P_3 x P_n.

Independent from-scratch solver (frozenset states, connected-component splitting,
memoization). Reproduces all published terms a(1..16) and gives the new a(17)=2.
On 3 x 17 = 51 vertices this takes ~1-2 min.
"""
import sys
from functools import lru_cache

sys.setrecursionlimit(1_000_000)

PUBLISHED = [2, 1, 1, 0, 3, 3, 2, 2, 2, 3, 3, 5, 2, 4, 1, 3]  # a(1..16)


def grid_adj(n, rows=3):
    """Adjacency dict for the rows x n grid; vertex (r,c) -> c*rows + r."""
    def vid(r, c):
        return c * rows + r
    adj = {vid(r, c): set() for r in range(rows) for c in range(n)}
    for r in range(rows):
        for c in range(n):
            v = vid(r, c)
            if r + 1 < rows:
                w = vid(r + 1, c); adj[v].add(w); adj[w].add(v)
            if c + 1 < n:
                w = vid(r, c + 1); adj[v].add(w); adj[w].add(v)
    return adj


def grundy_node_kayles(n, rows=3):
    adj = grid_adj(n, rows)
    memo = {}

    def components(state):
        seen, comps = set(), []
        for s in state:
            if s in seen:
                continue
            stack, comp = [s], []
            seen.add(s)
            while stack:
                x = stack.pop(); comp.append(x)
                for y in adj[x]:
                    if y in state and y not in seen:
                        seen.add(y); stack.append(y)
            comps.append(frozenset(comp))
        return comps

    def g(state):
        if not state:
            return 0
        if state in memo:
            return memo[state]
        comps = components(state)
        if len(comps) > 1:                       # XOR of independent components
            val = 0
            for c in comps:
                val ^= g(c)
            memo[state] = val
            return val
        reach = set()
        for v in state:                          # one connected component
            removed = {v} | (adj[v] & state)
            reach.add(g(state - removed))
        m = 0
        while m in reach:
            m += 1
        memo[state] = m
        return m

    return g(frozenset(adj))


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 17
    ok = True
    for n in range(1, nmax + 1):
        v = grundy_node_kayles(n)
        if n <= len(PUBLISHED):
            match = (v == PUBLISHED[n - 1])
            ok = ok and match
            tag = "OK" if match else "MISMATCH"
            print(f"a({n:2d}) = {v}  [pub {PUBLISHED[n-1]}] {tag}")
        else:
            print(f"a({n:2d}) = {v}  <-- NEW")
    print("\nreproduces all published terms:", ok)
