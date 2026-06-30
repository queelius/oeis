"""
Sandpile groups (critical groups, Jacobians) of finite connected graphs.

For a connected graph G on n vertices, the sandpile group K(G) is the
cokernel of the reduced Laplacian matrix tilde L:

    K(G) = Z^{n-1} / image(tilde L)

where tilde L is obtained from the Laplacian L = D - A by deleting one row
and one column (the choice of "root" does not affect the abstract group).
The Smith normal form of tilde L yields the invariant-factor decomposition

    K(G) ~= Z/d_1 x Z/d_2 x ... x Z/d_{n-1},   d_1 | d_2 | ... | d_{n-1}.

Closed forms for selected families:

  * Cycles:           K(C_n) = Z/n
  * Complete:         K(K_n) = (Z/n)^{n-2}
  * Complete bipart.: K(K_{m,n}) = Z/(mn) x (Z/m)^{n-2} x (Z/n)^{m-2}
                      (for m, n >= 2; trivial factors omitted)
  * Trees:            K(T) = trivial (one spanning tree)
  * Wheels W_n        |K(W_n)| = L_{2n} - 2 (Lucas numbers); structure
                      involves Lucas/Fibonacci patterns (Biggs 1999).

By the Matrix-Tree Theorem,

    |K(G)| = number of spanning trees of G = det(tilde L).

References:
  - Lorenzini 1991 ("A finite group attached to the Laplacian of a graph")
  - Biggs 1997, 1999 (Chip-firing and the critical group of a graph)
  - Bacher, de la Harpe, Nagnibeda 1997 (lattice of integer flows)
  - Baker, Norine 2007 (Riemann-Roch and Abel-Jacobi theory on a graph)
  - Levine, Propp 2010 (Sandpiles and the harmonic representation)
  - Klivans 2018 (The Mathematics of Chip-Firing, CRC Press)

Implementation notes:
  * We work entirely over Z with Python big-ints; the Smith normal form
    routine never sees floating point.
  * For each graph, the SNF computation is O(n^4) bit operations in the
    worst case but typically far less for sparse Laplacians on n <= 8.
  * We use sympy's Matrix.smith_normal_form when available (exact, well
    tested) and fall back to a Bareiss-style integer SNF.
"""

from __future__ import annotations

from typing import Iterable, Sequence


# --- Adjacency / Laplacian -------------------------------------------------

def adjacency_matrix(n: int, edges: Iterable[tuple[int, int]]) -> list[list[int]]:
    """Return the n x n adjacency matrix of a simple undirected graph.

    Self-loops are ignored.  Duplicate edges are counted once.
    """
    A = [[0] * n for _ in range(n)]
    seen: set[tuple[int, int]] = set()
    for u, v in edges:
        if u == v:
            continue
        e = (min(u, v), max(u, v))
        if e in seen:
            continue
        seen.add(e)
        A[u][v] = 1
        A[v][u] = 1
    return A


def laplacian(adj: Sequence[Sequence[int]]) -> list[list[int]]:
    """Return the (combinatorial) Laplacian L = D - A as a 2D list of ints."""
    n = len(adj)
    L = [[-int(adj[i][j]) if i != j else 0 for j in range(n)] for i in range(n)]
    for i in range(n):
        L[i][i] = sum(int(adj[i][j]) for j in range(n) if j != i)
    return L


def reduced_laplacian(adj: Sequence[Sequence[int]], root: int = 0) -> list[list[int]]:
    """Return the reduced Laplacian (delete row and column `root`) as ints."""
    n = len(adj)
    if not (0 <= root < n):
        raise ValueError(f"root {root} out of range for n={n}")
    L = laplacian(adj)
    return [
        [L[i][j] for j in range(n) if j != root]
        for i in range(n) if i != root
    ]


# --- Smith normal form ------------------------------------------------------

def _gcd(a: int, b: int) -> int:
    a, b = abs(int(a)), abs(int(b))
    while b:
        a, b = b, a % b
    return a


def _smith_fallback(M: list[list[int]]) -> list[int]:
    """Compute Smith normal form invariant factors of an integer matrix.

    This is a row/column reduction over Z that does not require sympy.
    The result is the list of diagonal entries (without trailing zeros)
    in divisibility order d_1 | d_2 | ... | d_r.

    Algorithm: at each pivot position (k, k), reduce row k and column k
    to have a single nonzero entry equal to the gcd of the remaining
    minor (this is Bareiss-style integer Gauss-Jordan).  Then iterate
    on the (k+1)-th block.  Finally, propagate divisibility by the
    standard "divisibility sweep" (replace consecutive pairs (a, b)
    not satisfying a | b by (gcd, lcm)).
    """
    if not M:
        return []
    M = [row[:] for row in M]
    nrows = len(M)
    ncols = len(M[0])
    k = 0
    while k < min(nrows, ncols):
        # Find a nonzero pivot in the (k:, k:) submatrix; if none, stop.
        # We pick the entry with smallest absolute nonzero value to avoid
        # blowup and rotate it to position (k, k).
        pi, pj, pv = -1, -1, 0
        for i in range(k, nrows):
            for j in range(k, ncols):
                v = M[i][j]
                if v != 0 and (pv == 0 or abs(v) < abs(pv)):
                    pi, pj, pv = i, j, v
        if pv == 0:
            break
        if pi != k:
            M[k], M[pi] = M[pi], M[k]
        if pj != k:
            for r in range(nrows):
                M[r][k], M[r][pj] = M[r][pj], M[r][k]
        # Reduce row k and column k modulo the pivot until both have only
        # the pivot nonzero in the (k:, k:) block.  Each such reduction
        # strictly decreases either the pivot magnitude or the count of
        # nonzero entries off the pivot, so the loop terminates.
        progress = True
        while progress:
            progress = False
            # Column k: zero out below-pivot entries.
            for i in range(k + 1, nrows):
                if M[i][k] == 0:
                    continue
                q = M[i][k] // M[k][k]
                for j in range(k, ncols):
                    M[i][j] -= q * M[k][j]
                if M[i][k] != 0:
                    # Swap to use this as new pivot if it's smaller.
                    if abs(M[i][k]) < abs(M[k][k]):
                        M[k], M[i] = M[i], M[k]
                    progress = True
            # Row k: zero out right-of-pivot entries.
            for j in range(k + 1, ncols):
                if M[k][j] == 0:
                    continue
                q = M[k][j] // M[k][k]
                for i in range(k, nrows):
                    M[i][j] -= q * M[i][k]
                if M[k][j] != 0:
                    if abs(M[k][j]) < abs(M[k][k]):
                        for r in range(nrows):
                            M[r][k], M[r][j] = M[r][j], M[r][k]
                    progress = True
        # Make the pivot positive.
        if M[k][k] < 0:
            for j in range(k, ncols):
                M[k][j] = -M[k][j]
        # Ensure pivot divides every entry in the rest of the matrix; if
        # not, mix it in via row operation and repeat.
        bad = False
        for i in range(k + 1, nrows):
            for j in range(k + 1, ncols):
                if M[i][j] % M[k][k] != 0:
                    # Bring row i into row k.
                    for jj in range(k, ncols):
                        M[k][jj] += M[i][jj]
                    bad = True
                    break
            if bad:
                break
        if bad:
            continue  # repeat reduction at position k
        k += 1

    # Collect diagonal entries.
    r = min(nrows, ncols)
    diag = [abs(M[i][i]) for i in range(r) if M[i][i] != 0]
    # Final divisibility sweep (insurance; should already hold).
    changed = True
    while changed:
        changed = False
        for i in range(len(diag) - 1):
            a, b = diag[i], diag[i + 1]
            if b % a != 0:
                g = _gcd(a, b)
                lcm = (a // g) * b
                diag[i], diag[i + 1] = g, lcm
                changed = True
    return diag


def smith_normal_form(M: Sequence[Sequence[int]]) -> list[int]:
    """Return the (sorted, divisibility-ordered) invariant factors of M.

    Uses sympy's exact Smith normal form when available; otherwise falls
    back to the integer reduction in `_smith_fallback`.  The output is the
    list of nonzero diagonal entries d_1 | d_2 | ... | d_r of SNF(M).
    """
    if not M or not M[0]:
        return []
    try:
        import sympy  # type: ignore
    except Exception:
        sympy = None

    if sympy is not None:
        S = sympy.Matrix([[int(x) for x in row] for row in M]).rref()  # noqa: F841
        # Use Matrix.smith_normal_form via matrices.normalforms.
        from sympy.matrices.normalforms import smith_normal_form as _snf
        Mi = sympy.Matrix([[int(x) for x in row] for row in M])
        S = _snf(Mi, domain=sympy.ZZ)
        diag = []
        nrows = S.shape[0]
        ncols = S.shape[1]
        for i in range(min(nrows, ncols)):
            v = int(abs(S[i, i]))
            if v != 0:
                diag.append(v)
        return diag

    M_list = [[int(x) for x in row] for row in M]
    return _smith_fallback(M_list)


# --- Sandpile group ---------------------------------------------------------

def sandpile_group(adj: Sequence[Sequence[int]], root: int = 0) -> dict:
    """Return the sandpile group of a connected graph in a structured form.

    The result is a dict with keys:
        n             -- number of vertices
        order         -- |K(G)| = number of spanning trees
        invariant_factors -- list of integers d_1 | ... | d_{n-1}, with
                             trailing 1s included so the length is n-1.
        nontrivial_factors -- the same list with the leading 1's stripped,
                              i.e. only the d_i > 1
        rank          -- number of nontrivial invariant factors
    """
    n = len(adj)
    if n == 0:
        return {"n": 0, "order": 1, "invariant_factors": [],
                "nontrivial_factors": [], "rank": 0}
    if n == 1:
        return {"n": 1, "order": 1, "invariant_factors": [],
                "nontrivial_factors": [], "rank": 0}

    Lr = reduced_laplacian(adj, root=root)
    diag = smith_normal_form(Lr)

    # SNF of an (n-1) x (n-1) full-rank integer matrix has exactly n-1
    # nonzero invariant factors.  Pad with leading 1s if `_smith_fallback`
    # returned fewer (which would only happen for non-connected graphs).
    expected = n - 1
    while len(diag) < expected:
        diag.insert(0, 1)
    # Trim if SNF returned more (should not happen with full rank).
    diag = diag[-expected:]

    order = 1
    for d in diag:
        order *= d
    nontrivial = [d for d in diag if d > 1]
    return {
        "n": n,
        "order": order,
        "invariant_factors": diag,
        "nontrivial_factors": nontrivial,
        "rank": len(nontrivial),
    }


def num_spanning_trees(adj: Sequence[Sequence[int]]) -> int:
    """Return the number of spanning trees via det(tilde L) (Matrix-Tree).

    Uses an integer-arithmetic determinant (Bareiss) so the result is
    exact for graphs of any size.  Mainly a sanity check for the SNF
    computation: |K(G)| must equal num_spanning_trees(G).
    """
    n = len(adj)
    if n == 0:
        return 0
    if n == 1:
        return 1
    M = reduced_laplacian(adj, root=0)
    return _det_bareiss(M)


def _det_bareiss(M: list[list[int]]) -> int:
    """Exact integer determinant via Bareiss' fraction-free elimination."""
    n = len(M)
    if n == 0:
        return 1
    A = [row[:] for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        # Pivot: find first nonzero in column k from row k onward.
        if A[k][k] == 0:
            for r in range(k + 1, n):
                if A[r][k] != 0:
                    A[k], A[r] = A[r], A[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
            A[i][k] = 0
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


# --- Group description helpers ---------------------------------------------

def group_signature(invariant_factors: Sequence[int]) -> str:
    """Pretty-print a sandpile group, e.g. 'Z/2 x Z/12'.

    Trivial factors (1) are omitted.  An empty signature is shown as 'Z/1'
    (the trivial group).
    """
    nz = [d for d in invariant_factors if d > 1]
    if not nz:
        return "Z/1"
    return " x ".join(f"Z/{d}" for d in nz)


def is_cyclic(invariant_factors: Sequence[int]) -> bool:
    """A finite abelian group is cyclic iff its invariant-factor list has
    at most one entry > 1."""
    return sum(1 for d in invariant_factors if d > 1) <= 1


# --- Convenience constructors for named families ---------------------------

def cycle_graph(n: int) -> list[list[int]]:
    """Adjacency matrix of the cycle C_n (n >= 3)."""
    if n < 3:
        raise ValueError("cycle requires n >= 3")
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        j = (i + 1) % n
        A[i][j] = A[j][i] = 1
    return A


def path_graph(n: int) -> list[list[int]]:
    """Adjacency matrix of the path P_n (n >= 1)."""
    A = [[0] * n for _ in range(n)]
    for i in range(n - 1):
        A[i][i + 1] = A[i + 1][i] = 1
    return A


def complete_graph(n: int) -> list[list[int]]:
    """Adjacency matrix of K_n."""
    A = [[1 if i != j else 0 for j in range(n)] for i in range(n)]
    return A


def complete_bipartite(m: int, n: int) -> list[list[int]]:
    """Adjacency matrix of K_{m,n}.  Vertex order: m left then n right."""
    N = m + n
    A = [[0] * N for _ in range(N)]
    for i in range(m):
        for j in range(m, N):
            A[i][j] = A[j][i] = 1
    return A


def wheel_graph(n: int) -> list[list[int]]:
    """Adjacency matrix of the wheel W_n: a hub connected to a C_n."""
    if n < 3:
        raise ValueError("wheel requires n >= 3 for the rim")
    N = n + 1
    A = [[0] * N for _ in range(N)]
    # Rim cycle on 0..n-1; hub = n.
    for i in range(n):
        j = (i + 1) % n
        A[i][j] = A[j][i] = 1
        A[i][n] = A[n][i] = 1
    return A


def prism_graph(n: int) -> list[list[int]]:
    """Adjacency matrix of the prism Y_n = C_n x K_2 on 2n vertices."""
    if n < 3:
        raise ValueError("prism requires n >= 3")
    N = 2 * n
    A = [[0] * N for _ in range(N)]
    for i in range(n):
        # Two parallel cycles.
        A[i][(i + 1) % n] = A[(i + 1) % n][i] = 1
        A[n + i][n + (i + 1) % n] = A[n + (i + 1) % n][n + i] = 1
        # Rungs.
        A[i][n + i] = A[n + i][i] = 1
    return A


# --- Catalog enumeration ----------------------------------------------------

def _adj_from_edges(n: int, edges: Iterable[tuple[int, int]]) -> list[list[int]]:
    return adjacency_matrix(n, edges)


def _is_connected_adj(adj: Sequence[Sequence[int]]) -> bool:
    n = len(adj)
    if n == 0:
        return True
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in range(n):
            if adj[u][v] and v not in seen:
                seen.add(v)
                stack.append(v)
    return len(seen) == n


def enumerate_connected_graphs(n: int):
    """Yield (n, sorted_edges) for each isomorphism class of connected
    simple graph on n vertices.  Uses pynauty canonical certificates for
    deduplication.  Iterative-extension strategy keeps memory modest.

    Returns adjacency matrices via _adj_from_edges; consumers can call
    sandpile_group directly on the yielded edges.
    """
    try:
        import pynauty  # type: ignore
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "pynauty is required for graph enumeration; install with "
            "`pip install -e \".[enumerate]\"`"
        ) from e

    if n <= 0:
        return
    if n == 1:
        yield (1, [])
        return

    def cert(n_, edges):
        adj_dict = {i: [] for i in range(n_)}
        for u, v in edges:
            adj_dict[u].append(v)
            adj_dict[v].append(u)
        adj_dict = {k: sorted(set(vs)) for k, vs in adj_dict.items()}
        g = pynauty.Graph(n_, directed=False, adjacency_dict=adj_dict)
        return pynauty.certificate(g)

    # Build all (not necessarily connected) graphs on n vertices by extension
    # from n-1, then filter to connected.  This is faster than enumerating
    # 2^C(n,2) graphs for n=8 because the n-1 representatives are few.
    if n == 2:
        small_reps: list[tuple[int, list[tuple[int, int]]]] = [(1, [])]
    else:
        small_reps = list(enumerate_all_graphs(n - 1))

    seen: set = set()
    new_v = n - 1
    for (_, e_small) in small_reps:
        for mask in range(1 << (n - 1)):
            extra = [(u, new_v) for u in range(n - 1) if (mask >> u) & 1]
            edges = sorted(set(tuple(sorted(e)) for e in list(e_small) + extra))
            adj = _adj_from_edges(n, edges)
            if not _is_connected_adj(adj):
                continue
            c = cert(n, edges)
            if c in seen:
                continue
            seen.add(c)
            yield (n, edges)


def enumerate_all_graphs(n: int):
    """Yield (n, sorted_edges) for every isomorphism class of simple graph
    on n vertices (connected or not)."""
    try:
        import pynauty  # type: ignore
    except Exception as e:  # pragma: no cover
        raise RuntimeError("pynauty required") from e

    if n == 0:
        yield (0, [])
        return
    if n == 1:
        yield (1, [])
        return

    def cert(n_, edges):
        adj_dict = {i: [] for i in range(n_)}
        for u, v in edges:
            adj_dict[u].append(v)
            adj_dict[v].append(u)
        adj_dict = {k: sorted(set(vs)) for k, vs in adj_dict.items()}
        g = pynauty.Graph(n_, directed=False, adjacency_dict=adj_dict)
        return pynauty.certificate(g)

    seen: set = set()
    if n == 1:
        yield (1, [])
        return
    if n == 2:
        for edges in [[], [(0, 1)]]:
            yield (2, edges)
        return

    small = list(enumerate_all_graphs(n - 1))
    new_v = n - 1
    for (_, e_small) in small:
        for mask in range(1 << (n - 1)):
            extra = [(u, new_v) for u in range(n - 1) if (mask >> u) & 1]
            edges = sorted(set(tuple(sorted(e)) for e in list(e_small) + extra))
            c = cert(n, edges)
            if c in seen:
                continue
            seen.add(c)
            yield (n, edges)


def catalog_record(n: int, edges: Sequence[tuple[int, int]]) -> dict:
    """Build one JSON-friendly catalog row for a connected graph."""
    adj = _adj_from_edges(n, edges)
    K = sandpile_group(adj, root=0)
    return {
        "n": n,
        "m": len(edges),
        "edges": [[int(u), int(v)] for u, v in edges],
        "order": K["order"],
        "invariant_factors": K["invariant_factors"],
        "nontrivial_factors": K["nontrivial_factors"],
        "rank": K["rank"],
        "signature": group_signature(K["invariant_factors"]),
        "is_cyclic": is_cyclic(K["invariant_factors"]),
        "is_tree": K["order"] == 1,
    }
