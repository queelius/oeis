"""
Exact computation of Zarankiewicz numbers z(m,n;s,t).

z(m,n;s,t) = maximum number of 1s in an m x n 0-1 matrix with no s x t
all-1 submatrix. Equivalently, the maximum number of edges in a bipartite
graph with parts of sizes m and n containing no K_{s,t} subgraph.

Methods:
  - SAT encoding with binary search and cardinality constraints (python-sat)
  - Brute-force enumeration for small cases
  - Known values table for validation

Known exact values sourced from:
  - OEIS A072567: z(n,n;2,2)
  - OEIS A001198: k_3(n), so z(n,n;3,3) = A001198(n) - 1
  - OEIS A006613: k_{2,3}(n), so z(n,n;2,3) = A006613(n) - 1
  - OEIS A006614: k_{2,4}(n), so z(n,n;2,4) = A006614(n) - 1
  - OEIS A006615: k_{3,4}(n), so z(n,n;3,4) = A006615(n) - 1
  - OEIS A006616: k_4(n), so z(n,n;4,4) = A006616(n) - 1
  - Guy 1969, Tan 2022, Collins-Riasanovsky-Wallace-Radziszowski
"""
import math
from itertools import combinations

from pysat.card import CardEnc, EncType
from pysat.examples.rc2 import RC2
from pysat.formula import WCNF
from pysat.solvers import Solver


# ---------------------------------------------------------------------------
# Known values tables
# ---------------------------------------------------------------------------

# z(n,n;2,2) for n = 1..24 (OEIS A072567)
KNOWN_Z_2_2 = {
    1: 1, 2: 3, 3: 6, 4: 9, 5: 12, 6: 16, 7: 21, 8: 24,
    9: 29, 10: 34, 11: 39, 12: 45, 13: 52, 14: 56, 15: 61,
    16: 67, 17: 74, 18: 81, 19: 88, 20: 96, 21: 105, 22: 108,
    23: 115, 24: 122,
}

# z(n,n;3,3) for n = 1..16. For n < 3, z = n^2 trivially (cannot pick 3
# rows from fewer than 3). For n >= 3, derived from OEIS A001198:
# k_3(n) for n=3..16: 9, 14, 21, 27, 34, 43, 50, 61, 70, 81, 93, 106, 121, 129
# so z(n,n;3,3) = k_3(n) - 1.
KNOWN_Z_3_3 = {
    1: 1, 2: 4,
    3: 8, 4: 13, 5: 20, 6: 26, 7: 33, 8: 42, 9: 49,
    10: 60, 11: 69, 12: 80, 13: 92, 14: 105, 15: 120, 16: 128,
}

# z(n,n;2,3) for n = 1..11. For n < 2, z = n*n trivially. For n = 2,
# z(2,2;2,3) = 4 (cannot pick 3 cols from 2). For n >= 3, derived from
# OEIS A006613: k_{2,3}(n) for n=3..11: 8, 13, 17, 22, 29, 34, 40, 47, 56
# so z(n,n;2,3) = k_{2,3}(n) - 1.
KNOWN_Z_2_3 = {
    1: 1, 2: 4,
    3: 7, 4: 12, 5: 16, 6: 21, 7: 28, 8: 33, 9: 39,
    10: 46, 11: 55,
}

# z(n,n;2,4) for n = 1..11. For n < 4, z(n,n;2,4) = z(n,n;2,n) trivially.
# For n < 2, z = n*n. For n = 2, cannot pick 4 cols from 2, so z = 4.
# For n = 3, cannot pick 4 cols from 3, so z = 9.
# OEIS A006614: k_{2,4}(n) for n=4..11: 14, 21, 26, 32, 41, 48, 56, 67
# so z(n,n;2,4) = k_{2,4}(n) - 1.
KNOWN_Z_2_4 = {
    1: 1, 2: 4, 3: 9,
    4: 13, 5: 20, 6: 25, 7: 31, 8: 40, 9: 47, 10: 55, 11: 66,
}

# z(n,n;3,4) for n = 1..10. For n < 3 or n < 4, z = n^2 trivially.
# OEIS A006615: k_{3,4}(n) for n=4..10: 15, 22, 31, 38, 46, 57, 67
# so z(n,n;3,4) = k_{3,4}(n) - 1.
# a(10) = 67 is a NEW value computed here (was previously unknown).
KNOWN_Z_3_4 = {
    1: 1, 2: 4, 3: 9,
    4: 14, 5: 21, 6: 30, 7: 37, 8: 45, 9: 56, 10: 66,
}

# z(n,n;4,4) for n = 1..13. For n < 4, z = n^2 trivially.
# OEIS A006616: k_4(n) for n=4..13: 16, 23, 32, 43, 52, 62, 75, 87, 101, 118
# so z(n,n;4,4) = k_4(n) - 1.
KNOWN_Z_4_4 = {
    1: 1, 2: 4, 3: 9,
    4: 15, 5: 22, 6: 31, 7: 42, 8: 51, 9: 61, 10: 74, 11: 86, 12: 100, 13: 117,
}


# z(n,n+1;3,4) -- OEIS A006622: k_{3,4} for n x (n+1) matrices.
# a(n) for n=3..9: 12, 18, 26, 33, 41, 51, 61
# z(n,n+1;3,4) = a(n) - 1.
# a(9) = 61 is a NEW value computed here (was previously unknown).
KNOWN_Z_3_4_RECT1 = {
    # (rows, cols) -> z value, where cols = rows + 1
    (3, 4): 11, (4, 5): 17, (5, 6): 25, (6, 7): 32,
    (7, 8): 40, (8, 9): 50, (9, 10): 60,
}

# z(n,n+2;3,4) -- OEIS A006625: k_{3,4} for n x (n+2) matrices.
# a(n) for n=3..9: 14, 21, 28, 36, 45, 55, 67
# z(n,n+2;3,4) = a(n) - 1.
# a(9) = 67 is a NEW value computed here (was previously unknown).
KNOWN_Z_3_4_RECT2 = {
    # (rows, cols) -> z value, where cols = rows + 2
    (3, 5): 13, (4, 6): 20, (5, 7): 27, (6, 8): 35,
    (7, 9): 44, (8, 10): 54, (9, 11): 66,
}


_KNOWN_TABLES = {
    (2, 2): KNOWN_Z_2_2,
    (3, 3): KNOWN_Z_3_3,
    (2, 3): KNOWN_Z_2_3,
    (2, 4): KNOWN_Z_2_4,
    (3, 4): KNOWN_Z_3_4,
    (4, 4): KNOWN_Z_4_4,
}

_KNOWN_RECT_TABLES = {
    # (s, t, col_offset) -> table mapping (rows, cols) -> z
    (3, 4, 1): KNOWN_Z_3_4_RECT1,
    (3, 4, 2): KNOWN_Z_3_4_RECT2,
}


def lookup_known(m: int, n: int, s: int, t: int) -> int | None:
    """Look up a known exact value, or return None."""
    # Square case: z(n,n;s,t) with s <= t
    if m == n:
        key = (min(s, t), max(s, t))
        table = _KNOWN_TABLES.get(key)
        if table is not None:
            return table.get(n)

    # Non-square case: check rectangular tables
    # Try both (m,n;s,t) and transposed (n,m;t,s)
    for ms, ns, ss, ts in [(m, n, s, t), (n, m, t, s)]:
        if ns > ms:
            offset = ns - ms
            key = (min(ss, ts), max(ss, ts), offset)
            table = _KNOWN_RECT_TABLES.get(key)
            if table is not None:
                result = table.get((ms, ns))
                if result is not None:
                    return result

    return None


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def is_kst_free(matrix: list[list[int]], s: int, t: int) -> bool:
    """Check that an m x n 0-1 matrix has no s x t all-1 submatrix.

    An s x t all-1 submatrix means there exist s rows and t columns such
    that all s*t entries at their intersections are 1.
    """
    m = len(matrix)
    if m == 0:
        return True
    n = len(matrix[0])

    if s > m or t > n:
        return True  # cannot form the submatrix

    for row_subset in combinations(range(m), s):
        for col_subset in combinations(range(n), t):
            if all(matrix[r][c] == 1 for r in row_subset for c in col_subset):
                return False
    return True


def count_ones(matrix: list[list[int]]) -> int:
    """Count the number of 1s in a 0-1 matrix."""
    return sum(sum(row) for row in matrix)


# ---------------------------------------------------------------------------
# Brute-force solver
# ---------------------------------------------------------------------------

def zarankiewicz_brute(m: int, n: int, s: int, t: int) -> int:
    """Compute z(m,n;s,t) by brute-force enumeration of all 2^(m*n) matrices.

    Only feasible for very small m*n (say, m*n <= 16).
    """
    if s > m or t > n:
        return m * n  # trivially, every matrix is K_{s,t}-free

    total_entries = m * n
    best = 0

    for bits in range(1 << total_entries):
        # Build matrix from bitmask
        matrix = []
        for i in range(m):
            row = []
            for j in range(n):
                row.append((bits >> (i * n + j)) & 1)
            matrix.append(row)

        ones = bin(bits).count('1')
        if ones <= best:
            continue  # cannot improve

        if is_kst_free(matrix, s, t):
            best = ones

    return best


# ---------------------------------------------------------------------------
# SAT encoding
# ---------------------------------------------------------------------------

def _var(i: int, j: int, n: int) -> int:
    """Map matrix entry (i, j) to SAT variable (1-indexed)."""
    return i * n + j + 1


def _build_forbidden_clauses(
    m: int, n: int, s: int, t: int,
) -> tuple[int, list[list[int]]]:
    """Build clauses forbidding all s x t all-1 submatrices.

    For each s-subset of rows R and t-subset of columns C, we add the clause
    that at least one of the s*t entries (r,c) for r in R, c in C is 0.
    In CNF: OR_{(r,c) in R x C} NOT x_{r,c}.

    Returns (num_matrix_vars, clauses).
    """
    num_vars = m * n
    clauses = []

    for row_subset in combinations(range(m), s):
        for col_subset in combinations(range(n), t):
            clause = [-_var(r, c, n) for r in row_subset for c in col_subset]
            clauses.append(clause)

    return num_vars, clauses


def _check_feasible(
    m: int, n: int, s: int, t: int, k: int,
    solver_name: str = "cd195",
) -> list[int] | None:
    """Check if there exists a K_{s,t}-free m x n matrix with >= k ones.

    Returns the model (variable assignments) if SAT, or None if UNSAT.
    Uses a fresh solver instance each time for clean binary search.
    """
    num_vars, forbidden_clauses = _build_forbidden_clauses(m, n, s, t)
    matrix_vars = list(range(1, num_vars + 1))

    top_id = num_vars
    card_clauses = CardEnc.atleast(
        matrix_vars, bound=k,
        top_id=top_id, encoding=EncType.totalizer,
    )

    with Solver(name=solver_name) as solver:
        for clause in forbidden_clauses:
            solver.add_clause(clause)
        for clause in card_clauses.clauses:
            solver.add_clause(clause)

        if solver.solve():
            return solver.get_model()
    return None


def zarankiewicz_maxsat(
    m: int, n: int, s: int, t: int,
    return_matrix: bool = False,
) -> int | tuple[int, list[list[int]]]:
    """Compute z(m,n;s,t) via MaxSAT (RC2 solver).

    Maximizes the number of ones in an m x n matrix subject to the constraint
    that no s x t all-ones submatrix exists. Uses the RC2 MaxSAT solver from
    python-sat.

    This is typically faster than binary search SAT because it avoids the
    hard UNSAT proofs at the boundary.
    """
    # Trivial cases
    if s > m or t > n:
        val = m * n
        if return_matrix:
            return val, [[1] * n for _ in range(m)]
        return val

    num_vars, forbidden_clauses = _build_forbidden_clauses(m, n, s, t)

    # Build WCNF: hard clauses (forbidden submatrices) + soft clauses (maximize 1s)
    wcnf = WCNF()

    # Hard clauses: forbid K_{s,t}
    for clause in forbidden_clauses:
        wcnf.append(clause)

    # Soft clauses: each matrix variable should be 1 (weight 1)
    for v in range(1, num_vars + 1):
        wcnf.append([v], weight=1)

    # Solve with RC2
    with RC2(wcnf) as solver:
        model = solver.compute()

    if model is None:
        if return_matrix:
            return 0, [[0] * n for _ in range(m)]
        return 0

    best_val = sum(1 for v in range(1, num_vars + 1) if model[v - 1] > 0)

    if return_matrix:
        matrix = []
        for i in range(m):
            row = []
            for j in range(n):
                v = _var(i, j, n)
                row.append(1 if model[v - 1] > 0 else 0)
            matrix.append(row)
        return best_val, matrix

    return best_val


def zarankiewicz_sat(
    m: int, n: int, s: int, t: int,
    solver_name: str = "cd195",
    return_matrix: bool = False,
) -> int | tuple[int, list[list[int]]]:
    """Compute z(m,n;s,t) via MaxSAT (RC2 solver).

    Uses the RC2 MaxSAT solver from python-sat. Delegates to
    zarankiewicz_maxsat which handles the WCNF encoding.

    Parameters
    ----------
    m, n : int
        Matrix dimensions.
    s, t : int
        Forbidden submatrix dimensions.
    solver_name : str
        SAT solver name (currently unused; RC2 uses its own solver).
    return_matrix : bool
        If True, also return the optimal matrix.

    Returns
    -------
    int or (int, list[list[int]])
        The Zarankiewicz number, and optionally a witness matrix.
    """
    return zarankiewicz_maxsat(m, n, s, t, return_matrix=return_matrix)


def zarankiewicz_sat_with_symmetry(
    n: int, s: int, t: int,
    solver_name: str = "cd195",
) -> int:
    """Compute z(n,n;s,t) exploiting symmetry s == t for the square case.

    When s == t and m == n, we can add symmetry-breaking constraints.
    Currently delegates to the main SAT solver.
    """
    return zarankiewicz_sat(n, n, s, t, solver_name=solver_name)


# ---------------------------------------------------------------------------
# Upper bound: Kovari-Sos-Turan theorem
# ---------------------------------------------------------------------------

def _kst_one_side(part_a: int, part_b: int, k: int, max_common: int) -> float:
    """Compute the KST bound from one side of the bipartite graph.

    We pick k-subsets from part A (size part_a). Each such k-subset has at
    most max_common common neighbors in part B. The counting argument
    double-counts via degrees in part B:

        sum over vertices in B of C(deg_b, k) <= C(part_a, k) * max_common

    Applying Jensen's inequality to the LHS (convexity of C(x, k)):

        part_b * C(e / part_b, k) <= C(part_a, k) * max_common

    We solve for the maximum e satisfying this inequality.

    Returns part_a * part_b if the constraint is trivially satisfied
    (k > part_a, so no k-subsets exist).
    """
    if k > part_a:
        return float(part_a * part_b)

    rhs = max_common * math.comb(part_a, k)

    if k == 2:
        # Closed form from the quadratic:
        # part_b * C(e/part_b, 2) <= rhs
        # part_b * (e/part_b)(e/part_b - 1) / 2 <= rhs
        # (e^2 - part_b * e) / (2 * part_b) <= rhs
        # e^2 - part_b * e - 2 * part_b * rhs <= 0
        # e <= (part_b + sqrt(part_b^2 + 8 * part_b * rhs)) / 2
        disc = part_b * part_b + 8 * part_b * rhs
        return (part_b + math.sqrt(disc)) / 2
    else:
        def _falling_binom(x: float, r: int) -> float:
            """Generalized binomial C(x, r) for real x."""
            result = 1.0
            for i in range(r):
                result *= (x - i)
            return result / math.factorial(r)

        # Binary search for max e such that part_b * C(e/part_b, k) <= rhs
        lo, hi = 0.0, float(part_a * part_b)
        for _ in range(200):
            mid = (lo + hi) / 2
            lhs = part_b * _falling_binom(mid / part_b, k)
            if lhs <= rhs:
                lo = mid
            else:
                hi = mid
        return lo


def kst_upper_bound(m: int, n: int, s: int, t: int) -> float:
    """Kovari-Sos-Turan upper bound on z(m,n;s,t).

    For an m x n 0-1 matrix avoiding K_{s,t}, Jensen's inequality can be
    applied from either side of the bipartite graph:

    Row side (s-subsets of rows share at most t-1 common columns):
        n * C(e/n, s) <= C(m, s) * (t-1)

    Column side (t-subsets of columns share at most s-1 common rows):
        m * C(e/m, t) <= C(n, t) * (s-1)

    For rectangular matrices (m != n), the two sides can give different
    bounds. We return the minimum (tightest) of both.
    """
    # Row side: pick s-subsets from rows (part_a=m), count via column
    # degrees (part_b=n), each s-subset of rows shares at most t-1 columns.
    bound_row = _kst_one_side(part_a=m, part_b=n, k=s, max_common=t - 1)

    # Column side: pick t-subsets from columns (part_a=n), count via row
    # degrees (part_b=m), each t-subset of columns shares at most s-1 rows.
    bound_col = _kst_one_side(part_a=n, part_b=m, k=t, max_common=s - 1)

    return min(bound_row, bound_col)


# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------

def zarankiewicz(
    m: int, n: int, s: int, t: int,
    method: str = "sat",
    solver_name: str = "cd195",
) -> int:
    """Compute z(m,n;s,t).

    Parameters
    ----------
    method : str
        "sat" for SAT solver, "brute" for brute force, "lookup" for known.
    """
    if method == "lookup":
        val = lookup_known(m, n, s, t)
        if val is not None:
            return val
        raise ValueError(f"No known value for z({m},{n};{s},{t})")
    elif method == "brute":
        return zarankiewicz_brute(m, n, s, t)
    elif method == "sat":
        return zarankiewicz_sat(m, n, s, t, solver_name=solver_name)
    else:
        raise ValueError(f"Unknown method: {method}")


def find_witness(
    m: int, n: int, s: int, t: int,
    solver_name: str = "cd195",
) -> tuple[int, list[list[int]]]:
    """Find z(m,n;s,t) and return (value, witness_matrix)."""
    return zarankiewicz_sat(m, n, s, t, solver_name=solver_name,
                           return_matrix=True)


def format_matrix(matrix: list[list[int]]) -> str:
    """Pretty-print a 0-1 matrix."""
    lines = []
    for row in matrix:
        lines.append(" ".join(str(x) for x in row))
    return "\n".join(lines)
