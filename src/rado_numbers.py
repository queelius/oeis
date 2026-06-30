"""
SAT-based computation of Rado numbers for linear equations.

A Rado number R(L, k) for a linear equation L and k colors is the minimum N
such that every k-coloring of {1, ..., N} contains a monochromatic solution
to L.

We encode the linear equation as a list of integer coefficients [a1, a2, ..., am]
representing a1*x1 + a2*x2 + ... + am*xm = 0, where xi are positive integers
in {1, ..., N} and all xi receive the same color.  Variables are NOT required
to be distinct (following the standard Rado-theory convention).

For the SAT encoding, since we only care about which *elements* appear, we
convert each solution to the set of *distinct* elements involved.  The
monochromatic constraint then says: forbid every color c from being assigned
to ALL elements in that set simultaneously.

Encoding details:
  - Boolean variable x_{i,c} for element i in {1,...,N}, color c in {1,...,k}.
  - Each element gets exactly one color.
  - For every solution's distinct-element set S and every color c, add clause
    OR_{s in S} (NOT x_{s,c}).

The Rado number is the smallest N where the formula is UNSAT.
"""

from pysat.solvers import Solver


def _var(i: int, c: int, k: int) -> int:
    """Map element i (1-based) and color c (0-based) to SAT variable (1-based)."""
    return (i - 1) * k + c + 1


def _enumerate_solutions(coefficients: list[int], n: int,
                         distinct: bool = False) -> set[frozenset[int]]:
    """Enumerate all solutions to sum(a_i * x_i) = 0 with x_i in {1,...,n}.

    Returns a set of frozensets, each being the set of distinct elements
    participating in a solution.  De-duplicating at the set level avoids
    adding redundant SAT clauses.

    Args:
        coefficients: equation coefficients.
        n: universe size.
        distinct: if True, require all variables to take distinct values.
    """
    m = len(coefficients)
    result: set[frozenset[int]] = set()

    if m <= 1:
        return result

    last_coeff = coefficients[-1]
    if last_coeff == 0:
        _enumerate_brute(coefficients, n, distinct, result)
        return result

    def _recurse(depth: int, used: set[int] | None, partial_sum: int,
                 assignment: list[int]):
        if depth == m - 1:
            # Solve for the last variable: partial_sum + last_coeff * x_m = 0
            remainder = -partial_sum
            if remainder % last_coeff != 0:
                return
            x_last = remainder // last_coeff
            if x_last < 1 or x_last > n:
                return
            if distinct and x_last in used:
                return
            result.add(frozenset(assignment + [x_last]))
            return

        for x in range(1, n + 1):
            if distinct and x in used:
                continue
            new_sum = partial_sum + coefficients[depth] * x
            assignment.append(x)
            if distinct:
                used.add(x)
            _recurse(depth + 1, used, new_sum, assignment)
            if distinct:
                used.discard(x)
            assignment.pop()

    _recurse(0, set() if distinct else None, 0, [])
    return result


def _enumerate_brute(coefficients: list[int], n: int, distinct: bool,
                     result: set[frozenset[int]]):
    """Brute-force enumeration for degenerate cases (last coeff is 0)."""
    m = len(coefficients)

    def _recurse(depth, used, partial_sum, assignment):
        if depth == m:
            if partial_sum == 0:
                result.add(frozenset(assignment))
            return
        for x in range(1, n + 1):
            if distinct and x in used:
                continue
            assignment.append(x)
            if distinct:
                used.add(x)
            _recurse(depth + 1, used,
                     partial_sum + coefficients[depth] * x, assignment)
            if distinct:
                used.discard(x)
            assignment.pop()

    _recurse(0, set() if distinct else None, 0, [])


def _check_n(coefficients: list[int], k_colors: int, n: int,
             solver_name: str = "cd195",
             distinct: bool = False) -> tuple[bool, list[int] | None]:
    """Check whether N admits a k-coloring with no monochromatic solution.

    Returns (sat, coloring) where:
      - sat=True means a valid coloring exists (N is too small for the Rado number)
      - sat=False means no valid coloring exists (N >= the Rado number)
    """
    k = k_colors
    solution_sets = _enumerate_solutions(coefficients, n, distinct=distinct)

    with Solver(name=solver_name) as solver:
        # Each element gets exactly one color
        for i in range(1, n + 1):
            vars_i = [_var(i, c, k) for c in range(k)]
            solver.add_clause(vars_i)            # at-least-one
            for a in range(k):                   # at-most-one (pairwise)
                for b in range(a + 1, k):
                    solver.add_clause([-_var(i, a, k), -_var(i, b, k)])

        # Forbid each monochromatic solution
        for elem_set in solution_sets:
            elems = sorted(elem_set)
            for c in range(k):
                clause = [-_var(x, c, k) for x in elems]
                solver.add_clause(clause)

        if solver.solve():
            model = solver.get_model()
            coloring = [0] * (n + 1)
            for i in range(1, n + 1):
                for c in range(k):
                    if model[_var(i, c, k) - 1] > 0:
                        coloring[i] = c
                        break
            return True, coloring
        else:
            return False, None


def rado_number(coefficients: list[int], k_colors: int = 2,
                max_n: int = 500, solver_name: str = "cd195",
                distinct: bool = False,
                verbose: bool = False) -> tuple[int, list[int] | None]:
    """Find the Rado number R(L, k) for equation L with k colors.

    Args:
        coefficients: List of integer coefficients [a1, ..., am] for the
            homogeneous equation a1*x1 + ... + am*xm = 0.
        k_colors: Number of colors (default 2).
        max_n: Maximum N to search up to (default 500).
        solver_name: SAT solver name (default "cd195" for CaDiCaL).
        distinct: If True, require all variables in a solution to be distinct.
        verbose: Print progress if True.

    Returns:
        (rado_num, witness_coloring) where rado_num is the Rado number
        (smallest UNSAT N), or -1 if not found within max_n.
        witness_coloring is a valid coloring of {1,...,rado_num-1}.
    """
    # Phase 1: exponential probing to find first UNSAT
    lo = 0  # largest known SAT (0 means none found yet)
    hi = None
    witness = None

    probe = 1
    while probe <= max_n:
        sat, coloring = _check_n(coefficients, k_colors, probe, solver_name,
                                 distinct=distinct)
        if verbose:
            print(f"  N={probe}: {'SAT' if sat else 'UNSAT'}")
        if sat:
            lo = probe
            witness = coloring
            probe = min(probe * 2, max_n)
            if probe == lo:
                # We've hit max_n and it's still SAT
                return -1, None
        else:
            hi = probe
            break

    if hi is None:
        return -1, None

    # Phase 2: binary search between lo and hi
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        sat, coloring = _check_n(coefficients, k_colors, mid, solver_name,
                                 distinct=distinct)
        if verbose:
            print(f"  N={mid}: {'SAT' if sat else 'UNSAT'}")
        if sat:
            lo = mid
            witness = coloring
        else:
            hi = mid

    # hi is the Rado number (smallest UNSAT)
    # Make sure we have a witness for lo
    if witness is None and lo >= 1:
        _, witness = _check_n(coefficients, k_colors, lo, solver_name,
                              distinct=distinct)

    return hi, witness


def verify_coloring(coefficients: list[int], coloring: list[int], n: int,
                    distinct: bool = False) -> bool:
    """Verify that a coloring of {1,...,n} has no monochromatic solution.

    Args:
        coefficients: Equation coefficients.
        coloring: coloring[i] is the color of element i (1-indexed).
        n: Universe size.
        distinct: Whether solutions require distinct variables.

    Returns:
        True if the coloring avoids all monochromatic solutions.
    """
    solution_sets = _enumerate_solutions(coefficients, n, distinct=distinct)
    for elem_set in solution_sets:
        colors = {coloring[x] for x in elem_set}
        if len(colors) == 1:
            return False
    return True


def equation_from_form(pos_coeffs: list[int], rhs_coeff: int) -> list[int]:
    """Convert equation x1 + x2 + ... = rhs_coeff * z to coefficient form.

    Example: x + y = 3z becomes equation_from_form([1, 1], 3) -> [1, 1, -3]
    """
    return pos_coeffs + [-rhs_coeff]


# ---- Convenience wrappers for common equation families ----

def schur_number(k_colors: int, max_n: int = 500, **kwargs) -> tuple[int, list[int] | None]:
    """Compute Schur number S(k): minimum N s.t. every k-coloring of
    {1,...,N} has monochromatic x + y = z (variables may repeat)."""
    return rado_number([1, 1, -1], k_colors, max_n, distinct=False, **kwargs)


def rado_xy_kz(k_mult: int, k_colors: int = 2, max_n: int = 500,
               **kwargs) -> tuple[int, list[int] | None]:
    """Compute Rado number for x + y = k*z."""
    return rado_number([1, 1, -k_mult], k_colors, max_n, **kwargs)


def rado_sum_eq(num_vars_lhs: int, rhs_coeff: int = 1, k_colors: int = 2,
                max_n: int = 500, **kwargs) -> tuple[int, list[int] | None]:
    """Compute Rado number for x1 + ... + x_m = rhs_coeff * z."""
    coeffs = [1] * num_vars_lhs + [-rhs_coeff]
    return rado_number(coeffs, k_colors, max_n, **kwargs)


if __name__ == "__main__":
    import time

    print("=" * 60)
    print("Rado Number Computations")
    print("=" * 60)

    # Validation
    print("\n--- Validation: Known Rado Numbers ---")

    for k in range(1, 5):
        t0 = time.time()
        val, w = schur_number(k, max_n=200)
        dt = time.time() - t0
        print(f"S({k}) = {val}  [{dt:.2f}s]")

    t0 = time.time()
    val, w = rado_xy_kz(2, k_colors=2, max_n=100)
    dt = time.time() - t0
    print(f"R(x+y=2z, 2) = {val}  [{dt:.2f}s]  (W(3;2)=9)")

    # New computations
    print("\n--- New Computations: R(x+y=kz, 2) ---")
    for k_mult in range(3, 20):
        t0 = time.time()
        val, w = rado_xy_kz(k_mult, k_colors=2, max_n=500)
        dt = time.time() - t0
        if val > 0:
            print(f"R(x+y={k_mult}z, 2) = {val}  [{dt:.2f}s]")
        else:
            print(f"R(x+y={k_mult}z, 2) > 500  [{dt:.2f}s]")
