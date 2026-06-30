import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))


def _solver_available(solver_name):
    try:
        from pysat.solvers import Solver
        with Solver(name=solver_name) as s:
            s.add_clause([1])
            return s.solve()
    except Exception:
        return False


requires_sat = pytest.mark.skipif(
    not _solver_available("g4"),
    reason="python-sat with Glucose4 not available",
)
