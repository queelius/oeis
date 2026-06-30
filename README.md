# OEIS

My working home for OEIS-focused research: a reusable toolkit (`src/`) for
finding, extending, and verifying integer sequences; the portfolio of
contributions to date; and a backlog of new targets to attack. This repo is
**self-contained and primary**. The value-generating code originated in my
research repos ([open-problems](https://github.com/queelius/open-problems),
[computational-explorations](https://github.com/queelius/computational-explorations))
and now lives and evolves here independently; the two do not sync. Author
search: [oeis.org/search?q=alex+towell](https://oeis.org/search?q=alex+towell).
A public, narrative version of this portfolio lives at
[metafunctor.com/research/oeis-contributions](https://metafunctor.com/research/oeis-contributions/).

**Why these.** These sequences have nothing in common at the level of subject:
they range over Ramsey theory, extremal graph theory, knot theory, additive
combinatorics, Boolean function complexity, combinatorial game theory, and the
algebra of graphs. What they share is how they were
found. Each is a quantity with no known formula, computed exactly by SAT
solving or exhaustive search, pushed one value past the published frontier, and
checked against all prior art before it was submitted. The spread across fields
is deliberate: it is evidence that the method travels. The computed values are
where I go looking for structure; the Rado numbers began as a column of
computed integers and ended as a proved theorem, and the rest are trailheads.

## Layout

- `contributions.yaml` : single source of truth (one record per sequence).
- `sequences/<id>/`     : per-sequence `metadata.yaml` and a standalone
  `reproduce.py` that recomputes the terms.
- `src/`                : copied value-generating modules (see `src/PROVENANCE.md`).
- `bfiles/`             : canonical b-files.
- `scripts/`            : `validate.py` (lint the index), `render_readme.py`
  (regenerate the tables below from `contributions.yaml`).

## Reproduce

```
python3 scripts/validate.py                 # lint the index
python3 sequences/A006622/reproduce.py      # recompute one sequence's terms
pytest                                       # run the copied generators' tests
```

## Contributions

<!-- AUTOGEN:contributions BEGIN -->

### Authored (new sequences)

| Sequence | Description | Provenance | Status |
|----------|-------------|------------|--------|
| [A394445](https://oeis.org/A394445) | Distinct-variable 2-color Rado numbers for x+y=n*z | [open-problems/src/rado_numbers.py](https://github.com/queelius/open-problems/blob/main/src/rado_numbers.py) | published |
| [A394661](https://oeis.org/A394661) | Triangle T(n,k): number of prime knots with n crossings and three-genus k | [open-problems/discoveries/topology/knot-invariants](https://github.com/queelius/open-problems/blob/main/discoveries/topology/knot-invariants) | published |
| [A395521](https://oeis.org/A395521) | Number of distinct sandpile groups of connected graphs on n vertices | [open-problems/discoveries/algebraic-graph-theory/sandpile-groups](https://github.com/queelius/open-problems/blob/main/discoveries/algebraic-graph-theory/sandpile-groups) | published |
| [A395644](https://oeis.org/A395644) | Number of fibered prime knots with n crossings | [open-problems/discoveries/topology/knot-invariants](https://github.com/queelius/open-problems/blob/main/discoveries/topology/knot-invariants) | published |
| [A397074](https://oeis.org/A397074) | NPN classes of n-variable Boolean functions with sensitivity = n | [open-problems/discoveries/boolean-measures-n5](https://github.com/queelius/open-problems/blob/main/discoveries/boolean-measures-n5) | published |
| [A397075](https://oeis.org/A397075) | NPN classes of n-variable Boolean functions with (real) degree = n | [open-problems/discoveries/boolean-measures-n5](https://github.com/queelius/open-problems/blob/main/discoveries/boolean-measures-n5) | published |
| [A397076](https://oeis.org/A397076) | NPN classes of n-variable Boolean functions with decision-tree complexity = n | [open-problems/discoveries/boolean-measures-n5](https://github.com/queelius/open-problems/blob/main/discoveries/boolean-measures-n5) | published |
| [A396815](https://oeis.org/A396815) | Maximum size of a subset of Z/nZ partitionable into two sum-free sets | [computational-explorations/docs/oeis_candidates.md](https://github.com/queelius/computational-explorations/blob/main/docs/oeis_candidates.md) | published |
| [A396816](https://oeis.org/A396816) | Maximum size of a subset of Z/nZ partitionable into three sum-free sets | [computational-explorations/docs/oeis_candidates.md](https://github.com/queelius/computational-explorations/blob/main/docs/oeis_candidates.md) | published |

### Extended (new terms on classical sequences)

| Sequence | Description | Provenance | Status |
|----------|-------------|------------|--------|
| [A006615](https://oeis.org/A006615) | z(n,n;3,4): least k forcing an all-ones 3x4 submatrix in an n X n 0/1 matrix a(10)=67, a(11)=79 | [open-problems/src/zarankiewicz.py](https://github.com/queelius/open-problems/blob/main/src/zarankiewicz.py) | published |
| [A006622](https://oeis.org/A006622) | z(n,n+1;3,4): least k forcing an all-ones 3x4 submatrix in an n X (n+1) matrix a(9)=61, a(10)=73, a(11)=85 | [open-problems/src/zarankiewicz.py](https://github.com/queelius/open-problems/blob/main/src/zarankiewicz.py) | published |
| [A006625](https://oeis.org/A006625) | z(n,n+2;3,4): least k forcing an all-ones 3x4 submatrix in an n X (n+2) matrix a(9)=67, a(10)=79 | [open-problems/src/zarankiewicz.py](https://github.com/queelius/open-problems/blob/main/src/zarankiewicz.py) | published |

### Candidate / in review

| Sequence | Description | Provenance | Status |
|----------|-------------|------------|--------|
| [A006672](https://oeis.org/A006672) | a(n) = Ramsey number r(C_4, K_{1,n}) a(11)=16 | [open-problems/discoveries/oeis-extensions-jun2026](https://github.com/queelius/open-problems/blob/main/discoveries/oeis-extensions-jun2026) | in-review |
| [A316632](https://oeis.org/A316632) | Sprague-Grundy value for Node-Kayles on the 3 X n grid graph P_3 X P_n a(17)=2 | [open-problems/discoveries/oeis-extensions-jun2026](https://github.com/queelius/open-problems/blob/main/discoveries/oeis-extensions-jun2026) | in-review |
| [A046671](https://oeis.org/A046671) | Nim-values G(3,n) for Sylver coinage a(17)=14 (submitted), a(18..30) computed, held | [open-problems/discoveries/oeis-extensions-jun2026](https://github.com/queelius/open-problems/blob/main/discoveries/oeis-extensions-jun2026) | in-review |

### Other dispositions

| Sequence | Description | Provenance | Status |
|----------|-------------|------------|--------|
| [A042950](https://oeis.org/A042950) | Comment: Lucas-triangle row sums = Boolean-measure tuple counts | [open-problems/discoveries/boolean-measures-n5](https://github.com/queelius/open-problems/blob/main/discoveries/boolean-measures-n5) | withdrawn |
| [A069001](https://oeis.org/A069001) | 3-row Chomp Grundy values | [open-problems/src/chomp.py](https://github.com/queelius/open-problems/blob/main/src/chomp.py) | not-submittable |

<!-- AUTOGEN:contributions END -->

## Reviewers

Submissions reviewed and approved by OEIS editors including Michel Marcus,
Sean A. Irvine, Max Alekseyev, Manfred Scheucher, and Jon E. Schoenfield.

---

The OEIS has been a shared ledger of integer sequences since 1964. Adding a
term is a small, permanent contribution to a commons that outlasts any one
project: a number that was unknown is now known, checked, and citable.
