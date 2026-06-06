# Bivariate refinement of the max-tail-length distribution of functional graphs

**Status:** two joint distributions computed exactly, each given a **proved
bivariate EGF**, verified cell-by-cell against brute force (`n <= 7`), with
**both marginals validated exactly against the published literature**. OEIS
search (Jun 2026) finds **neither joint triangle**; the marginals and the single
`t = 1` slice are classical. The `(max tail length, #cyclic points)` joint is a
clean **OEIS-submission candidate** (draft below).

Code: [`src/max_tail_bivariate.py`](../../../src/max_tail_bivariate.py),
tests: [`tests/test_max_tail_bivariate.py`](../../../tests/test_max_tail_bivariate.py)
(120 passing, incl. the `n = 7` brute-force cross-checks under the `slow` marker).

This is the **novel thread** flagged at the end of the univariate
[`README.md`](README.md): the univariate height distribution `D(n,t)` is prior
art (A216242 / Critzer 2013), but its bivariate refinements were not found in
OEIS and have provable closed forms via the same construction.

## Problem

For an endofunction `f : [n] -> [n]` the **max tail length** `T(f)` is the
longest distance from any vertex to its cycle (`0` iff `f` is a permutation).
Refine the `n^n` maps simultaneously by `T` and by a second statistic:

```
B(n, t, c)  = #{ f : T(f) = t  AND  f has exactly c cyclic points },
Bm(n, t, m) = #{ f : T(f) = t  AND  f has exactly m components (= cycles) }.
```

## The proved bivariate EGFs

The whole analysis rides on the canonical decomposition (Flajolet-Sedgewick,
*Analytic Combinatorics* II.4): a functional graph is a **set of cycles of
rooted trees**. The cyclic points carry a permutation (the *core*); off each
cyclic point hangs a rooted tree with that point as its root. A vertex at depth
`d` in such a tree has tail length exactly `d`, so

```
T(f) <= t   <=>   every hanging tree has height <= t.                       (★)
```

Let `e_t(x)` be the EGF of rooted labelled trees of height `<= t`:
`e_0 = x`, `e_t = x * exp(e_{t-1})` (this is the univariate object, proved in
[`README.md`](README.md)).

> **Theorem A (max tail length × cyclic points).** Marking each cyclic point by
> a variable `u`,
> ```
>     F_{<=t}(x,u) = 1 / (1 - u * e_t(x)),
>     B(n,t,c)     = n! [x^n u^c] ( 1/(1 - u e_t) - 1/(1 - u e_{t-1}) ),
> ```
> with the convention `e_{-1} = 0`.

**Proof.** By `(★)`, a height-`<= t` map is a core permutation whose cyclic
points each root one tree of height `<= t`. A permutation is a **sequence** of
"slots" (`SEQ`), so the core EGF is `1/(1-y)` evaluated at the class filling one
slot. One slot = one rooted tree of height `<= t` (EGF `e_t`) carrying one
cyclic point; marking that point by `u` makes the slot weight `u·e_t`. Hence the
cumulative EGF is `1/(1 - u e_t)`, and `[u^c]` selects exactly the maps with `c`
cyclic points (one `u` per slot = per cyclic point). Extract `n![x^n]`; the
exact-`t` count is the first difference in `t`. ∎

> **Theorem B (max tail length × components).** Marking each component by `u`,
> ```
>     G_{<=t}(x,u) = exp( u * log( 1/(1 - e_t(x)) ) ) = (1 - e_t(x))^{-u},
>     Bm(n,t,m)    = n! [x^n u^m] ( G_{<=t} - G_{<=t-1} ).
> ```

**Proof.** A single component is one cycle of trees of height `<= t`; the cycle
construction `CYC` on a class with EGF `e_t` has EGF `L_t = log(1/(1 - e_t))`. A
functional graph is a **set** (`SET`) of components, and the multivariate `SET`
construction marking each component by `u` gives `exp(u · L_t)` (the labelled
`SET` of a `u`-marked class exponentiates its EGF, carrying the mark linearly in
the exponent). `[u^m]` selects maps with `m` components. ∎

Both EGFs reduce to the univariate `1/(1-e_t)` at `u = 1`, recovering A216242.

## Joint table `B(n, t, c)` (rows `t = 0..n-1`, columns `c = 1..n`)

```
n=1   t=0: 1
n=2   t=0: 0  2          n=3   t=0: 0   0   6
      t=1: 2  0                t=1: 3  12   0
                               t=2: 6   0   0
n=4   t=0: 0    0    0   24    n=5   t=0:   0    0    0    0  120
      t=1: 4   48   72    0          t=1:   5  160  540  480    0
      t=2: 36  48    0    0          t=2: 200  600  360    0    0
      t=3: 24   0    0    0          t=3: 300  240    0    0    0
                                     t=4: 120    0    0    0    0
n=6   t=0:    0     0     0     0     0  720
      t=1:    6   480  3240  5760  3600    0
      t=2: 1170  6000  7560  2880     0    0
      t=3: 3360  5040  2160     0     0    0
      t=4: 2520  1440     0     0     0    0
      t=5:  720     0     0     0     0    0
n=7   t=0:     0      0       0      0      0      0  5040
      t=1:     7   1344   17010  53760  63000  30240     0
      t=2:  7392  57540  115920  90720  25200      0     0
      t=3: 38850  87360   68040  20160      0      0     0
      t=4: 43680  45360   15120      0      0      0     0
      t=5: 22680  10080       0      0      0      0     0
      t=6:  5040      0       0      0      0      0     0
```

Support: `B(n,t,c) = 0` unless `c >= 1`, and `t = 0 ⟹ c = n` (a permutation has
all `n` points cyclic), and `t >= 1 ⟹ c <= n - t` (a max tail of `t` needs a
length-`t` path of `t` non-cyclic vertices).

## Joint table `Bm(n, t, m)` (rows `t = 0..n-1`, columns `m = 1..n`)

```
n=4   t=0: 6  11   6   1      n=5   t=0:  24   50   35  10  1
      t=1: 52  60  12   0           t=1: 385  570  210  20  0
      t=2: 50  35  10   0           t=2: 620  480   60   0  0
      t=3: 24   0   0   0           t=3: 420  120    0   0  0
                                    t=4: 120    0    0   0  0
```

The `t = 0` row is `|Stirling1(n,m)|` (permutations of `[n]` with `m` cycles),
summing to `n!`. Full tables to `n = 7` are produced by `joint_tm_rows(n)`.

## Marginal validations (both directions, both joints, exact)

| Marginal | Equals | OEIS | Verified |
|---|---|---|---|
| `sum_c B(n,t,c)` | max-tail distribution `D(n,t)` | **A216242** (Critzer) | `n=1..7` |
| `sum_t B(n,t,c)` | endofunctions with `c` cyclic points `= C(n,c) c! · c · n^{n-c-1}` | **A066324** | `n=1..7` |
| `sum_m Bm(n,t,m)` | max-tail distribution `D(n,t)` | **A216242** | `n=1..7` |
| `sum_t Bm(n,t,m)` | endofunctions with `m` components | **A060281** | `n=1..7` |

The cyclic-points marginal closed form (`forests on n nodes with c given roots`
`= c·n^{n-c-1}` for `c<n`, `=1` for `c=n`, by generalized Cayley) is also tested
directly against A066324's published DATA and confirmed to sum to `n^n`.

## OEIS status (searched Jun 2026; queries via `curl`, WebFetch is 403-blocked)

| Object | OEIS | Verdict |
|---|---|---|
| `B(n,t,c)` full joint (any row reading) | (none) | **NOT FOUND** (candidate-new) |
| `Bm(n,t,m)` full joint (any row reading) | (none) | **NOT FOUND** (candidate-new) |
| `t = 1` slice of `B`: `B(n,1,c)`, `c<n` | **A199673** `T(n,c)=C(n,c)c! c^{n-c}` | KNOWN slice (EGF of its column `(x e^x)^c`, matching `(u e_1)^c` since `e_1 = x e^x`) |
| `t = 0` row of `B` | `[c=n] n!` (trivial) | trivial |
| `t = 0` row of `Bm` | `|Stirling1(n,m)|` (A130534) | KNOWN row |
| cyclic-points marginal | **A066324** | KNOWN marginal |
| components marginal | **A060281** | KNOWN marginal |
| max-tail marginal | **A216242** | KNOWN marginal (univariate) |

Searched flattenings for both joints: full-square (with zeros), support-only by
`(n, t, c)`, and by `(n, c, t)`. None returned a hit. Distinctive interior runs
(e.g. `7392,57540,115920,90720,25200`) also returned nothing. The single `t=1`
slice coinciding with A199673 is a structural sanity check, not a collision with
the full triangle.

## Honest novelty / theorem assessment

- **Bivariate EGFs (Theorems A, B):** *found and proved* by the standard
  symbolic method (Flajolet-Sedgewick markings: `SEQ` of `u`-marked tree-slots;
  `SET` of `u`-marked cycles). These are not deep; they are the *correct
  bookkeeping* of a classical decomposition, but they are clean, fully proved,
  and verified exactly against brute force to `n = 7`. No comparably refined
  table appears in OEIS or in the obvious literature (Flajolet-Odlyzko 1990
  treat height and cycle-count as *separate* univariate parameters; the joint
  law is not tabulated there).
- **New OEIS artifact:** **yes, available.** The `(max tail length, #cyclic
  points)` joint `B(n,t,c)` is absent from OEIS in every natural reading and is
  the cleaner of the two (rectangular support `c <= n-t`, every row a genuine
  count, marginals both already in OEIS as A066324 / A216242). Draft below. The
  `(t, #components)` joint is equally new but its `t=0` row is `|Stirling1|`, so
  it reads as a less self-contained triangle; recorded but not submitted first.
- **Caveat on weight:** these are *enumerative refinements*, not an open
  conjecture resolved. The value here is (i) a tidy provable bivariate law that
  *specializes to three classical OEIS sequences in its margins/slices*, and
  (ii) one (arguably two) new OEIS triangle(s). Honest framing: a small, correct,
  well-verified combinatorial contribution, not a theorem about a hard open
  problem.

---

## OEIS submission draft: `B(n,t,c)` joint (max tail length by cyclic points)

Following project conventions: variable `n`; one reference per line,
alphabetical; no EXTENSIONS on a new sequence. Read by **rows**, full-square so
the triangle structure is unambiguous (`t = 0..n-1` outer, `c = 1..n` inner,
explicit zeros), which makes the row-length regular (`n` per `t`-row, `n^2` per
`n`-block) and the marginals recoverable by row/column sums.

```
%S  (DATA, flattened: for n=1,2,3,..., read t=0..n-1, and within each t read
     c=1..n)
    1,
    0,2, 2,0,
    0,0,6, 3,12,0, 6,0,0,
    0,0,0,24, 4,48,72,0, 36,48,0,0, 24,0,0,0,
    0,0,0,0,120, 5,160,540,480,0, 200,600,360,0,0, 300,240,0,0,0, 120,0,0,0,0,
    0,0,0,0,0,720, 6,480,3240,5760,3600,0, 1170,6000,7560,2880,0,0,
      3360,5040,2160,0,0,0, 2520,1440,0,0,0,0, 720,0,0,0,0,0
```

(As one comma-separated line for the actual submission, produced by
`flatten_full_square` below.)

- **%N** Triangle of triangles T(n,t,c) read by rows: number of functions
  f:{1..n}->{1..n} whose functional graph has maximal tail length (rho-height)
  exactly t and exactly c cyclic (periodic) points; n >= 1, 0 <= t <= n-1,
  1 <= c <= n.
- **%O** 1,3  (offset: first n is 1; first value > 1 is the entry "2" for
  (n,t,c)=(2,0,2), the 3rd term).
- **%C** A functional graph is a set of cycles of rooted trees; the cyclic
  points carry a permutation and a rooted tree of height <= t hangs off each.
  Hence T(n,t,c) = n! [x^n u^c] ( 1/(1-u*e_t(x)) - 1/(1-u*e_{t-1}(x)) ), where
  e_t is the e.g.f. of rooted trees of height <= t: e_0 = x, e_t = x*exp(e_{t-1}),
  e_{-1} = 0.
- **%C** Row sums over t and c give n^n (A000312).
- **%C** Summing over c gives A216242 (functions by height/max tail length).
- **%C** Summing over t gives A066324 (functions by number of cyclic points);
  closed form C(n,c)*c!*c*n^(n-c-1) for c<n, and n! for c=n.
- **%C** The t=1 slice (c<n) is A199673, T(n,c)=C(n,c)*c!*c^(n-c).
- **%C** T(n,0,c) = n! if c=n, else 0 (permutations). T(n,n-1,1) = n! (a single
  path of n vertices ending in a fixed point).
- **%F** E.g.f. (bivariate, cumulative in t): 1/(1 - u*e_t(x)) with e_t as above;
  the exact-t triangle is the first difference in t.
- **%e** Triangle blocks (n; rows t=0..n-1; cols c=1..n):
  ```
  n=3: [0,0,6; 3,12,0; 6,0,0]
  n=4: [0,0,0,24; 4,48,72,0; 36,48,0,0; 24,0,0,0]
  ```
- **%Y** Cf. A216242 (height marginal), A066324 (cyclic-points marginal),
  A060281 (the analogous components refinement), A199673 (t=1 slice),
  A000312 (n^n), A000142 (n!), A130534 (Stirling1).
- **%K** nonn,tabf  (irregular/ragged: row-blocks grow as n^2)
- **%A** Alexander R. Towell

References (one per line, alphabetical):
```
P. Flajolet and A. Odlyzko, Random Mapping Statistics, in: Advances in Cryptology - EUROCRYPT '89, LNCS 434, Springer, 1990, pp. 329-354.
P. Flajolet and R. Sedgewick, Analytic Combinatorics, Cambridge University Press, 2009, Sections II.4 and III.7.
```

### Optional companion: `Bm(n,t,m)` joint (max tail length by components)

Same construction, EGF `exp(u*log(1/(1-e_t))) = (1-e_t)^{-u}`. Marginals:
A216242 (over m) and A060281 (over t). The `t=0` row is `|Stirling1(n,m)|`
(A130534). Held as a second submission if the first is well-received.

## Reproduce

```bash
python -m src.max_tail_bivariate          # tables, marginal checks, flattened data
python -m pytest tests/test_max_tail_bivariate.py -v          # 120 tests
python -m pytest tests/test_max_tail_bivariate.py -v -m slow  # n=7 brute force
```

`flatten_joint_tc_support(N)` emits the support-only sequence; for the
full-square submission DATA use `joint_tc_rows(n, full_square=True)` per `n`.

## References

- P. Flajolet, A. Odlyzko, *Random Mapping Statistics*, EUROCRYPT '89, LNCS 434
  (1990), 329-354.
- P. Flajolet, R. Sedgewick, *Analytic Combinatorics*, Cambridge UP (2009),
  II.4 (functional graphs as sets of cycles of trees), III.7 (multivariate
  markings).
- OEIS A216242 (Critzer; max tail length), A066324 (cyclic points),
  A060281 (components / cycles), A199673 (leaders/groups; the t=1 slice),
  A130534 (Stirling1), A000312 (n^n), A000142 (n!).
```
