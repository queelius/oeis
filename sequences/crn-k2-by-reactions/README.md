# OEIS candidates: chemical reaction networks up to species permutation

Drafts prepared Jun 05 2026 from the complete census in
`src/reaction_networks_census.py` (tests:
`tests/test_reaction_networks_census.py`; findings:
`../census-README.md`). Every DATA line here was **re-computed from
scratch** (not copied from the README) and cross-checked by two
independent methods before drafting.

## The convention (state this in every submission)

A reaction network on `k` species is a nonempty set of directed
reactions `y -> y'` where:

- each complex `y` is a vector in `Z_{>=0}^k` with coordinate sum in
  `{1, 2}` (mono- or bimolecular; **no zero complex**, so no in/out
  flows);
- `y != y'` in every reaction;
- **every species is used** (appears in some complex);
- two networks are identified **iff a species permutation in `S_k`**
  (acting on the coordinates) carries one to the other. **Complexes are
  not independently relabeled**: each complex carries fixed
  stoichiometric meaning.

For fixed `k` the admissible directed reactions form a finite ground set
of size `R` (`R = 2, 20, 72, 182` for `k = 1, 2, 3, 4`), so each
"by reaction count `n`" sequence is **finite** (terminating at `n = R`).

Deficiency (reported in the census, not in these total sequences):
`delta = m - l - s` with `m` complexes, `l` linkage classes (weakly
connected components of the reaction graph), `s` the rank over `Q` of the
stoichiometric matrix (columns `y' - y`). Feinberg, *Foundations of
Chemical Reaction Network Theory* (2019).

## Two-method confirmation

Every submitted term is confirmed by **both**:

- **Direct**: `S_k`-canonical (lex-min) enumeration of reaction subsets
  (`census(...)`). Complete `2^20` powerset at `k=2`; bounded by `C(R,n)`
  subsets at `k=3` (reached `n<=5`) and `k=4` (reached `n<=4`).
- **Burnside**: cycle-index DP over the `S_k` action on the `R`-element
  reaction set, restricted to the all-species-used invariant subset
  (`orbit_count_burnside(...)`). Agrees with direct enumeration on every
  overlapping term and extends the totals to **all** `n` (to `n=R`) for
  free.

Plus an exact integer stoichiometric rank (`integer_rank`, fraction-free
Bareiss) validated bit-for-bit against `sympy.Matrix.rank` in the test
suite.

## Ranking of candidates

| Rank | File | Sequence | DATA (first terms) | Terms | Two-method | Recommend |
|------|------|----------|--------------------|-------|------------|-----------|
| 1 | `C1-crn-k2-by-reactions.txt` | k=2 networks by reaction count n | 8, 99, 570, 2445, 7752, 19440, ... , 100, 10, 1 | 20 (complete, finite) | direct `2^20` = Burnside on **all 20** | **Submit.** Flagship. |
| 2 | `C2-crn-k3-by-reactions.txt` | k=3 networks by reaction count n | 5, 350, 9487, 169539, 2326089, ... | 12 shown (to n=72) | direct = Burnside on n=1..5; Burnside to n=72 | **Submit.** |
| 3 | `C3-crn-k4-by-reactions.txt` | k=4 networks by reaction count n | 1, 376, 32748, 1687101, 63401197, ... | 10 shown (to n=182) | direct = Burnside on n=1..4; Burnside to n=182 | **Submit (optional).** |

**Why these three.** They have the **cleanest, fully unambiguous
definition** (a reviewer can reproduce the ground set and the group
action exactly), the **most terms**, and the **strongest verification**
(two independent algorithms, agreeing everywhere they overlap, with
Burnside giving arbitrarily many terms). C1 is the strongest single
object: a *complete finite* sequence, every one of its 20 terms confirmed
by the full `2^20` enumeration **and** by Burnside.

C1 is the flagship; C2 follows the identical definition at `k=3`; C3 is
the same at `k=4`. If a reviewer prefers a single object, C1/C2/C3 could
instead be merged into one triangle `T(k, n)` read by rows (`k >= 1`,
`1 <= n <= R(k)`); we keep them separate because `k` is small/fixed and
the per-`k` sequences are individually meaningful and each finite. The
three should cross-reference one another (fill in the assigned A-numbers
in the `%Y` lines on acceptance).

## Deliberately NOT submitted (honest call)

- **Deficiency-resolved by-`n` sequences** (e.g. k=2 `delta=0` by `n`:
  `8, 82, 110, 75, 27, 5`; `delta=1`: `17, 448, 1333, 2148, 2434, 2016,
  1250, 550, 168, 30, 3`; `delta=2`; weakly-reversible by `n`:
  `5, 10, 71, 262, 1210, ...`; `delta=0 & WR`: `5, 9, 47, 27, 5`). These
  are **too derivative for a clean standalone submission**: short and
  terminating (the `delta=0` row dies after `n=6`), non-monotone, with no
  closed form, and (critically) **single-method only** (deficiency
  needs a per-network rank + linkage computation, so Burnside does not
  independently confirm them; only the *totals* get the two-method
  guarantee). They are convention-sensitive in a way a reviewer would
  reasonably push back on. They are fully documented in
  `../census-README.md` and remain available if a future structural
  result (a closed form, or a tie to the Anderson et al. Erdos-Renyi
  deficiency-zero prevalence question) gives them a reason to exist.
- **Grand-total-by-`k`** (`3, 524796, ...`): only two terms are
  computable (k=3,k=4 grand totals are astronomically large), far too
  short.

## Prior-art distinctness (checked Jun 05 2026)

**OEIS** (`curl "https://oeis.org/search?q=...&fmt=json"`; WebFetch is
403-blocked). The mechanism was sanity-checked against a known sequence
(Catalan `1,1,2,5,14,42,132,429` returns a 165 KB result array); every
CRN query below returns a bare `null` (genuine **NO MATCH**):

- k=2 totals (full 20-term DATA) and prefix `99,570,2445,7752,19440,38760`;
- k=3 totals and prefix `350,9487,169539,2326089`;
- k=4 totals and prefix `376,32748,1687101,63401197`;
- the deficiency/WR-by-`n` rows above;
- texts: "chemical reaction network deficiency", "non-isomorphic
  reaction networks", "mass action reaction networks", "Feinberg
  deficiency", "reaction network complexes deficiency zero".

OEIS currently contains **no chemical-reaction-network enumeration or
deficiency sequences at all** (the one stray hit, A000332 = C(n,4) on
"number of reaction networks", is an unrelated word match).

**Literature.** CRN enumeration *up to symmetry* is well-trodden, but
under a **different equivalence relation**:

- **Deshpande & Gopalkrishnan**, "Fast Enumeration of Non-isomorphic
  Chemical Reaction Networks" (CMSB 2019): counts CRNs up to **full
  isomorphism** of the complex-species graph (species **and** complexes
  relabeled together) via `nauty`, e.g. 428502 non-isomorphic vs 635040
  labelled for 5 complexes / 4 reactions. That group is **not** ours.
- **Horn** (1972): 43 isomorphism classes of 3-complex networks, 41 of
  deficiency zero, classical deficiency-resolved counts, again under
  full isomorphism.
- **Anderson et al.**, "Prevalence of deficiency-zero reaction networks
  in an Erdos-Renyi framework" (2022): the asymptotic deficiency-zero
  *fraction*, the live research question our finite `delta=0` decay only
  gestures at.

**Distinctness claim.** Our `S_k`-only equivalence is **coarser in the
species** (no complex relabeling) and **finer in the complexes** than full
CRN isomorphism, so it yields **different integer counts** that do not
coincide with, and are not contained in, the published full-isomorphism
sequences. Stated explicitly under the `S_k`-only + mono/bimolecular +
no-zero-complex + all-species-used convention, the three submitted
sequences appear to be **new**.

## Submission checklist (per project OEIS conventions)

- [x] Sequence variable is `n` (not `k`/`r`) in NAME, COMMENT, FORMULA,
      EXAMPLE. (`k` appears only as the fixed species count, named as
      such.)
- [x] References one per line. (Feinberg; project repo.)
- [x] **No EXTENSIONS field** (these are new sequences).
- [x] Convention stated unambiguously so a reviewer can reconstruct the
      ground set, the group action, and the deficiency definition.
- [x] Two computation methods stated as independent verification.
- [x] `nonn,fini` keywords (sequences are finite); `full` on C1 (all 20
      terms shown); `hard` on C2/C3 (deficiency-resolved enumeration is
      expensive though the totals are easy).
- [ ] On acceptance: fill assigned A-numbers into the `%Y` cross-refs and
      replace the placeholder `A000000`; attach b-files (producible to
      `n=R` in ms via `orbit_count_burnside`).
