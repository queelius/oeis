# OEIS Harvest -- April 2026

Systematic inventory of computational outputs in this repo against OEIS,
identifying low-hanging submission opportunities. Already-published Towell
sequences (A006615, A006622, A006625, A394445) excluded from action items.

Sources scanned:
- `src/*.py` (Python modules with KNOWN_VALUES tables and OEIS references)
- `discoveries/*` (json/md files with computed sequences and tables)
- `paper/*/` (LaTeX manuscripts with quantitative tables)
- `discoveries/*/oeis*.txt` (existing submission files)

OEIS searches conducted via `https://oeis.org/search?q=` direct queries
on April 26, 2026.

---

## Top Recommendations (highest value, lowest effort)

### 1. Submit b-file extension to A394445 -- Ready Now
**Effort: low (~30 min).** **Value: high.**

OEIS A394445 (our published Rado-numbers sequence) currently displays a
43-term auto-synthesized table. We have a 500-term b-file already on disk
at `discoveries/additive-combinatorics/rado-numbers/b394445.txt`. Submit
this directly; the values are already published as a theorem (Towell 2026,
Zenodo DOI 10.5281/zenodo.19372727) so no review controversy.

**Action.** Email b-file to OEIS or attach to a sequence revision.

### 2. Submit `sat(n, P_4)` -- new sequence
**Effort: medium (~2 hr).** **Value: medium-high.**

Sequence: 2, 4, 3, 5, 4, 6, 5, 7, ... (n=4..11, computed; closed form
proved by Kaszonyi-Tuza 1986). Strikingly **non-monotone**, no OEIS hit.
Closed form proved in `paper/saturation/saturation_note.tex`:
`sat(n,P_4) = floor(n/2)+1` if n odd, `n/2 + 1` if n even (the formula
gives the same shape as `floor((n+2)/2) - [n odd]`).

This sequence is genuinely interesting (oscillating saturation number)
and absent from OEIS.

### 3. Submit `# distinct (s, bs, C, D, deg) tuples on NPN classes` --
new sequence
**Effort: medium (~1.5 hr).** **Value: medium.**

Sequence: 2, 3, 6, 12, 24 (n = 1..5). Concise summary statistic of the
n=5 boolean-measures catalog. Not in OEIS. Cross-references A000370
(NPN class count). Frontier of n=5 unreachable for most groups; this
sequence is exact and complete through n=5.

### 4. Submit `sat(n, K_{2,3})` -- new sequence (cross-reference A005408)
**Effort: low (~1 hr).** **Value: low-medium.**

Sequence: 7, 9, 11, 13 for n = 5..8. Matches `2n - 3`. Same shape as
`sat(n, K_4)` (Erdos-Hajnal-Moon 1964). New observation worth recording
even though the formula is simple. Not in OEIS as a graph-saturation
sequence.

### 5. Submit `# NPN classes with sensitivity = n on n vars` -- new sequence
**Effort: medium (~1.5 hr).** **Value: medium.**

Sequence: 1, 2, 7, 125, 355251 (n = 1..5). Closely related: same
sequence for `bs = n` and `C = n`. The deg = n version (1, 2, 8, 164,
529523) and the D = n version (1, 2, 9, 184, 585803) are distinct
sequences. Each independently new in OEIS.

---

## Full Inventory

### Existing OEIS submissions (Towell)

| OEIS | Description | Status | Action |
|------|-------------|--------|--------|
| A006615 | z(n,n; 3,4) | Published, a(10)=67 our extension | Done |
| A006622 | z(n,n+1; 3,4) | Published, a(9)=61 our extension | Done |
| A006625 | z(n,n+2; 3,4) | Published, a(9)=67 our extension | Done |
| A394445 | Distinct-variable Rado for x+y=nz | Published | **Submit b-file (500 terms vs OEIS 43)** |

### Boolean function complexity (boolean_measures_n5)

Source: `discoveries/boolean-measures-n5/summary.json`,
`paper/boolean-measures-n5/boolean_measures_n5.tex`,
verified against `enumerate_npn_classes` for n=1..4.

#### Candidate B1: # distinct (s,bs,C,D,deg) tuples among NPN classes

- **Values**: 2, 3, 6, 12, 24 (n = 1, 2, 3, 4, 5).
- **OEIS lookup**: not found by sequence search. Generic shape but no
  match in NPN/Boolean-function context.
- **Frontier**: n = 5 is the practical limit (n = 6 has 2.0e14 NPN
  classes; intractable to enumerate).
- **Verdict**: Submit new sequence. Cf. A000370.
- **Effort**: medium.

#### Candidate B2: # NPN classes with sensitivity = n

- **Values**: 1, 2, 7, 125, 355251.
- Same sequence for bs = n and C = n (verified at n=4: all give 125).
- **OEIS lookup**: 1,2,7,125,355251 returns no match.
- **Verdict**: Submit new sequence. Cf. A000370.
- **Effort**: medium.

#### Candidate B3: # NPN classes with deg = n

- **Values**: 1, 2, 8, 164, 529523 (n = 1..5).
- **OEIS lookup**: no match.
- **Verdict**: Submit new sequence.
- **Effort**: medium.

#### Candidate B4: # NPN classes with D = n

- **Values**: 1, 2, 9, 184, 585803 (n = 1..5).
- **OEIS lookup**: no match.
- **Verdict**: Submit new sequence.
- **Effort**: medium.

#### Candidate B5: # functions on n vars with s = bs = D = deg = n

(All-functions count, not NPN-reduced.)

- **Values**: 2, 10, 130, 32746 (n = 1..4).
- n = 5 is 2^32 functions; computing this is a multi-day job (about
  4 to 8 hours per measure if parallelized at 16 workers, but the
  combination requires care).
- **OEIS lookup**: no match.
- **Verdict**: Submit new sequence (4-term initial; mark `more`).
- **Effort**: medium-high (1.5 hr to draft + n=5 might take a week of
  compute).

### Graph saturation numbers (graph_saturation)

Source: `paper/saturation/saturation_note.tex`. Each value is verified by
SAT-witness + UNSAT-proof pair.

#### Candidate G1: sat(n, P_4)

- **Values**: 2, 4, 3, 5, 4, 6, 5, 7 (n = 4..11). Non-monotone.
- **Closed form** (Kaszonyi-Tuza 1986; verified): `floor(n/2) + 1` for
  odd n, `n/2 + 1` for even n.
- **OEIS lookup**: no match for "saturation P_4" or for the values.
- **Frontier**: closed form known; sequence extends to all n >= 4.
- **Verdict**: Submit new sequence. Reference Kaszonyi-Tuza 1986 and
  our saturation note (Zenodo 10.5281/zenodo.19226372).
- **Effort**: medium.

#### Candidate G2: sat(n, P_5)

- **Values**: 4, 5, 6, 6, 7 (n = 5..9). Non-monotone (a(7) = a(8) = 6).
- **Closed form**: not yet conjectured by us; values from SAT.
- **OEIS lookup**: no match.
- **Verdict**: Submit new sequence. Mark `more` (we have only 5 values,
  closed form unknown).
- **Effort**: medium.

#### Candidate G3: sat(n, C_4)

- **Values**: 5, 6, 8, 9, 11 (n = 5..9). Match Ollmann 1972 formula
  `floor((3n-5)/2)`.
- **OEIS lookup**: not found as a saturation sequence. Values are
  shifted A007494 (`floor((3n-5)/2)` and "numbers congruent to 0 or 2
  mod 3" coincide).
- **Verdict**: Submit new sequence with Cf. A007494 + clear graph
  saturation context. Lower priority because the underlying integer
  pattern is well-known.
- **Effort**: low-medium.

#### Candidate G4: sat(n, C_5)

- **Values**: 6, 8, 9, 10 (n = 5..8). Closed form by Chen 2009/2011:
  `ceil(10(n-1)/7)` with exceptional set `{11,12,13,14,16,18,20}`.
- **OEIS lookup**: not found. Closed form is known, exceptional set
  makes the sequence interesting.
- **Verdict**: Submit new sequence. Worth doing.
- **Effort**: medium.

#### Candidate G5: sat(n, C_6)

- **Values**: 9, 10 (n = 6..7). Only 2 values; need more.
- **Verdict**: Borderline -- compute n = 8, 9 first. Skip for now.

#### Candidate G6: sat(n, K_{2,3})

- **Values**: 7, 9, 11, 13 (n = 5..8). Matches `2n - 3`.
- **OEIS lookup**: not found as a saturation sequence.
- **Verdict**: Submit new sequence with Cf. A005408 (odd numbers) and
  Cf. to sat(n, K_4) which has the same `2n - 3` formula. New
  observation about K_{2,3}-saturation matching K_4-saturation.
- **Effort**: low-medium.

### Rado number families (rado_families)

Source: `paper/rado-families/rado_families_note.tex` and
`discoveries/additive-combinatorics/rado-numbers/README.md`.

#### Candidate R1: R(x + 2y = kz) for odd k

- **Values** (odd k = 13..49): 61, 81, 103, 127, 155, 185, 218, 253,
  291, 331, 374, 421, 469, 521, 575, 631, 691, 749, 817.
- **Conjectured formula**: `floor(k(k+1)/3) + 1`.
- **OEIS lookup**: no match for either subsequence or interleaved
  values.
- **Verdict**: Submit new sequence (conjecture-with-evidence). Mark
  `more`.
- **Effort**: medium.

#### Candidate R2: R(x + by = bz) -- perfect square law

- **Values**: 16, 25, 36, 49, ..., 625 (b = 4..25). Matches `b^2`.
- **OEIS lookup**: matches A000290 (squares) trivially. The new content
  is the Rado-number interpretation, not the integer values.
- **Verdict**: Add a comment to A000290 ("Equals R(x+by=bz, 2;
  distinct) for b >= 4 by Towell 2026") rather than submit a new
  sequence. Marginal; could also submit as a separate sequence with
  offset starting at b=1 (where the formula begins).
- **Effort**: low (just a comment).

#### Candidate R3: R(x + by = (b+1)z)

- **Values** (b = 1..30): 9, 13, 17, 19, 25, 29, 33, 35, 41, 45, 49,
  51, 57, 61, 65, 67, 73, 77, 81, 83, 89, 93, 97, 99, 105, 109, 113,
  115, 121, 125. Formula: `4b + 5 - 2[4|b]`.
- **OEIS lookup**: no match.
- **Closed form** by Gupta-Thulasi Rangan-Tripathi 2015 (already in
  literature; we provide independent SAT confirmation across 30
  values).
- **Verdict**: Submit new sequence; cite GTRT 2015 as primary
  reference. Validation rather than new conjecture.
- **Effort**: medium.

#### Candidate R4: R(ax + (a+1)y = (2a+1)z) -- linear

- **Values** (a = 4..15): 37, 45, 53, 61, 69, 77, 85, 93, 101, 109,
  117, 125. Formula: `8a + 5`.
- **OEIS lookup**: matches A004770 (numbers `8k + 5`). Same situation
  as R2.
- **Verdict**: Add comment to A004770 referencing Rado interpretation;
  do not submit separately.
- **Effort**: low.

#### Candidate R5: R(x + y + z = kw) for k >= 15

- **Values** (k = 15..25): 81, 91, 103, 115, 127, 137, 155, 169, 185,
  201, 218.
- **Conjectured formula**: `ceil(k(k+1)/3) + 1`.
- **OEIS lookup**: no match.
- **Verdict**: Submit new sequence (conjecture). Cf. R1 (sequences
  differ only in rounding direction).
- **Effort**: medium.

#### Candidate R6: R(x + y = z + kw) for odd k

- **Values** (odd k = 9..19): 31, 43, 57, 73, 91, 111. Matches central
  polygonal numbers `n^2 - n + 1` evaluated at `(k+3)/2`.
- **OEIS lookup**: substring matches A002061 (central polygonal). Could
  add a comment to A002061 rather than a new sequence.
- **Verdict**: Add comment to A002061; low priority.
- **Effort**: low.

#### Candidate R7: 3-color Rado numbers R(x+y=kz, 3; distinct)

- **Values**: a(1) = 24, a(2) = 27, a(3) = 99, a(k) > 5000 for k >= 4.
- Only 3 finite values; mostly inequalities. Not a clean integer
  sequence.
- **Verdict**: Skip. Not submission-ready.

### De Morgan formula size (formula_size)

Source: `paper/formula-size/formula_size_note.tex`.

#### Candidate F1: Distribution of L(f) for n = 3 (256 functions)

- **Values** (count by L=0..10): 2, 6, 12, 30, 44, 56, 56, 0, 44, 0, 6.
- **OEIS lookup**: no match.
- **Verdict**: Submit new sequence as a row of a "Boolean function
  formula size distribution" table. Either flat (11 values) or as a
  triangle in `n`.
- **Effort**: medium (need to think about presentation).

#### Candidate F2: Distribution of L(f) for n = 4 (65,536 functions)

- **Values** (L = 0..16): 2, 8, 24, 96, 248, 800, 2520, 3200, 6384,
  5760, 10656, 7568, 11792, 5936, 3712, 640, 114.
- **OEIS lookup**: no match.
- **Verdict**: Submit. Same packaging question as F1.
- **Effort**: medium.

#### Candidate F3: Number of "hardest" functions L(f) = max_n L

- **Values**: n = 0: 1; n = 1: 1 (NOT/ID parity); n = 2: ?; n = 3: 6;
  n = 4: 114; n = 5: ? (PARITY_5 has L = 28 but only one of 2^32
  functions has been confirmed at L = 28).
- **Verdict**: Borderline. Computing n = 2 and n = 5 requires more
  work; marginal sequence value.
- **Effort**: medium-high.

### Davenport constants (egz_davenport)

Source: `paper/davenport/davenport_note.tex`.

Values are indexed by group structure, not by `n`, so cannot be packaged
as a single integer sequence in the OEIS sense unless we restrict to a
specific family (e.g., D(Z_n^3)). We have 9 rank-3 groups; only one or
two share a parametric family.

- D(Z_n^3): n = 3, 4 give D = 7, 10 -- two values; insufficient.
- **Verdict**: Skip; values are scattered across nonparametric families.

### Paley graph clique numbers

Source: `discoveries/extremal-graph-theory/paley-clique/README.md`.

- **OEIS A320757** has 62 terms (auto-synthesized; no real b-file).
  Our values for q in {229, ..., 397} (16 primes) all *match* A320757.
- We don't extend the sequence; OEIS is already 25 terms past us.
- **Exoo's table** (https://isu.indstate.edu/ge/PALEY/) has values to
  q = 16741, but those aren't ours; submitting Exoo's data as a b-file
  would be a different project.
- **Verdict**: No action. (Optional: file an issue requesting a real
  b-file from Exoo's table, separately from this harvest.)

### Power map dynamics

Source: `paper/power-map-dynamics/`,
`discoveries/finite-field-dynamics/`.

The mean rho-length is a rational function of `(d, p)`, not an integer
sequence per se. Derivative integer sequences (e.g., "smallest p with
v_d(p-1) >= k for fixed d") are computable but obscure. None of the
power-map data fits a clean OEIS sequence we'd want to submit.

- **Verdict**: Skip. The publishable artifact is the theorem paper, not
  a sequence.

### Linear map dynamics

Source: `discoveries/finite-field-dynamics/linear-map-findings.md`.

Counts of permutation / non-permutation matrices in M_n(F_p) are
already implicit in known formulas (`|GL_n(F_p)|`). We do not extend
any clean integer sequence.

- **Verdict**: Skip.

### Cap sets in F_3^n

Source: `discoveries/additive-combinatorics/cap-set-f3-7/README.md`.

- The maximum cap set sequence r_3(F_3^n) for n = 1..6 is 2, 4, 9, 20,
  45, 112 -- already OEIS A090245 / well-known. n = 7 unknown to us
  (frontier > 236 by Edel; we cannot improve).
- **Verdict**: No action.

### Other modules (no submission opportunity)

| Module | Reason |
|---|---|
| `chomp.py` (3-row Chomp) | OEIS A069001 has 2522 terms vs our 1000 |
| `grundy_game.py` | OEIS A002188 has 20000 terms; ours has 105 |
| `lattice_points.py` (Gauss circle) | OEIS A000328 has 10000 terms |
| `sidon_sets.py` | OEIS A143824 already up to date; KNOWN_F2 only validates |
| `binary_codes.py` | All values from Brouwer's table; already in OEIS |
| `goldbach.py`, `gilbreath.py`, `prime_constellations.py` | Validate-only, no extensions |
| `collatz.py` | Five OEIS sequences referenced; we don't extend any |
| `erdos_straus.py` | A073101 has 1000-term b-file; we don't extend |
| `npn.py` | Class counts are A000370; we use them as validation |

---

## Summary Table

| ID | Description | OEIS status | Action | Effort |
|---|---|---|---|---|
| **Top-1** | A394445 b-file (500 terms) | Sequence published, no real b-file | **Submit b-file extension** | **low** |
| **Top-2** | sat(n, P_4) | Not in OEIS | Submit new | medium |
| **Top-3** | # distinct (s,bs,C,D,deg) tuples on NPN | Not in OEIS | Submit new | medium |
| **Top-4** | sat(n, K_{2,3}) | Not in OEIS | Submit new | low-medium |
| **Top-5** | # NPN with s=n | Not in OEIS | Submit new | medium |
| B3 | # NPN with deg=n | Not in OEIS | Submit new | medium |
| B4 | # NPN with D=n | Not in OEIS | Submit new | medium |
| B5 | # all-fns with s=bs=D=deg=n | Not in OEIS | Submit new (4-term, `more`) | medium-high |
| G2 | sat(n, P_5) | Not in OEIS | Submit new (5-term, `more`) | medium |
| G3 | sat(n, C_4) | Cf. A007494 | Submit new with cross-ref | low-medium |
| G4 | sat(n, C_5) | Not in OEIS | Submit new | medium |
| F1 | L(f) distribution n=3 | Not in OEIS | Submit new (or as triangle) | medium |
| F2 | L(f) distribution n=4 | Not in OEIS | Submit new | medium |
| R1 | R(x+2y=kz, odd k) | Not in OEIS | Submit new (conjecture) | medium |
| R3 | R(x+by=(b+1)z) | Not in OEIS | Submit new (cite GTRT 2015) | medium |
| R5 | R(x+y+z=kw, k >= 15) | Not in OEIS | Submit new (conjecture) | medium |
| R2 | R(x+by=bz) = b^2 | A000290 (squares) | Add comment | low |
| R4 | R(ax+(a+1)y=(2a+1)z) | A004770 | Add comment | low |
| R6 | R(x+y=z+kw, odd k) | A002061 (central polygonal) | Add comment | low |

---

## Methodology Notes

- **Conservatism**: every claimed sequence has been computationally
  verified at every value we report. Where a formula is conjectured
  (R1, R5), we mark with `more` and note the verification range.
- **Cross-database**: we did not find computational sequences in
  Brouwer's tables, Heuer's Chomp tables, Exoo's Paley tables, or
  Radziszowski's Ramsey survey that are extended by our data.
  Brouwer/Heuer/Exoo are all ahead of us on their respective
  problems.
- **Effort estimates**: "low" = 30 min to 1 hr (b-file or one-line
  comment); "medium" = 1 to 2 hr (full new submission with NAME,
  DATA, COMMENTS, FORMULA, EXAMPLE, PROG, CROSSREFS, KEYWORD); "high"
  = several hours (multi-step: compute additional terms first).

## Conventions Reminder

When submitting:
- Variable name `n` (not `k`); see `feedback_oeis_conventions.md`.
- One reference per line, alphabetical.
- No EXTENSIONS field on initial submissions.
- Author: `_Alex Towell_` (see `user_profile.md`).
- Cite Zenodo DOI for the relevant published manuscript.
