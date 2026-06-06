# OEIS Extension Candidates: Discovery Survey

A curated list of OEIS sequences with `keyword:more` and/or `keyword:hard` that are amenable to computational extension via SAT solvers, exhaustive search, or computer algebra. Excludes sequences already targeted by the existing portfolio (Rado, Zarankiewicz, saturation, Davenport, cap sets in F_3^n, NPN classes, formula sizes, Paley cliques, knot invariants, sandpiles, power maps, Sidon, sunflower-free, Steiner gaps, Goldbach, Erdős-Straus, Gauss circle, Boolean measures, Chomp, Grundy).

All A-numbers verified against `oeis.org` text-format pages.

---

## Top 5 to Attack First

1. **A179183 / A230380 (Optimal binary codes for edit distance).** 8 known terms each, well-defined SAT/CP encoding (binary words with pairwise edit distance >= d), no symmetry-breaking nightmare. Plausible 1-2 new terms.
2. **A005864 / A005865 / A005866 (A(n,d) coding-theoretic functions).** Classical extremal coding theory, table of bounds at codetables.de. Tight bounds within reach via clique on Hamming-distance graphs plus LP relaxation.
3. **A287644 / A287648 (Maximum (diagonal) transversals in diagonal Latin squares).** 9 known terms, n=10 attackable with constraint propagation; this is the SAT-friendly subproblem in a deeply studied area.
4. **A000983 / A029866 (Minimal binary covering codes for radius 1, 2).** 9 / 8 known terms, classic covering code problem; ILP and structured search attack the next term.
5. **A335203 (Packing chromatic number of n-hypercube).** 8 known terms; recent SAT attacks have pushed similar problems forward; attractive because Q_n has rich automorphism group for symmetry breaking.

---

## Combinatorial Designs

| A-number | Description | Known | Next | Why open |
|---|---|---|---|---|
| **A030129** | Nonisomorphic Steiner triple systems S(2,3,n) | nonzero at n=7,9,13,15,19,21 (known: 1,1,2,80, 11084874829, 14796207517873771) | n=25 | n=21 took years; n=25 currently infeasible. **LOW**. |
| **A124119** | Nonisomorphic Steiner quadruple systems S(3,4,v) | 4 terms (v=8,10,14,16 -> 1,1,4,1054163) | v=20 | Massive jump expected; v=20 is frontier. **LOW**. |
| **A296086** | Binary self-dual codes of length 2n with highest min distance | 18 terms; sporadic gaps | parameters where best-known bound not tight | Classification problem; canonical-form filtering. **MEDIUM**. |
| **A007299** | Hadamard matrices of order 4n | 9 terms | n=9 (order 36) | Frontier; n=9 ongoing. **LOW** for full count, **MEDIUM** for inequivalence-class enumeration. |
| **A001230** | Closed knight's tours on 2nx2n | 4 terms (zeros + 9862 + 13267364410532) | 2n=10 | n=5 known to be ~10^16; very large. **LOW**. |

## Latin Squares

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A001438** | Maximum MOLS of order n | 8 terms | n=10 (conjectured 2; proved >=2, <=6) | Famous open: a(10). Established frontier, no quick win. **LOW** for definitive, **MEDIUM** for tightening upper bound. |
| **A287644** | Max transversals in diagonal Latin square of order n | 9 terms | n=10 | SAT/CP with diagonal plus transversal constraints. Highly tractable structure. **HIGH**. |
| **A287648** | Max diagonal transversals (variant) | 9 terms | n=10 | Same as above plus diagonal constraint. **HIGH**. |
| **A287647** | Min diagonal transversals | 9 terms | n=10 | Min variant; CP attack. **HIGH**. |
| **A287695** | Max diagonal Latin squares orthogonal to a given DLS | 9 terms | n=10 | Constraint solver. **MEDIUM-HIGH**. |
| **A309210** | Main classes of extended self-orthogonal DLS | 7 terms | n=8 or n=9 | Classification; canonical augmentation. **MEDIUM**. |
| **A003090** | Main classes of Latin squares of order n | 11 terms | n=12 | n=11 enumerated by McKay; n=12 a frontier. **LOW**. |

## Coding Theory Tables

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A005864** | A(n,4): max binary code length n, min dist 4 | 16 terms | many gaps within range, plus n>26 | Best-known bounds at codetables.de; clique formulation. **HIGH** for n=27,28; **MEDIUM** beyond. |
| **A005865** | A(n,6) | 16 terms | n=22+ | Same. **HIGH/MEDIUM**. |
| **A005866** | A(n,8) | 24 terms | n=25+ | Same. **MEDIUM**. |
| **A005851** | A(n,8,5) constant-weight | 22 terms | n=27+ | Constant-weight code clique attack. **HIGH**. |
| **A005852** | A(n,8,6) constant-weight | 17 terms | n=22+ | Same. **HIGH**. |
| **A005853** | A(n,8,7) constant-weight | 17 terms | n=24+ | Same. **HIGH**. |
| **A179183** | Optimal binary code, edit distance 3 | 8 terms | n=10 | Edit-distance graph clique; SAT viable. **HIGH**. |
| **A230380** | Optimal binary code, edit distance 4 | 8 terms | n=11 | Same. **HIGH**. |
| **A057657** | Max binary code correcting one transposition | 11 terms | n=12 | Clique on transposition graph. **HIGH**. |
| **A085680** | Constant-weight 2 transposition-correcting | several | n where last term added | Small; **HIGH**. |
| **A265032** | Max code length n, min dist 3 over alphabet 4 | 5 terms | n=6 | Quaternary clique. **HIGH**. |

## Covering Codes

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A000983** | Min binary covering code, radius 1 | 9 terms | n=10 | Classic. ILP or hitting-set. **HIGH**. |
| **A029866** | Min binary covering code, radius 2 | 8 terms | n=9 or n=10 | Same. **HIGH**. |
| **A004044** | Football pool: min covering in {0,1,2}^n radius 1 | 6 terms | n=7,8 | Famous; n=6 was recent. ILP or structured search. **MEDIUM-HIGH**. |
| **A060438** | Triangle: max binary code length n, covering radius k | 9 rows | row 10 | Tabular, builds on A000983. **HIGH**. |

## Crossing Numbers

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A000241** | Crossing number of K_n | 13 terms (verified through n=12) | n=13 | Known up to n=12; n=13 SAT/QUBO has been attempted. **MEDIUM**. |
| **A128422** | Projective plane crossing number K_{4,n} | several | next n | Easier than K_n. **HIGH**. |

## Sprague-Grundy / Game Theory

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A316533** | SG of Node-Kayles on Petersen P(n,2) | 22 terms | next n | Game-tree memoization; periodicity check. **HIGH**. |
| **A316632** | SG of Node-Kayles on 3 x n grid | 16 terms | n=17+ | Transfer-matrix or memoized DAG. **HIGH**. |
| **A344227** | SG of Node-Kayles on n-queens graph | 14 terms | n=14+ | Game tree on small graph; attack with alpha-beta plus memoization. **HIGH**. |
| **A006016** | Nim value of Sym game with n tails, 1 head | 26 terms | n=26+ | Standard SG recursion. **HIGH**. |
| **A006018** | Periods for game of Third One Lucky | 29 terms | next | SG recursion + period detection. **HIGH**. |
| **A046671** | Nim values G(3,n) Sylver coinage | 13 terms | n>16 | Combinatorial; very deep but small extensions possible. **MEDIUM**. |
| **A363934** | SG for Heat Toggle on n x k grid | tabular | larger grids | New game; small n,k tractable. **HIGH**. |
| **A075274** | Reachable configurations in Blet | 11 terms | n=12+ | BFS; **HIGH**. |

## Spectral Graph Theory

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A064731** | Connected integral graphs on n vertices | 12 terms | n=13 | Generate-and-test on graph database; nauty plus char poly. **HIGH** through n=13, **MEDIUM** beyond. |
| **A077027** | All simple integral graphs on n vertices | 11 terms | n=12 | Same. **HIGH**. |
| **A363064** | Connected Laplacian-integral graphs | 11 terms | n=12 | Same with Laplacian spectrum. **HIGH**. |
| **A243332** | Integral triangle-free connected graphs | several | next n | Same with triangle-free filter. **HIGH**. |

## Snarks / Cubic Graph Enumeration

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A130315** | Nontrivial snarks on 2n nodes | 18 terms | 2n=38 | Brinkmann/Goedgebeur snarkhunter codebase exists; current frontier 2n=36. **MEDIUM**. |
| **A216834** | Weak snarks on 2n nodes | 18 terms | 2n=38 | Same. **MEDIUM**. |
| **A175847** | Cyclically 4-connected cubic graphs on 2n | 13 terms | 2n=28 | nauty; **MEDIUM**. |

## Ramsey / Additive Combinatorics

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A006672** | Ramsey r(C_4, K_{1,n}) | 10 terms | n=11 | Small Ramsey; SAT-tractable. **HIGH**. |
| **A123938** | Ramsey r(K_{2,2}, K_{2,n}) | 13 terms | n=14 | Same. **HIGH**. |
| **A123939** | Ramsey r(K_{2,2}, K_{3,n}) | 9 terms | n=11 | Same. **HIGH**. |
| **A267295/6** | Circulant Ramsey RC_2(3,n) etc | 11 terms each | next n | Restricted to circulants; tractable SAT. **HIGH**. |
| **A005346** | Van der Waerden W(2,n) | 6 terms | n=7 | a(7) is at frontier; massive SAT attacks ongoing. **LOW** but worth tracking. |
| **A030126** | Schur numbers (version 1) | 5 terms | n=6 | a(5)=160 was huge SAT proof; a(6) frontier. **LOW**. |
| **A072842** | Largest m partitionable into n weakly sum-free sets | 4 terms | n=5 | a(5) is frontier (related to Schur). **LOW-MEDIUM**. |
| **A394031** | Sidon set in GF(2)^n | 10 terms | n=11 | Affine Sidon; clique on F_2^n. **HIGH**. |

## Domination & Independence

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A075458** | Queens domination Q(n) | 25 terms | n=26+ | Constraint solver; well-studied. **MEDIUM**. |
| **A075324** | Independent queens domination | 31 terms | n=32+ | Same. **MEDIUM**. |
| **A251419** | Domination on triangle grid TG_n | 40 terms | n=41+ | Likely periodic structure; transfer matrix. **HIGH**. |
| **A279402** | Toroidal queens domination | 22 terms | n=23+ | More symmetric than Q(n). **MEDIUM-HIGH**. |
| **A220840/1** | k-cop-win graphs of order n | 10 / 10 terms | n=11 | Graph database plus cops-robbers DP; **MEDIUM**. |

## Costas Arrays

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A008404** | Costas arrays of order n (with rotations) | 29 terms | n=30,31,32 | Frontier currently n=29; permutation backtracking, **MEDIUM**. |
| **A001441** | Inequivalent Costas arrays | 29 terms | n=30+ | Same. **MEDIUM**. |
| **A001440** | Symmetric inequivalent | 32 terms | n=33+ | Restricted symmetry; **MEDIUM-HIGH**. |
| **A213270** | Costas with involution permutation | 29 terms | n=30+ | **MEDIUM-HIGH**. |
| **A213271** | Costas with derangement | 29 terms | n=30+ | **MEDIUM-HIGH**. |

## Packing Chromatic / Coloring

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A335203** | Packing chromatic of n-cube Q_n | 8 terms | n=9 | SAT; Q_n has rich automorphism group. **HIGH**. |
| **A362580** | Packing chromatic of n x n grid | 10 terms | n=11+ | SAT; current SAT bounds tight. **HIGH-MEDIUM**. |

## Geometric / Packing

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A084617** | Max diameter-1 circles in n x n square | 16 terms | n=17+ | Continuous geometry; nonlinear opt. **MEDIUM**. |
| **A308578** | Max radius-1/n circles in unit square | 12 terms | n=13+ | Same. **MEDIUM**. |
| **A367557** | Fixed kissing polyominoes with n cells | 12 terms | n=13+ | Polyomino enumeration with adjacency constraint; **HIGH** at n=13. |

## Boolean Functions / Threshold

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A002077** | N-equivalence classes of self-dual threshold | 7 terms | n=8 | Threshold class enumeration; **MEDIUM**. |
| **A001206** | Self-dual monotone Boolean of n variables | 9 terms | n=10 | Dedekind-like; **LOW** for n=10 (frontier), **MEDIUM** for self-dual variant. |
| **A107765** | Inequivalent self-dual nondegenerate monotone | 6 terms | n=7 | **MEDIUM**. |

## Diophantine / Number Theory

| A-number | Description | Known | Next | Attack |
|---|---|---|---|---|
| **A332201** | Sum of three cubes: smallest |x| with n=x^3+y^3+z^3 | many; gaps remain at sporadic n | known frontier (n=42, 114, 165, ... were resolved 2019-2021) | Heuristic search; integer programming. **MEDIUM** for closing remaining gaps. |
| **A018899** | Smallest non-representable as sum of n distinct 2^a*3^b | 13 terms | n=14+ | Number-theoretic plus DP. **HIGH**. |
| **A118771** | Sum-set partition Schur-like | 4 terms | n=5 | Same as Schur. **LOW**. |

---

## Notes on Selection

- Sequences with `keyword:hard,more` and 5-30 terms are the sweet spot.
- Avoided sequences whose next term is known to require petascale resources (e.g., a(13) of A000241 in pure exhaustive enumeration; Schur a(6); van der Waerden W(2,7)). Listed as LOW with brief commentary so the user can decide.
- Many Latin-square and Costas-array sequences cluster around the same parameter horizon (n=10, n=30 respectively); a single SAT/CP framework could attack several at once.
- The Sprague-Grundy game-theoretic candidates are the cleanest "score quick wins" targets: short SG recursions, period detection, no symmetry-breaking required.
- For coding-theoretic A(n,d), the codetables.de tables list the best-known bounds; computing matching constructions provides the OEIS terms.

## Suggested Bundling

- **One-week SAT push (HIGH feasibility, related)**: A179183, A230380, A057657 (edit/transposition codes) plus A287644/8 (DLS transversals).
- **Game-theoretic memoization push**: A316533, A316632, A344227, A006016, A006018.
- **Coding-table sweep**: A005864/5/6, A005851/2/3, A265032 with a unified clique-on-Hamming-graph solver.
- **Spectral-graph sweep**: A064731, A077027, A363064, A243332 from a single nauty-based pipeline.
