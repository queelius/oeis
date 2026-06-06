---
title: "OEIS Contributions"
date: 2026-06-05
summary: "Integer sequences contributed to the On-Line Encyclopedia of Integer Sequences, spanning Ramsey theory, extremal graph theory, knot theory, and algebraic graph theory."
tags: ["OEIS", "combinatorics", "research"]
---

# OEIS Contributions

**Why these.** These sequences have nothing in common at the level of
subject. They range over Ramsey theory, extremal graph theory, knot
theory, and the algebra of graphs. What they share is how they were found.
Each is a quantity with no known formula, computed exactly by SAT solving
or exhaustive search, pushed one value past the published frontier, and
checked against all prior art before it was submitted. The spread across
fields is deliberate: it is evidence that the method travels. And the
method has a direction. The computed values are where I go looking for
structure. The Rado numbers began as a column of computed integers and
ended as a proved theorem; the rest are trailheads, some already climbed,
others marked for whoever arrives next. The OEIS is where the results
live, as permanent reference points anyone can build on or check.

Sequences I have authored or extended in the [On-Line Encyclopedia of
Integer Sequences](https://oeis.org), validated against all published
prior art before submission (the working method of the
[open-problems](https://github.com/queelius) project). Author search:
[oeis.org/search?q=alex+towell](https://oeis.org/search?q=alex+towell).

## Authored sequences

These sequences were created and first populated by me.

| Sequence | Description | Domain |
|----------|-------------|--------|
| [A394445](https://oeis.org/A394445) | Distinct-variable 2-color Rado numbers for x+y=nz: the least k such that every 2-coloring of {1,...,k} has a monochromatic distinct-variable solution | Ramsey theory |
| [A394661](https://oeis.org/A394661) | Triangle T(n,k): number of prime knots with n crossings and three-genus k | Knot theory |
| [A395521](https://oeis.org/A395521) | Number of non-isomorphic abelian groups appearing as the sandpile group K(G) over graphs on n vertices | Algebraic graph theory |
| [A395644](https://oeis.org/A395644) | Number of fibered prime knots with n crossings | Knot theory |

A394445 is backed by a closed-form theorem (a proof that the Rado number
R(x+y=kz, 2; distinct) follows an explicit formula for all k >= 8), with a
500-term b-file. A394661 and A395644 come from a census of the KnotInfo
prime-knot tables; A395521 from an exhaustive sandpile-group computation.

## Extended sequences

These are classical Zarankiewicz-problem sequences (originally by N. J. A.
Sloane) where I computed and added new terms past the known frontier.

| Sequence | Description | My terms |
|----------|-------------|----------|
| [A006615](https://oeis.org/A006615) | z(n,n;3,4): least k forcing an all-ones 3x4 submatrix in an n x n 0/1 matrix | a(10)=67; a(11)=79 (in review) |
| [A006622](https://oeis.org/A006622) | z(n,n+1;3,4): same for n x (n+1) matrices | a(9)=61; a(10)=73 (in review) |
| [A006625](https://oeis.org/A006625) | z(n,n+2;3,4): same for n x (n+2) matrices | a(9)=67; a(10)=79 (in review) |

The Zarankiewicz extensions use SAT: a satisfying assignment exhibits a
dense matrix (a lower bound), and an unsatisfiable instance proves the
matching upper bound. The most recent terms (a(11)=79 for A006615, a(10)=73
for A006622, a(10)=79 for A006625) required a double-lex symmetry-breaking
encoding to make the upper-bound proof tractable: the plain encoding could
not settle z(10,11;3,4) in seven days, while the symmetry-broken version
proved it in 27 minutes.

## In review (June 2026)

- A006615 a(11)=79, A006622 a(10)=73, A006625 a(10)=79 (the three exact
  Zarankiewicz values above), proposed June 5 2026.

The A394445 b-file (500 terms) is already approved and live. An A042950
comment was proposed and then withdrawn (the matching terms were judged
insufficient evidence without a structural reason), and so is not counted
here.

## Reviewers

These submissions were reviewed and approved by OEIS editors including
Michel Marcus, Sean A. Irvine, Max Alekseyev, and Jon E. Schoenfield.

---

The OEIS has been a shared ledger of integer sequences since 1964. Adding
a term is a small, permanent contribution to a commons that outlasts any
one project: a number that was unknown is now known, checked, and citable
by anyone who needs it.
