# Saturation OEIS Harvest Assessment (May 2026)

Candidates from the published saturation work (paper 5, Zenodo
10.5281/zenodo.19226372) assessed for OEIS submission. Ranked by
(novelty x cleanliness) / collision-risk.

## Tier 1 -- submit: sat(n, P_4)

- DATA (n=4..): 2, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 9, ...
- Closed form: a(n) = n/2 (even), (n+3)/2 (odd); Kaszonyi-Tuza 1986.
- **Non-monotone** -- the genuinely interesting feature. No OEIS hit
  as a saturation sequence (Apr 2026 search).
- Submission file: `sat-P4-submission.txt`.
- Risk: a pure-integer interleaving like {2,4,3,5,4,6,...} could collide
  with a non-graph sequence. Re-search OEIS at submission time. If a
  collision exists, add a graph-saturation COMMENT + the KT1986 ref to
  the existing entry rather than creating a duplicate.

## Tier 2 -- compute small-n first, then maybe submit: sat(n, C_5)

- Formula (Chen 2009/2011): ceil(10(n-1)/7) - 1 for n in
  {11,12,13,14,16,18,20}, else ceil(10(n-1)/7). Proved only for n >= 11.
- Formula-extrapolated DATA (n=5..): 6, 8, 9, 10, 12, 13, 14, 15, 17, ...
- The irregular exception set makes this a genuinely interesting
  quasi-polynomial -- good OEIS material IF the small-n terms are
  verified by direct MaxSAT (the formula is NOT proved for n < 11).
- ACTION: compute sat(n, C_5) for n = 5..10 via the module, confirm or
  correct the extrapolated head, then submit. Pending compute.

## Tier 3 -- low priority, high collision risk: sat(n, K_{2,3})

- DATA (n=5..): 7, 9, 11, 13, 15, 17, 19, 21 = 2n - 3.
- Same shape as sat(n, K_4) (Erdos-Hajnal-Moon). Just the odd numbers
  >= 7; will collide with many sequences. OEIS likely rejects as "too
  simple" unless framed strictly as a saturation sequence.
- ACTION: hold. At most, add as a COMMENT/cross-ref to an existing
  odd-number sequence, noting the saturation interpretation.

## Also available (from harvest-oeis.md, non-saturation)

- A394445 b-file: 500 terms ready vs OEIS 43-term auto table. Pure win,
  already a published theorem. Submit the b-file. (Tracked separately in
  paper/STATUS.md as "Proposed Apr 27 2026".)
- Boolean-measures n=5 tuple-count sequences (2,3,6,12,24) and
  sensitivity=n class counts -- see boolean-measures-n5/oeis-candidates.

## Bottom line

One clean submit now (P_4), one worth a short compute then submit (C_5),
one to hold (K_{2,3}). The P_4 non-monotonicity is the headline.
