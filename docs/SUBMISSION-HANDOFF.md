# OEIS submission handoff

Resume point for submitting the candidate-new sequences in this repo. Five
sequences are drafted, reproducible, and OEIS-novelty-checked; they are not yet
submitted. Work through them when you open a fresh session here.

## Before you submit (per sequence)

1. Re-run `python3 sequences/<id>/reproduce.py` to reconfirm the DATA.
2. Re-do the OEIS novelty check (the database changes): search the DATA and a
   distinctive sub-run at `https://oeis.org/search?q=...&fmt=text`. The
   `fmt=json` endpoint throttles to false negatives; use `fmt=text`.
3. Confirm conventions: variable `n` (not k), one reference per line
   (alphabetical), no EXTENSIONS field on a NEW sequence, correct OFFSET.
4. Submit via the OEIS web form logged in as `queelius` (same flow used for the
   Zarankiewicz terms: edit -> fill NAME/DATA/OFFSET/COMMENTS/KEYWORD/AUTHOR ->
   Save Changes -> "These changes are ready for review by an OEIS Editor").
5. After acceptance: set `status: published` in `contributions.yaml`, add the
   A-number to the sequence's `metadata.yaml`, run `python3 scripts/render_readme.py`.

## The five candidates (ready)

| id | title | offset | DATA | draft | keyword |
|----|-------|--------|------|-------|---------|
| **S-Zn-2** | Maximum 2-colorable sum-free subset of Z/nZ | 2 | 1,2,3,4,4,4,6,6,8,8,9,8,9,12,12,12,12,12,16 | `sequences/S-Zn-2/oeis-submission.txt` | nonn |
| **S-Zn-3** | Maximum 3-colorable sum-free subset of Z/nZ | 2 | 1,2,3,4,5,6,7,8,9,10,11,12,13,13 | `sequences/S-Zn-3/oeis-submission.txt` | nonn,more |
| **sat-P4** | sat(n, P_4): min edges in a P_4-saturated graph on n vertices | 4 | 2,4,3,5,4,6,5,7,6,8,7,9,8,10,9,11,10,12,11,13 | `sequences/sat-P4/sat-P4-submission.txt` | nonn,easy |
| **crn-k2** | # chemical reaction networks on 2 species (up to species perm.) by reaction count | 1 | 8,99,570,2445,7752,19440,38760,63090,83980,92504,83980,63090,38760,19440,7752,2445,570,100,10,1 | `sequences/crn-k2-by-reactions/C1-crn-k2-by-reactions.txt` | nonn,fini,full |
| **max-tail-bivariate** | Triangle: functional graphs on [n] by (max tail length, # cyclic points) | 1 (read by rows) | see draft (flattened triangle) | `sequences/max-tail-bivariate/bivariate-README.md` | nonn,tabf |

Notes per sequence:
- **S-Zn-2 / S-Zn-3**: candidate-new (checked Jun 2026; no full-prefix match, and
  the distinctive non-monotone dip 9,8,9 at n=13 for k=2 is absent from OEIS).
  S-Zn-2 is order-invariant across abelian groups of order n (verified n<=20);
  S-Zn-3 is NOT (state this in the comment). Reproduce ~20s / ~2s.
- **sat-P4**: candidate-new; non-monotone; closed form Kaszonyi-Tuza 1986
  (n/2 if n even, (n+3)/2 if n odd). The genuinely interesting one is the
  oscillation. From computational-explorations? No: from this repo's
  graph_saturation.
- **crn-k2**: complete 20-term finite sequence, two-method confirmed
  (direct census + Burnside). Companion drafts C2 (k=3) and C3 (k=4) are in the
  same dir as OPTIONAL follow-ups (longer, same definition). OEIS has no CRN
  deficiency/enumeration sequences, so this convention (species-permutation only)
  is genuinely new; state the equivalence precisely so a reviewer can reproduce.
- **max-tail-bivariate**: the joint triangle is candidate-new; its marginals are
  A216242 (max tail length) and A066324 (cyclic points). Submit as `tabf`.

## Do NOT submit

- **Sidon-squares** (max Sidon subset of {1,4,9,...,n^2}): already in OEIS as
  **A390813**, which holds the correct larger values (the local heuristic
  undercounts past n~25). Ruled out.
- **A042950 comment**: withdrawn (Irvine: 5 matching terms insufficient without a
  structural reason).

## Reviewer context

Prior approvers for this author: Michel Marcus, Sean A. Irvine, Max Alekseyev,
Manfred Scheucher, Jon E. Schoenfield. The recently-accepted Zarankiewicz terms
(A006615 a(11)=79, A006622 a(10)=73, A006625 a(10)=79) went through cleanly.

## Also worth doing from this repo (new-work backlog)

`docs/harvest-oeis.md` and `docs/oeis-candidates-survey.md` list further
candidate targets to compute and submit. The toolkit in `src/` is the starting
point for hammering out new sequence attacks.
