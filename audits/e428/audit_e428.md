# audit_e428 — independent check of the e428 exposure audit

> Reviewer phase, 2026-10-09, fresh context. Script `experiments/audit_e428.py`
> (`python3 -I`, offline, cached files only, pure Python). I wrote and ran it
> before opening `experiments/e428_exposure_audit.py` (never opened) or
> `docs/e428_results.md` (opened only for the comparison below).

## Scope and method

1. **T-P1 from raw inputs.** I parsed the unzipped workbook XML myself
   (sharedStrings + `PRED-RMSDLIST` via the workbook rels; xlsx SHA-256 matches
   the registration). Rows 8–109: AP = ID, AR = r.m.s.d. pLDDT from
   `Summary` B/H (215 IDs, no duplicates). Classes taken from
   `t_structures.json` after checking its IDs, r.m.s.d. and pLDDT against my
   parse. Bootstrap: within-group, 10,000 resamples, seed 12345. Permutation:
   10,000, seed 54321. My own RNG and code, not numpy.
2. **Classes.** I re-derived every class from `data/e428/http/` (all 2,233
   body SHA-256s verified), using RCSB entry/entity records (T) and SIFTS
   author-chain mapping plus the `seqres` sequence (A). I matched the dated
   search requests by sequence, identity cutoff and accession. The registered
   seed-7 sample (`random.Random(7).sample`, ceil 20 %) was 21 T and 17 A. I
   also rechecked all 184 because it was cheap.
3. Train-split overlap: I intersected the cached undated search hits with the
   `atlas_train.csv` PDB IDs and SIFTS accessions. I also checked the release
   dates, the fusion constructs and 7DRH.

## Recomputed vs reported

| quantity | audit | reported |
|---|---|---|
| T n / overall median AR | 102 / 0.954 | 102 / 0.954 |
| JSON vs xlsx (ID, r.m.s.d., pLDDT) | 0 mismatches | — |
| P1 n close / other | 24 / 77 | 24 / 77 |
| P1 medians | 0.847 / 0.956 | 0.847 / 0.956 |
| P1 diff | −0.109 | −0.109 |
| P1 95 % CI | [−0.426, 0.527] | [−0.412, 0.534] |
| P1 perm p | 0.381 | 0.376 |
| P1 verdict (CI upper < 0) | not supported | not supported |
| S2 medians / diff | 93.19 / 93.95 / −0.76 | same |
| S2 CI, p | [−2.82, 0.93], 0.316 | [−2.89, 0.91], 0.320 |
| S3 coef close (CI) | −0.254 [−0.518, 0.010] | −0.254 [−0.532, 0.011] |
| S4 r.m.s.d. medians SA/95/30/NOV | 0.757/0.873/1.2505/0.943 | same |
| T classes 9/15/28/49/1 | all 102 match (sample 21/21) | — |
| A classes 28/4/25/25 | all 82 match (sample 17/17) | — |
| A released ≤ 2020-05-01 | 40 | 40 |
| A train overlap SA / ≥95 / ≥30 | 3 / 5 / 10 | 3 / 5 / 10 |
| Fusions (2 SIFTS accessions) | 6xds_A MBP+TREM2, 6xrx_A MBP+Q8T5C5, 7p41_D T4-lysozyme (D9IEF7)+MARC1 | same |
| 7DRH | cached RCSB entry 404 "No data found" | 404 |
| No-accession chains | T: 7DTR, 7S5L; A: 7 chains | same |

All T entries use protein entity 1. No T chain has more than one accession.

## Discrepancies

1. **Fusion bound is mis-stated (conservative, not inflating).** Using only
   TREM2 (Q9NZC2), 6xds_A has 3 training-era same-accession entries, so it is
   SAME_ACCESSION without its MBP tag. Only 6xrx_A (MBP only) and 7p41_D
   (D9IEF7 = T4 lysozyme; MARC1 has 0 accession hits) depend on the tag for
   the accession criterion. The lower bound should be 30/82 (37 %), not 29.
   Caveat 2 and the falsification line ("may come from the tag" for all three)
   should be corrected. The train-split reading (2 genuine overlaps) stands.
2. **The 7DRH removal/8ITE replacement has no source in the e428 artefacts.**
   No cached response, receipt or `run.log` line supports it. The 404 itself
   is genuine. Either cache that holdings record or mark the statement as
   unverified.
3. **HTTP counts.** The 2,258 figure counts receipts and includes 25 cache
   re-hits. There are 2,233 unique cached responses (200 ×1,872, 204 ×340,
   404 ×21). "SHA-256 of all 2,258 cached responses" should say receipts.
4. **Minor wording.** "No difference in confidence" (S2) should read "no clear
   difference" (CI [−2.9, 0.9]). The bottom line "exposure does not explain
   the low r.m.s.d." does not follow from P1, which is uninformative (CI ≈
   ±0.5 Å). It is supported descriptively by the NOVEL median (0.943 Å, n 49),
   which is close to the overall 0.954. Cite that number instead.

Small CI/p differences come from the different seeds and RNG, as expected.

## Verdict

**PASS-with-notes.** Every load-bearing number and all 184 classes reproduce
exactly from raw and cached inputs. The registered T-P1 verdict (not supported)
is unchanged under an independent seed and implementation. Before promotion,
correct notes 1–3.
