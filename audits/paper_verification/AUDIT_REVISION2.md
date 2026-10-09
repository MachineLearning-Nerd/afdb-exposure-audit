# Independent audit: revised temporal-leakage paper

**Verdict: FAIL**

The P5 arithmetic reproduces, but the residue labels used by that evaluation are not reliable. I independently confirmed a residue-numbering bug in the e422 mapping path on the paper’s K7PQ54 example, 9qj6: the stored low lDDT values are reproduced by the mixed-number mapping, while a mapping in the coordinate system returned by SIFTS gives near-perfect AFDB/reference agreement. Until affected labels are corrected and P5 is rerun, the registered coverage failure is only a result on the current ledger; it cannot support a scientific conclusion about coverage on valid reference labels. A ledger-wide screen also finds other unresolved candidate anomalies.

## Scope and audit boundary

I audited the current working-copy revision of `paper/arxiv/temporal_leakage.tex` and `references.bib`, with priority on task 4b (K7PQ54) and task 4 (independent P5 reproduction). I then checked the other verification tasks in `paper/arxiv/verification/CHECKER_BRIEF.md`: V1 recount, V2 raw-cache sample, V2b e421 release dates, numeric claim tracing, claims and source support, bibliography resolution, adversarial review, and the TeX build.

At audit start, `HEAD` was `923948836d9242701a76fe7074ba4574e5336d0e`. The paper and bibliography were modified in the working tree, and the revision log and verification materials were not committed. This report therefore audits the working-copy paper snapshot, not a clean committed release. Snapshot checksums are listed below. No paper, `docs/`, `results/`, or e422 source/result file was written by this audit. Audit scripts, caches, snapshots, and outputs are confined to this `audit/` directory; the TeX build and regenerated figures went to `/tmp`.

The required Herdr environment check did not pass (`HERDR_ENV` was not `1`). I did not inspect or control any panes or contact other sessions. The task’s write boundary precluded writing a status file outside this audit directory; this execution note records the blocker here.

## Priority findings

### F1 — Confirmed mapping-coordinate defect for K7PQ54 / 9qj6

**Severity: critical; verdict: artifact confirmed for 9qj6, scope beyond this entry unresolved.**

The raw PDBe SIFTS response maps each 9qj6 chain from UniProt 35–507 using `start.residue_number = 1` and `end.residue_number = 473`. These values are in the mmCIF `label_seq_id` coordinate space. The first observed residues are label 2 / author 36 / UniProt 36 for chain A and label 3 / author 37 / UniProt 37 for chain B. AFDB’s SIFTS cross-reference positions agree with these UniProt positions.

The e422 runner at `experiments/e422_runner.py:263–268` takes `author_residue_number` when present and falls back to `residue_number`, then constructs one continuous residue map. The e421 label parser at `experiments/e421_pipeline.py:62–79` looks up atoms using `auth_seq_id`. For 9qj6, the mapping therefore starts in `label_seq_id` space and is applied to `auth_seq_id`. In the runner’s mixed-number map, chain A author residue 36 is assigned UniProt 70; chain B author residue 37 is assigned UniProt 71. This is the offset reported in the paper’s low-lDDT rows.

I independently parsed the cached PDBe mmCIF, SIFTS response, and AFDB v6 CIF and aligned Cα atoms using the SIFTS label-space mapping. Results:

| Chain | Aligned Cα | Sequence mismatches | Median pLDDT | Median Kabsch residual | Median lDDT |
|---|---:|---:|---:|---:|---:|
| A | 472 | 0 | 98.50 | 0.3833 Å | 0.992908 |
| B | 471 | 0 | 98.50 | 0.3433 Å | 0.995000 |

Replaying the runner’s mixed-number map independently reproduces the stored batch-4 low-lDDT groups to within `0.00005` median lDDT: chain A has median lDDT `0.197674` and median displacement `20.7202 Å`; chain B has `0.199074` and `20.7092 Å`. This directly explains the apparent structural disagreement for the tested entry. Input hashes and residue-level calculations are in [k7pq54_9qj6_audit.json](k7pq54_9qj6_audit.json); the source files are cached under `cache/`.

This confirms a pipeline-coordinate mismatch for 9qj6, rather than only an inference from the internal chain inconsistency. It does not by itself establish the correct mapping for every K7PQ54 chain in the 19-entry batch. The runner needs an explicit coordinate-space-aware mapping, and affected labels need to be rebuilt and independently checked before the complete K7 result is described as repaired.

### F2 — P5 is numerically reproduced but is not scientifically interpretable on the present labels

**Severity: critical; verdict: arithmetic reproduced, registered data outcome not validated.**

I wrote a fresh calculation from the e422 registration’s §3A.1 specification, using 200 protein-disjoint splits, seeds 0–199, the registered calibration quantile, and inclusive interval coverage. I did not import project code. The independent means match the registered evaluation closely:

| Nominal level | Registered mean (MCSE) | Independent mean (MCSE) | Difference | Registered threshold | Ledger outcome |
|---|---:|---:|---:|---:|---|
| 0.90 | 0.857965 (0.008007) | 0.858036 (0.008005) | +0.000071 | 0.875980 | Fail |
| 0.95 | 0.915196 (0.005718) | 0.915276 (0.005715) | +0.000080 | 0.932845 | Fail |

The paper’s rounding and threshold arithmetic are correct. The independent computation also reproduces the post-hoc K7PQ54-excluded residue-weighted means (`0.879204` and `0.934207`) and the reported calibration-width effect. In 44 splits with K7PQ54 in calibration, the mean half-width is `0.825739`; the other test proteins’ protein-averaged coverage is `0.97235` at 0.90. The miss shares are `49.84%` at 0.90 and `55.68%` at 0.95 in the current ledger.

Those matches validate the calculation against the ledger and registration, not the truth of the labels. F1 establishes that at least the cited 9qj6 low scores result from a mapping artifact. Because K7PQ54 is 14% of the ledger and about half the misses, this is consequential for P5. The current ledger failure can remain recorded as the registered computation on the submitted ledger, but the paper must not present it as an outcome on correctly mapped reference labels until the relevant labels are repaired and the registered computation is rerun. The exclusion analysis is post hoc and is not a substitute for that rerun; it also does not meet the registered nominal-minus-three-MCSE criterion at the displayed point estimates.

The independent output and method details are in [p5_independent_reproduction.json](p5_independent_reproduction.json). No P6/P7 conditional-coverage or Mondrian/CQR comparison was independently computed; the current audit output explicitly records that limitation.

### F3 — More suspicious e422 mappings remain unresolved

**Severity: major; verdict: unresolved.**

I screened 117,905 valid e422 rows across 138 accessions using a predeclared transparent flag: at least 30 rows, median pLDDT ≥90, and median lDDT <0.4 at chain level. The screen flagged 116 chain groups and 28 accession groups. A second screen found 13 accession/entry groups with eligible chains whose median lDDT differs by at least 0.4. These are audit candidates, not independently confirmed mapping errors.

The author-flagged examples all remain conspicuous. Q5JF22 / 9sq1 is especially informative: chain A has median lDDT 0.9853, while chains B–D are around 0.185–0.189, despite similarly high pLDDT. Other flagged patterns include low lDDT across both chains for P07342 / 9y8y and A0A4Q0WMG2 / 9wvj. Some may reflect biological construct/domain differences, genuine structural disagreement, or another artifact; the ledger alone cannot distinguish these explanations. The full candidate list, summary rules, and source hashes are in [e422_signature_screen.json](e422_signature_screen.json).

These unresolved rows prevent attributing residual P5 undercoverage after K7 exclusion to protein-level nonexchangeability, exposure, or a method limitation. The revised paper appropriately labels several summaries exploratory, but its causal interpretation must remain suspended pending mapping checks and a corrected-label rerun.

## Remaining verification tasks

| Task | Verdict | Independent result and limitation |
|---|---|---|
| V1 channel recount | **PASS** | Recounted unique `(accession, PDB entry)` pairs from the e420 temporal artifacts. Admitted pool: 599 accessions / 646 pairs; training channel 553 accessions / 595 pairs; template bound 599 / 646. Overflow: 269 / 289; training 238 / 256; template bound 264 / 284. Union: 868 / 935; training 791 / 851 (91.0% of accessions); template bound 863 / 930 (99.5% of pairs). The 72 training/template classification differences are all template-only. These match the paper’s table. See [v1_channel_recount.json](v1_channel_recount.json). |
| V2 raw-cache sample | **PASS-with-notes** | Re-derived classifications from request-keyed raw RCSB cache bodies for 31 accessions: 25 uniformly sampled from the sorted V2 union, two required proteins, and five V2 training-NOVEL proteins. No class or raw-count disagreements. Table 2 aggregates were independently recounted from the supplied 199-accession V2 output; the raw queries were not replayed for all 199. See [v2_cache_sample_audit.json](v2_cache_sample_audit.json). |
| V2b e421 release dates | **PASS-with-notes** | All 168 qualifying entries had a cached RCSB release date later than their recorded `modelCreatedDate`; none were missing or at/before that date. Of 96 successful entries (24,457 residue rows), none were at/before it. This verifies the stated date comparison, not that `modelCreatedDate` is the date model inference ran. See [v2b_e421_audit.json](v2b_e421_audit.json). |
| E421 numeric recount | **PASS-with-note** | Recounted all error-bin counts/rates and Spearman values from the labeled residue artifact. For pLDDT 90–100, raw `n=20,778` includes one nonfinite lDDT; the valid denominator is 20,777 and the rate is `121/20,777 = 0.005824`. The rounded paper rate remains 0.6%. Spearman correctly uses 24,456 finite pairs. See [e421_numeric_recount.json](e421_numeric_recount.json). |
| Claim/source check | **FAIL** | The main date-channel distinction and eligibility caveat are appropriately stated, but two AFDB date statements need correction/qualification. The cited human-proteome paper reports PDB and PDB70 snapshot dates in its Methods ([Tunyasuvunakool et al.](https://www.nature.com/articles/s41586-021-03828-1)); saying the cited source does not state the template cutoff is too broad. Also, the official AFDB metadata documentation defines `modelCreatedDate` as the date of creation for the entry, not a documented inference or v6 relabel date ([AFDB metadata documentation](https://github.com/google-deepmind/alphafold/blob/main/afdb/README.md)). The date comparison should be presented as a metadata-based bound/proxy, not a verified run date. PDBe’s release notes support the statement that unchanged v4 coordinates were carried into the v6 relabel ([PDBe release notes](https://www.ebi.ac.uk/pdbe/news/alphafold-database-release-notes)). The official AlphaFold training cutoff source supports the 2018-04-30 cutoff ([DeepMind technical note](https://github.com/google-deepmind/alphafold/blob/main/docs/technical_note_v2.3.0.md)). |
| Bibliography and references | **PASS-with-notes** | All 20 DOI-bearing entries resolved to the cited works through DOI/Crossref/publisher checks; official URL/arXiv references also resolve. `durairaj2024` and `vovk2005` have no DOI or URL in the bibliography although stable identifiers are available (PLINDER DOI `10.1101/2024.07.17.603955`; Vovk book DOI `10.1007/b106715`). The `bertoni2026` and `haas2018` author lists place `and others` before a final named author, so the built bibliography literally prints “others, and [author]”. Put `others` last or list the authors. See [bibliography_doi_audit.json](bibliography_doi_audit.json). |
| TeX build and figures | **PASS-with-notes** | Rebuilt the exact paper snapshot and regenerated all three figures in `/tmp`; regenerated figure hashes exactly match the paper’s checked-in PNGs. Tectonic produced a 9-page PDF with resolved citations and figures. Final log had no errors, unresolved citations, or overfull/underfull boxes. It retained a cross-reference rerun warning and changed two `[h]` placements to `[ht]`. |
| Adversarial pass | **FAIL** | A credible pipeline artifact is now confirmed for the central anomaly example, and a ledger-wide scan finds other unexplained high-pLDDT/low-lDDT and within-entry chain discrepancies. Their alternatives remain open. No defensible conclusion about the full P5 residual survives until they are checked. |

### Claim language and interpretation

The draft correctly avoids claiming that historical eligibility proves a specific structure was used for training, and it labels the K7 exclusion and exposure-stratified analysis post hoc. Keep those qualifications. After this audit, the 9qj6 diagnosis itself should be described as **confirmed for the tested entry**, while the extent of the bug across K7PQ54 and e422 remains unresolved. The methods discussion should identify the `label_seq_id`/`auth_seq_id` convention directly.

The post-K7 shortfalls are 2.08 percentage points at 0.90 and 1.58 points at 0.95, not approximately two points at both levels if precision matters. More importantly, those are current-ledger exclusion summaries, not corrected-label results. The exposure-stratified point estimates are descriptive; they do not establish either an exposure effect or equivalence. The conclusion should not explain the residual shortfall through protein-level nonexchangeability until valid labels and an appropriate protein-level analysis support that explanation.

The current sentence that `modelCreatedDate` records when a prediction was “produced or relabelled” is not supported by the AFDB metadata definition. Replace it with a qualified description of the recorded metadata date and cite the exact source supporting any claimed date semantics. Similarly, describe the Tunyasuvunakool database snapshots rather than saying the source gives no template date information.

## Required resolution before a passing revision

1. Fix the mapping code to preserve coordinate-space identity; do not apply a `label_seq_id` offset to `auth_seq_id` values. Add a focused regression fixture for 9qj6 and author/label numbering offsets.
2. Rebuild and independently inspect all K7PQ54 chains and the other flagged e422 candidates, retaining explicit unresolved/error states rather than silently keeping or dropping rows.
3. Freeze corrected labels and provenance, then rerun the registered P5 calculation unchanged. Report the original ledger result as the original ledger result and corrected-label result separately.
4. Revise claims about the cause of residual undercoverage and exposure dependence only after the corrected-label result and adversarial review.
5. Correct the AFDB date wording and the two bibliography author-list formats; add stable identifiers for the PLINDER and Vovk references.

## Audit artifacts and reproducibility

Independent scripts and machine-readable output in this directory:

- `k7pq54_mapping_audit.py` / [k7pq54_9qj6_audit.json](k7pq54_9qj6_audit.json)
- `e422_signature_screen.py` / [e422_signature_screen.json](e422_signature_screen.json)
- `p5_independent_reproduction.py` / [p5_independent_reproduction.json](p5_independent_reproduction.json)
- `v1_channel_recount.py` / [v1_channel_recount.json](v1_channel_recount.json)
- `v2_cache_sample_audit.py` / [v2_cache_sample_audit.json](v2_cache_sample_audit.json)
- `v2b_e421_audit.py` / [v2b_e421_audit.json](v2b_e421_audit.json)
- `e421_numeric_recount.py` / [e421_numeric_recount.json](e421_numeric_recount.json)
- `bibliography_doi_audit.py` / [bibliography_doi_audit.json](bibliography_doi_audit.json)

The scripts were run with the repository virtualenv Python in isolated mode and from `/tmp`; project modules were not imported by the independent K7, P5, V1, V2, or V2b calculations. The bibliography check queried DOI/Crossref records and checked publisher pages for the Crossref rate-limited entries. The V2 and V2b analyses used cached, request-keyed response bodies and rechecked their hashes.

Snapshot/input hashes:

| Input | SHA-256 |
|---|---|
| Paper snapshot `temporal_leakage.tex` | `0f61eba0a3654b10f1323b0bbaa028ce9ee46eab26ff72526d6b16eaafb8a946` |
| Bibliography snapshot `references.bib` | `938c866994a145e11d6639516d8c6cf50ce0f59236355e2a90c0edd4e822635e` |
| Revision log snapshot | `f31d5c37db469bb7e40aaee182a8c2c7e6fdc5b2161d13055b565b2bc428726c` |
| PDBe 9qj6 mmCIF | `9e85341cb02c98cbfea6f22fba664e335247899dfbd3100f5df358090f5f639e` |
| PDBe 9qj6 SIFTS mapping response | `f242c6c1d15b1dbd87ad2245ccacb7714be409cbf2882187490d867ab0d94a8f` |
| AFDB K7PQ54 v6 model | `07406b15c127dc43b5c5febdd5e6766637ec15ed1e3a7cf3f6e43edf7fa546c9` |
