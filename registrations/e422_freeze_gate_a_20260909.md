# e422 gate-(a) freeze artifact (2026-09-09, w1:pG under delegated authority)

This document is the FULL enumeration the frozen registration text points
to ("each content-SHA-256-bound in the freeze artifact").  It is
committed at the gate-(a) completion commit; any later change to any
object below requires a new E-number.

## Frozen object of record

- **Text:** docs/e422_registration.md
- **git blob:** `6e1b9ce128ac3a1be4678463bdc777c25b50ca8d` (cat-file
  verified `blob`; this is the exact object the citation resolves)
- **content SHA-256 of the identical bytes:**
  `dc7a362773fbf3a03e759aa7d178125963a14d88d9bbd61b66389a04ba8d7b3b`
  (42,079 bytes)
- **Freeze commit (text):** `9dd2380b8ea73e71494358393a4fb5bb4cbbaa07`
- **Citation-bind commit (module constant):** `ca8e3ea32058cd204b6293374a3a2e0a911eb4c0`
- **T0:** `2026-09-11T00:00:00Z` (exact second-level UTC)
- **Batch windows:** batch k covers (T0 + (k−1)·7d, T0 + k·7d], UTC;
  batch 1 = (2026-09-11T00:00:00Z, 2026-09-18T00:00:00Z];
  batch-4 evaluation day-28 (2026-10-09), batch-8 day-56 (2026-11-06),
  batch-12 day-84 close (2026-12-04).

## Bound implementation set (registration §6)

Content SHA-256 of each file's committed bytes at HEAD
`ca8e3ea32058cd204b6293374a3a2e0a911eb4c0` (all five verified identical
to the working tree at freeze):

| module | content SHA-256 | bytes |
|---|---|---|
| experiments/e422_protocol.py | `9b773b7bc1611d58255917b5f52faf836d6cc1ea5f808a59bab11435c845ec33` | 34,530 |
| experiments/e422_manifest.py | `6142d362abe7a1e1762c15b0247f9dd53f87bff5081598b2526bb87a2da72b16` | 13,935 |
| experiments/e422_status.py | `e10c529a73299a84986524655a6d9941fc9818facdd66a300508435402a9e96e` | 3,933 |
| experiments/e422_runner.py | `ac7729a53013b83712521707776c340fbd24700190cd0b1b169511d138e0e523` | 26,306 |
| experiments/e422_evaldriver.py | `9307dd1a00bfb1b464e04d73a680e990abdec3c7d1c6c50eaa11ba8cab69179b` | 15,200 |

Additional committed pre-freeze wiring (informational, NOT part of the
bound analysis set; it carries no registered rule — the registered rules
live in the five modules above and the frozen text):
experiments/e422_dispatch.py content SHA-256
`12211ca7728d237d159462d22f1152248af6ba45dee45f059be42baabf039a6f`
(4,048 bytes) — the gate-(e) CLI wrapper around e422_runner.run_batch.

## Bound registered-rule objects

| object | identifier | verification |
|---|---|---|
| analysis-grid source (PLDDT_BINS, §3 grids) | git blob `f55ac83de71ac2b0fb65a64a0a2c660e4e6842a9` | cat-file verified; content SHA-256 `628b411d915c6bcaee24375bff2ce129c2c3a0f59fed02e2b256a63e0e24d3eb` (15,648 bytes); identical to the worktree bytes of experiments/e421_conformal.py at freeze |
| label rules (§3 eligibility/labels) | docs/e421_label_rules.json content SHA-256 `e614590e82b5414b8c2763b2031ca670c444598ec3b54117600984ad4485bebd` | 1,587 bytes at HEAD |
| source allowlist | docs/e421_source_allowlist.json content SHA-256 `c4aeccd6061aff021eb04b3b7615c18cbc410c0445ad5201319d727cc481005b` | 1,968 bytes at HEAD |
| analysis C design of record (e421-P §11 executes as written) | docs/e421_prospective_registration.md whole-file content SHA-256 `6f26f5bbdf0b1503c8fe4bd34cdbf96cbe4599a534dc66f1806611ea691d3fc8` | 26,447 bytes at HEAD |

## Pre-freeze review records (bound for provenance)

- gate-(c) detached audit record: docs/audit_e422_gate_c_detached_audit_20260909.md
  content SHA-256 `f01f1d8a1af78844cf1df5b1b94ef9cf21312cdb0d9be6e60f9b09d11cc5a843`
  (18 findings, all RESOLVED; commit `a559b437`).
- §3B joint-grid detached review record:
  docs/audit_e422_joint_grid_review_20260909.md content SHA-256
  `a02b8e0dab9f1ba4fac71d012cfe827bc654b15a7bbee9701621e7df1b4df06f`
  (workflow `wf_75c31e8b-fbc`, 2 auditors READY_WITH_FINDINGS, 13
  findings — 0 BLOCKER/MAJOR — all disposed; commit `aa6eb3e8`).

## Dependency lock (verified at freeze)

Shared `.venv`: python 3.12.13 / numpy 1.26.4 / scikit-learn 1.9.0 /
scipy 1.17.1.

## Test state at freeze

e422 suite **117/117 passing** at the citation-bind commit
(`tests/test_e422_protocol.py`, `test_e422_evaldriver.py`,
`test_e422_manifest.py`, `test_e422_runner.py`, `test_e422_status.py`,
`test_e422_dispatch.py`), including: the T0 freeze-enforcement test
`test_default_t0_is_frozen_bind` (parses the registration for the
`T0 = …Z` bind and asserts the dispatch default matches) and the frozen
citation test `test_registration_citation_is_frozen_gate_a_object`
(asserts both labeled identifiers).

## Post-freeze discipline

Any change to any object enumerated here — registration text, the five
bound modules, the grid, label rules, allowlist, or analysis C's design
of record — is a NEW E-number with its own registration, review, and
anchor.  Corrections to THIS artifact (e.g. transcription errors) are
dated amendments to docs/DECISIONS.md only; the hashes above are
independently verifiable from the git objects and are authoritative.
