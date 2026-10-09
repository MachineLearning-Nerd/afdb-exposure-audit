# Changes from the original project repository

This release was assembled from the authors' project repository at commit
`7c95d4fd3555c50967b663b9b8d589f274fe625c` (2026-10-09). Its git history is
not included.

## Summary

- **No copied file was modified.** Every file listed below is byte-identical
  to its source (SHA-256 shown). No registered/frozen code, registration,
  OpenTimestamps-stamped document or derived data file was edited — not even
  to remove absolute paths or internal references — so that every recorded
  digest and timestamp proof remains verifiable.
- **Layout.** Files moved to new top-level directories (source path in the
  table): `experiments/` → `code/experiments/`, `tests/` → `code/tests/`,
  `paper/arxiv/verification/*.py` → `code/verification/`, audit scripts →
  `code/audit/`, `docs/` → `registrations/`, data → `data/derived/` (which keeps
  the original relative paths). `code/tools/make_workroot.py` recreates the
  original layout in a scratch copy so that the registered code runs
  unchanged; no path in any registered file was rewritten.
- **Stamped earlier versions.** Where a record was amended after it was
  OpenTimestamps-stamped, the exact stamped version was recovered from git
  history (commit shown) and shipped next to its proof. The `.ots` files were
  renamed where needed so that `ots verify X.ots` finds `X`; `.ots` bytes are
  those of the source (19 proofs were later replaced by their upgraded
  versions; see "Final paper" below).
- **Request caches.** The request-keyed HTTP caches are copied in full, all
  files byte-identical: `data/paper_verification/http` (1,096: 897 RCSB,
  199 UniProt FASTA), `paper/arxiv/verification/cache/v2c/http` (136: 111
  RCSB, 25 UniProt) and `paper/arxiv/verification/cache/v5/http` (815: 403
  PDBe SIFTS mapping, 412 RCSB entry responses; 8.3 MB). The first release
  shipped only the 1,008 RCSB files; the 224 UniProt and 815 V5 files were
  added in the revision below. The V5 run logs in the source cache directory
  are not copied.
- **Withheld.** The batch-4 evaluation artifacts are replaced by the
  key-selected extract `data/derived/p5_batch04_outcomes.json` (see README).

## New files (written for this release)

`README.md`, `LICENSE`, `LICENSE-CC-BY-4.0` (Creative Commons legal code,
fetched from creativecommons.org), `CITATION.cff`, `requirements.txt`,
`.gitignore`, `CHANGES_FROM_ORIGINAL.md`, `MANIFEST.sha256`,
`code/reproduce_p5.py`, `code/paper_numbers.py`, `code/tools/make_workroot.py`,
`code/tools/extract_p5_outcomes.py`, `data/fetch/fetch_raw.py`,
`data/fetch/README.md`, `data/derived/p5_batch04_outcomes.json` (output of
`extract_p5_outcomes.py`), `audits/README.md`.

## Revision of 2026-10-09 (after three independent code reviews)

No registered, frozen or historical file was edited; all changes are new
files or files written for this release.

- `code/paper_numbers.py` rewritten: 160 checks, each labelled RECOMPUTED
  (114), RECONCILED (34) or RE-READ (12); recomputation from row- and
  request-level data instead of re-reading stored outputs; no silent skips
  (the number of checks is asserted; missing inputs exit with status 2).
- `code/reproduce_p5.py`: verifies the ledger checkpoints' canonical digests
  and `MANIFEST.sha256` entries before computing; compares mean, MCSE and
  threshold with the registered evaluator at 1e-12; explicit exit codes
  instead of `assert`.
- New: `code/tools/extract_v5_template_table.py` and its output
  `data/derived/paper/arxiv/verification/v5_afdb_template_table.json` (the
  template list, software and target of the 138 AFDB model files read by V5;
  the model files themselves, 38.6 MB, are still re-fetched).
- New: the UniProt and V5 cache files listed under "Request caches".
- `audits/census_correction/` and `audits/paper_verification/` (27 files)
  are byte-identical copies of the source repository's `results/e427/audit/`
  and `paper/arxiv/verification/audit/` files of the same names (checked
  2026-10-09; digests in `MANIFEST.sha256`); they were present in the first
  release's working tree but not listed here.
- New: `audits/label_crosscheck_biotite/` (biotite cross-check of the
  corrected labels: scripts, outputs, the SHA-256 of the PDBe updated mmCIF
  files used, README).
- Documentation: `README.md` (title, paper_numbers classes, V3/V6 usage,
  paths, blinding scope), `audits/README.md` (per-script table, boundary
  convention, verify_digests warning), `data/fetch/README.md`, `CITATION.cff`
  (title).
- Not changed, documented instead: `code/verification/v1_channel_recount.py`
  relies on the label-run receipt having `ledger_rows == 0`; it is kept
  byte-identical and `code/paper_numbers.py` checks the condition.

## Revision of 2026-10-09 (final paper)

No registered, frozen or historical file was edited.

- **OpenTimestamps proofs upgraded.** All `.ots` proofs of the source
  repository were upgraded with `ots upgrade` on 2026-10-09 (source commit
  `c261d0fe`), which adds Bitcoin block-header attestations; the stamped file
  digest recorded in every proof is unchanged (checked with `ots info` against
  the first-release proofs). The 19 proofs below were shipped in their
  pre-upgrade form and are now replaced by the upgraded source files, so their
  own SHA-256 differ from those of the first release; the table gives the new
  values. Every shipped `.ots` (36, including the two in
  `audits/census_correction/`) now carries at least one Bitcoin attestation
  and is byte-identical to the source at `c261d0fe`. First-release
  (pre-upgrade) SHA-256:
  - `data/derived/results/e427/batches/e422_batch02_checkpoint.json.ots`: `89af19f07d8f27b7e4ff654a1cb5733875a5bba6b0051cf58a3325f3c3eb270f`
  - `data/derived/results/e427/batches/e422_batch03_checkpoint.json.ots`: `743e5a37c95c68082235306068d45b487ac21cf0eeb99f0f3fb46eb6039807fe`
  - `data/derived/results/e427/batches/e422_batch04_checkpoint.json.ots`: `a4ffa7de3979501755edb1b7c6559cc4bf8af753036fc3472f3f3c2691043f63`
  - `data/derived/results/e427/batches/e422_binding_state.json.ots`: `c47e63e29012affc1a62459b5301fc2e6bcc95ed0003962d446fc77101deb63d`
  - `data/derived/results/e427/e421_relabel/census_results.json.ots`: `89c4cd8547ac8628726011f80ae5932fc627dca000f1ebc865cf1876e4549667`
  - `data/derived/results/e427/e421_relabel/receipts.json.ots`: `16cd3e5cd2be5c6595fcaca20f66a937147aff4d0e1dce8ac7cf0015d35f9098`
  - `data/derived/results/e427/e421_relabel/summary.json.ots`: `6b00aac4f1ccf88df7fcb935f698b84ee9186e4a964ec25939e3e274fe9b0e0d`
  - `data/derived/results/e427/evals/e422_batch04_evaluation.json.ots`: `026e927c93d35dabd66dec56abcfae11669b754008e6c590c918d4586ea7aac6`
  - `data/derived/results/e427/evals/e422_batch04_execution_receipt.json.ots`: `91c505d27845e0845b4e18303d94620261a285eba1179a169a95efa9963836ef`
  - `data/derived/results/e427/replay_report.json.ots`: `af3fdbc691b10ba68fa30b79f6a70f64d142929a047870fd6885d61fed132a29`
  - `registrations/e421_post_run_checkpoint_v2.json.ots`: `e092fabfab3e0e11b1e3a7a8fc33f3abccf3a6990fa2f1fdfc90a0b0d225219c`
  - `registrations/e421_post_run_checkpoint_v3.json.ots`: `8bcac02f528056f94861b76d1049ffcd327c2948a3f65df4e65c961253c95483`
  - `registrations/e421_pre_source_checkpoint.json.ots`: `32fbee6f98d4187fd88cf4facab5670b5de78d9fda6e549b56960e0bc10cfb9b`
  - `registrations/e422_correction_of_record_20261007.md.ots`: `6a0bd15350e734257cf1c70f13f82b1397a0c3d285616b1879dfb95ded1ac43d`
  - `registrations/e427_registration.md.ots`: `6839a50ee30171e005fd1cf11649f79096592b978a8d7b0cd8456887d3e853ae`
  - `registrations/stamped_versions/e422_correction_of_record_20261007_add20261008.md.ots`: `1bd59f18d92ea56598c8ae4dd31dbd354c580a169386485121caf7b2e37c7146`
  - `registrations/stamped_versions/e422_correction_of_record_20261007_add20261009a.md.ots`: `74cc1b8a07a019718d4618495a0d641116df1882f8a894408f963c8e929bce64`
  - `registrations/verification/REGISTRATION_A3_24a51a94.md.ots`: `36c7c8f90e2421a9c2424b0559d33d304461a4a36bbde1106c014a2a95f4d919`
  - `registrations/verification/REGISTRATION_A4_de0e7a32.md.ots`: `05070eb34657eb2a2040797e902c39dceff3012c090dbb602ad6c365483b1387`
- **Final paper.** `paper/temporal_leakage.tex` and `paper/references.bib`
  replaced by the final versions (source working tree, 2026-10-09; not yet
  committed there); new `paper/figures/rev4/f1_mechanism.pdf`,
  `paper/figures/rev4/f2_scope.pdf` and `code/verification/make_figures_rev4.py`
  (styling-only revision of `make_figures_rev3.py`, kept alongside it; run in
  a work root it reproduces `figures_rev3_numbers.json` byte-for-byte and the
  two PDFs pixel-identically at 150 dpi, only `/CreationDate` differing).
  `paper/temporal_leakage_paper.pdf` is built here from the tex with tectonic
  (no undefined references, no overfull boxes).
- `code/paper_numbers.py`: 208 checks (RECOMPUTED 147, RECONCILED 43,
  RE-READ 18), covering every quantitative statement of the final paper
  except the dates and design parameters listed in the README.
- New: `code/tools/extract_e421_model_table.py` and its output
  `data/derived/results/e427/e421_relabel/afdb_model_table.json` (model
  type, software and SHA-256 of the 93 AFDB model files downloaded by the
  post-2022 relabel; no coordinates or pLDDT).
- `.gitignore`: exception so that `audits/label_crosscheck_biotite/full_ledger_out.log`
  is shipped.
- **Revision 5 (2026-10-09).** `paper/` replaced by the source's revision 5
  (git `52345ba7`): wording revised section by section for readers outside
  the project, no number changed; the ATLAS step and the registration
  repository described more precisely. `code/verification/make_figures_rev4.py`:
  Figure 1 labels in plain words only (numbers JSON unchanged). README section
  map, `code/paper_numbers.py` check labels and docstrings,
  `code/reproduce_p5.py` and `audits/README.md` renumbered to the paper's
  tables and sections (homology Table 1, defect scope Table 2, coverage
  Table 3); no check logic changed, 208/208 still pass.
- **Revision 6 (2026-10-09).** `paper/` replaced by the source's revision 6
  (git `6c899fd9`), made after a two-round external referee review: methods
  describe the corrected mapping as implemented, the post-2022 census's
  chain-A scoring and unreleased selection step are stated, census flow rows
  added to Table 2, coverage and exposure claims narrowed.
  `code/paper_numbers.py`: new section `rev6()` with 12 checks for the numbers
  this revision adds (220 checks: 157 recomputed, 43 reconciled, 20 re-read).
  README: check counts, the unchecked-items list and a note on the
  unreleased post-2022 selection step.
- **Revision 7 (2026-10-09).** `paper/` replaced by the source's revision 7
  (git `b5fa9f00`): new Section 6.4 measures homology exposure in two published
  evaluations (registered analysis e428; focused referee review, two rounds).
  Added from the source (git `3d1e5ebb`, results doc at `b5fa9f00`): the e428
  registration and proof, code, independent audit and its script, results
  write-up, per-chain classes, statistics and request receipts, the AlphaFlow
  split files (MIT, with the upstream `LICENSE` fetched at commit `0408d7c8`),
  the RCSB record of the removed entry 7DRH, and the e428 response cache.
  `data/derived/results/e428/t_structures.json` is the source file with the
  `rmsd` and `plddt` fields removed (values from the Terwilliger et al.
  workbook, which is not redistributed). `code/paper_numbers.py`: new section
  `rev7()` with 16 checks (236 checks: 168 recomputed, 45 reconciled,
  23 re-read). README: e428 entries in the identifiers, check counts,
  repository map, withheld items, web sources and licences.

## Copied files

"State" compares the copied bytes with the source repository at the commit
above: `= HEAD` (tracked and identical), `not tracked` (file present in the
working tree but not under version control), or the commit from which an
earlier stamped version was recovered.

| File in this repository | Source path | SHA-256 | Bytes vs source | State |
|---|---|---|---|---|
| `code/audit/e427_phaseB/recompute_p5.py` | `results/e427/audit/phaseB_subagent/recompute_p5.py` | `e4e09d7b722c1e6d571452363e97fcd824727878eb70eceb093fb6c1043c8d0a` | identical | = HEAD |
| `code/audit/e427_phaseB/sensitivity_p5.py` | `results/e427/audit/phaseB_subagent/sensitivity_p5.py` | `902bba852442331d1b4bda5759b925428d8c95b7de288d2522a19d981fcd863a` | identical | = HEAD |
| `code/audit/independent_label_sample.py` | `results/e427/audit/independent_label_sample.py` | `6b40a958bb83b4569b55cbe6f97a04859a8118367f682d2ed35eaae0088211e4` | identical | = HEAD |
| `code/audit/k7pq54_mapping_audit.py` | `paper/arxiv/verification/audit/k7pq54_mapping_audit.py` | `914faaae8bd9d9f1ff3b8d03cccd9fafe8b8fe47060a782ede247f363fe6ea0c` | identical | = HEAD |
| `code/experiments/e420_census.py` | `experiments/e420_census.py` | `10d8adaf7e8d308070f39b816501ccebb979ab7d2dc76e1c4d6e9c43b05c98e9` | identical | = HEAD |
| `code/experiments/e420_grammar_validate.py` | `experiments/e420_grammar_validate.py` | `ac3a02c0542cc7dbda280b54418f0fff83c77389259d17a47f7809f165a7501d` | identical | = HEAD |
| `code/experiments/e420_label_run.py` | `experiments/e420_label_run.py` | `13fbadd69b9ea8d6fec0d897c63dc0cea677b860745c952acb8fd0bcb5d792cf` | identical | = HEAD |
| `code/experiments/e420_labels.py` | `experiments/e420_labels.py` | `a9b82c712b3d5e1b6d731e8a1ae593c3c5bc3d1a63571df4c86d1e77e6227dae` | identical | = HEAD |
| `code/experiments/e420_overflow_check.py` | `experiments/e420_overflow_check.py` | `e920e7449aa7b22b6625eecc17e4786822f8b3c2be1dbbb3f1a45c583c173b89` | identical | = HEAD |
| `code/experiments/e420_pool.py` | `experiments/e420_pool.py` | `a6f9079ce3e5060a36bf2f33db33f2dabcc1da1b6d8f3c2005a4381974139bc4` | identical | = HEAD |
| `code/experiments/e421_census_fixed.py` | `experiments/e421_census_fixed.py` | `6e0db8d6831ab5d484545f3ceb8deed86e1f5f9e09288bf2c904abca897be392` | identical | = HEAD |
| `code/experiments/e421_conformal.py` | `experiments/e421_conformal.py` | `628b411d915c6bcaee24375bff2ce129c2c3a0f59fed02e2b256a63e0e24d3eb` | identical | = HEAD |
| `code/experiments/e421_pipeline.py` | `experiments/e421_pipeline.py` | `f2b924bd4c877bfb4961e375e2e563cac943e2345913bd2eb2880151b9870145` | identical | = HEAD |
| `code/experiments/e422_dispatch.py` | `experiments/e422_dispatch.py` | `b95eef70514230663b24d4a625042755a0f826508bfc2e931ef6c8dabd99623f` | identical | = HEAD |
| `code/experiments/e422_evaldriver.py` | `experiments/e422_evaldriver.py` | `9307dd1a00bfb1b464e04d73a680e990abdec3c7d1c6c50eaa11ba8cab69179b` | identical | = HEAD |
| `code/experiments/e422_manifest.py` | `experiments/e422_manifest.py` | `6142d362abe7a1e1762c15b0247f9dd53f87bff5081598b2526bb87a2da72b16` | identical | = HEAD |
| `code/experiments/e422_protocol.py` | `experiments/e422_protocol.py` | `9b773b7bc1611d58255917b5f52faf836d6cc1ea5f808a59bab11435c845ec33` | identical | = HEAD |
| `code/experiments/e422_runner.py` | `experiments/e422_runner.py` | `ac7729a53013b83712521707776c340fbd24700190cd0b1b169511d138e0e523` | identical | = HEAD |
| `code/experiments/e422_status.py` | `experiments/e422_status.py` | `e10c529a73299a84986524655a6d9941fc9818facdd66a300508435402a9e96e` | identical | = HEAD |
| `code/experiments/e423_dispatch.py` | `experiments/e423_dispatch.py` | `296daef3ff1f014153615f3d00054bac448f3f309a2a5b38b6873ef2b6ac522e` | identical | = HEAD |
| `code/experiments/e423_runner.py` | `experiments/e423_runner.py` | `04ecfa47ead990093123894e51c17d0a93a9131a36eb81286895a6ddfc744a2f` | identical | = HEAD |
| `code/experiments/e424_evalgate.py` | `experiments/e424_evalgate.py` | `6c4ef7ba63bea03b1f763248512570d41645627eabb018ad78758c736edb63f3` | identical | = HEAD |
| `code/experiments/e425_execreceipt.py` | `experiments/e425_execreceipt.py` | `d0bee3e1132edea75357a3af2040a30cfb997921acebb7e8eeeb08edc22fe208` | identical | = HEAD |
| `code/experiments/e426_attested_empty.py` | `experiments/e426_attested_empty.py` | `0f9dfe98df3e1851bb9004e995b91e79b5d045e4cadaa595b208a8ea6be727b9` | identical | = HEAD |
| `code/experiments/e427_e421_relabel.py` | `experiments/e427_e421_relabel.py` | `c4cf94b12ec53488bce1ad230532734f8288d25ed653d7fcaaca366d7e119316` | identical | = HEAD |
| `code/experiments/e427_label_repair.py` | `experiments/e427_label_repair.py` | `8f7328e6d719e3b669e8df91724add7b69f28853ef80da367cdf7be758b6d6f6` | identical | = HEAD |
| `code/experiments/e427_offline_relabel.py` | `experiments/e427_offline_relabel.py` | `a282dd0110645524e27adeaa5309aad6977a14784caa55a07f0ed34cd6ed4da7` | identical | = HEAD |
| `code/experiments/e427_pause.py` | `experiments/e427_pause.py` | `53ec0568a47d92f57653e13a1295a15bf7e2b7508828375b2d91fda3b184137e` | identical | = HEAD |
| `code/experiments/e427_replay_relabel.py` | `experiments/e427_replay_relabel.py` | `e747439f0ea4437ea2eab1a0a400c7d253a70b6abefe50de2e8995ae90c96a63` | identical | = HEAD |
| `code/experiments/e427_rerun_p5.py` | `experiments/e427_rerun_p5.py` | `7965e129c0cac7fc21e2dfa4ac2c3c569e46eeebd82e8c365f0a3d9b39754635` | identical | = HEAD |
| `code/tests/conftest.py` | `tests/conftest.py` | `00277a3a00b858de173b72af4443680c0ea7e0ce2792a443b843970a6e731450` | identical | = HEAD |
| `code/tests/test_e422_dispatch.py` | `tests/test_e422_dispatch.py` | `1c793de6e3ea4160e9819f29529e37d7328115ae31b04712a8b1bb2f357fbc7b` | identical | = HEAD |
| `code/tests/test_e422_evaldriver.py` | `tests/test_e422_evaldriver.py` | `9e609a45f9c1a8aa501193a79315b8fa26e32741f6d1b94c2244d6d83bb8fb93` | identical | = HEAD |
| `code/tests/test_e422_manifest.py` | `tests/test_e422_manifest.py` | `969bb175f1cd0bd1720ee690d876592199bf9a7df5238491bb2c7fd269789381` | identical | = HEAD |
| `code/tests/test_e422_protocol.py` | `tests/test_e422_protocol.py` | `a97561f944673a4331becc734a6bbd969773a7e6eef5ea2f2efe5c7f1c67a57b` | identical | = HEAD |
| `code/tests/test_e422_runner.py` | `tests/test_e422_runner.py` | `aa551b4e3f71d4be8e7d0d7ee3680a8aac2a76cef196c850c919be2de8c6562c` | identical | = HEAD |
| `code/tests/test_e422_status.py` | `tests/test_e422_status.py` | `fa064e76326a356bb0f6d9d7528f86012b963a36e0ebf27c9efe87e4654bc082` | identical | = HEAD |
| `code/tests/test_e423_dispatch.py` | `tests/test_e423_dispatch.py` | `60eeef0a4e929f1976c8d11dff4e5370370e38f26dffa6d4a0bdc278ee9e953d` | identical | = HEAD |
| `code/tests/test_e423_runner.py` | `tests/test_e423_runner.py` | `d98091e6e7efb6fb0d2e16b1d6a9f3081fc3e5711df2c9d8cb5ca930d811c4d5` | identical | = HEAD |
| `code/tests/test_e424_evalgate.py` | `tests/test_e424_evalgate.py` | `ed4dd778f49972684388a87c067666a87d4d2a7ac567187200f80e9906e79543` | identical | = HEAD |
| `code/tests/test_e425_execreceipt.py` | `tests/test_e425_execreceipt.py` | `6c8d77089bf6bc96df46b0734ef357931988d8dacfa98e802379e1f8a70d5f1f` | identical | = HEAD |
| `code/tests/test_e426_attested_empty.py` | `tests/test_e426_attested_empty.py` | `72f0f0517a0538259992895d8ca0bb83a09a63d7a457f58ca5e3a31818945b8b` | identical | = HEAD |
| `code/tests/test_e427_label_repair.py` | `tests/test_e427_label_repair.py` | `addbddc45e390b698eefbe20c1bcac69bf0d1bf0bf7a298bf14eb1456b96d382` | identical | = HEAD |
| `code/verification/make_figures_rev3.py` | `paper/arxiv/verification/make_figures_rev3.py` | `651001efc3932396a90388f311883adc5459313b5a8137e4d1517d0569e0075f` | identical | = HEAD |
| `code/verification/make_figures_rev4.py` | `paper/arxiv/verification/make_figures_rev4.py` | `0f9e6c2a408002b3b1f37a96bfe7860cc04e1fac2e409a458c9ada8c0b6c6522` | identical | git `52345ba7` |
| `code/verification/v1_channel_recount.py` | `paper/arxiv/verification/v1_channel_recount.py` | `06ee861947492e424bb7863d1c0a148a0c1c074b0490ac4f3702772c08555655` | identical | = HEAD |
| `code/verification/v2_prior_exposure.py` | `paper/arxiv/verification/v2_prior_exposure.py` | `4c8359e422a4be33e4051c452784b39ec603fd775e886929afe2477cfc25d33f` | identical | = HEAD |
| `code/verification/v2c_coverage_extension.py` | `paper/arxiv/verification/v2c_coverage_extension.py` | `3a44267d22ee28f84d71a526365400d0da279d2dab22a8913898e655390de3d5` | identical | = HEAD |
| `code/verification/v3_stratified.py` | `paper/arxiv/verification/v3_stratified.py` | `48e8d285a3e5fe20a75a10b01245536bc5ad3218918b6d1018801bd6567bb6ab` | identical | = HEAD |
| `code/verification/v4_k7pq54_diagnostics.py` | `paper/arxiv/verification/v4_k7pq54_diagnostics.py` | `430d9634613c9e60afe784e89a3666c3e67d24380dc22042410ad4e2c08288bb` | identical | = HEAD |
| `code/verification/v5_templates.py` | `paper/arxiv/verification/v5_templates.py` | `0fb91482255d3a045e9fd53028d3eb838c12050a15bb52e8663ca9cb910481d0` | identical | = HEAD |
| `code/verification/v6_exposure_effect.py` | `paper/arxiv/verification/v6_exposure_effect.py` | `8818968d9cb034d03926913e376c94a74dab5792dec38172e149691475116b59` | identical | = HEAD |
| `code/verification/v7_trigger_prevalence.py` | `paper/arxiv/verification/v7_trigger_prevalence.py` | `a1ca79ed5c0656a6107c01c9d1707ff0ad48aa1ed887d2fec3cac504c053ec54` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/audit/k7pq54_9qj6_audit.json` | `paper/arxiv/verification/audit/k7pq54_9qj6_audit.json` | `039cdbd0d95de37cad63faa9a5f43a9033bde27b8b7ad38413eafdb8a3dbce2b` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/out/figures_rev3_numbers.json` | `paper/arxiv/verification/out/figures_rev3_numbers.json` | `43ee0bc7f3056a2448e321d3fbc8f3c02af38b1793274ba4779ee574cb0e8246` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/out/v1_channel_recount.json` | `paper/arxiv/verification/out/v1_channel_recount.json` | `1f7e0d323310e49aa3db53ec4edcd3cad46adc6d3b825d3ee54ccc04a612b044` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/out/v2_prior_exposure.json` | `paper/arxiv/verification/out/v2_prior_exposure.json` | `19fcabc7a36bc73f02fa5be0f963e53cc84cdc702c80285c7ad7fe88878ce720` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/out/v2c_extension.json` | `paper/arxiv/verification/out/v2c_extension.json` | `cd22c8ac3f13bc9ee15b0b9d186d75fe354600ccc503ed7b7ed5a359f71e4051` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/out/v3_stratified.json` | `paper/arxiv/verification/out/v3_stratified.json` | `e6d734c1b2bcb26d9c9845e1818b450efab42a1719b767d2b365563e65ea7d03` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/out/v4_k7pq54.json` | `paper/arxiv/verification/out/v4_k7pq54.json` | `54694397335fc6b4a18492f0fb506270b5f821e0fa4f35f5508d25776e2ade16` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/out/v5_templates.json` | `paper/arxiv/verification/out/v5_templates.json` | `d8b184c902819f2f7969ea19778b5bf3496438eb9de915b88f294f1c6ef7e69a` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/out/v6_exposure_effect.json` | `paper/arxiv/verification/out/v6_exposure_effect.json` | `d0c0f89cdf0271757480424feb67d748e374d283def1a894f375cbef03075cf1` | identical | = HEAD |
| `data/derived/paper/arxiv/verification/out/v7_trigger_prevalence.json` | `paper/arxiv/verification/out/v7_trigger_prevalence.json` | `37efe7281f13173493432aaf10cba595aa8c7b255c5d519ae856d9b265bf1e16` | identical | = HEAD |
| `data/derived/results/e420/label_run_receipt.json` | `results/e420/label_run_receipt.json` | `b2c36992ec965b77a0d040e44c25fcebbd1ff1fea03a7e3e9a1747d6f593a452` | identical | = HEAD |
| `data/derived/results/e420/overflow_set.json` | `results/e420/overflow_set.json` | `3fcec8e55949f335c7d4bd5374aeaf3bc981fecbf652123507b1ff10f72773fc` | identical | not tracked |
| `data/derived/results/e420/overflow_temporal_check.json` | `results/e420/overflow_temporal_check.json` | `6c8cf0a26757a252d44c500ca1619812d1c6349a0bf1b24313ad685bf06f95d7` | identical | = HEAD |
| `data/derived/results/e420/source_availability_manifest.json` | `results/e420/source_availability_manifest.json` | `8de716b8bd587cab2e8e65416e7b45a94a39311c085b5da0259d97260fe2b144` | identical | = HEAD |
| `data/derived/results/e421/census_failures.json` | `results/e421/census_failures.json` | `d9b912d9fdfe3142286751e4a1434c7e23ddde4b43e78259fd044dc5a15a7933` | identical | = HEAD |
| `data/derived/results/e421/census_results.json` | `results/e421/census_results.json` | `df5a6833734c8173dfc6c5cc9c452ecfc107701acfde66a189714b7dfaa182f1` | identical | = HEAD |
| `data/derived/results/e421/labeled_residues_with_plddt.json` | `results/e421/labeled_residues_with_plddt.json` | `9a771761d668b7250891602b91223787250f02a47ef1edb107fc674ce4626dfd` | identical | = HEAD |
| `data/derived/results/e421/qualifying_entries.json` | `results/e421/qualifying_entries.json` | `96f5a4ea2e735551c59f64e105c5d3154ad20ae424d264838e0db7ebe952cb6d` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch01_checkpoint.json` | `results/e422/batches/e422_batch01_checkpoint.json` | `526fd4c54603ae128a2b6313a09c10c7a18e5f5699af3a3d7fc46bfb81c5d029` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch01_checkpoint.json.ots` | `results/e422/batches/e422_batch01_checkpoint.json.ots` | `fbd98770594abe466ea28023106cc2279556df721aeb7a7c6caccf0271e05bab` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch01_freeze.json` | `results/e422/batches/e422_batch01_freeze.json` | `bf46fee816e41ab68457fc06a84eb96a6abccd87d940a7ca7a40ccec3f028011` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch02_checkpoint.json` | `results/e422/batches/e422_batch02_checkpoint.json` | `08cc71203d3830813f5b096e27a73e6d1e53cfc1776b20a0a04cb95c772a4aea` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch02_checkpoint.json.ots` | `results/e422/batches/e422_batch02_checkpoint.json.ots` | `faa4c191186c9c8494d60c6e3d096c06ddc399969c4ffa1cb6cd95373c3583b7` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch02_freeze.json` | `results/e422/batches/e422_batch02_freeze.json` | `c96a2b8b88bc2fde46375227b71a02ddb344bf5797242abcc83f59429d52c03b` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch03_checkpoint.json` | `results/e422/batches/e422_batch03_checkpoint.json` | `b31e0ad5a00ea6e12c8c809d51f59e5eba75966f3dd8f38f42b61f51b6d61c3c` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch03_checkpoint.json.ots` | `results/e422/batches/e422_batch03_checkpoint.json.ots` | `fe93741165e72011882560cc549675a62e70688abd843538d95bbab7ecb9ea49` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch03_freeze.json` | `results/e422/batches/e422_batch03_freeze.json` | `654e5b78a5d66290239cef2957fa7b02dca915db25acad3a0a42c28ea50abf48` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch04_checkpoint.json` | `results/e422/batches/e422_batch04_checkpoint.json` | `a08e538db7dd54f2a046f2bbc7ed3c9f37e6c2e1ba2c30f42403aa66f6fa934b` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch04_checkpoint.json.ots` | `results/e422/batches/e422_batch04_checkpoint.json.ots` | `0adfe6bd4c5ef88d54aad88e54896645b990ae1e0511210342dad23612a4e774` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_batch04_freeze.json` | `results/e422/batches/e422_batch04_freeze.json` | `5dbe67e1a15678e22ffb5545e2f52251c140c7e022de6002f88d1628c9f5775b` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_binding_state.json` | `results/e422/batches/e422_binding_state.json` | `709ac7f5e016dc1f07343f5fabe9a9ee559e26187253b0bef43d4c342382f63d` | identical | = HEAD |
| `data/derived/results/e422/batches/e422_binding_state.json.ots` | `results/e422/batches/e422_binding_state.json.ots` | `8624e49770a33583989fd6af4cc8f8842816889116fe2423b1891e26249decfd` | identical | = HEAD |
| `data/derived/results/e422/e426_batch1_empty_attestation.json` | `results/e422/e426_batch1_empty_attestation.json` | `7f470270d61c012699b1135889229f0adab8f3da2bef4cfdfcc638cb3ca02a92` | identical | = HEAD |
| `data/derived/results/e422/e426_batch1_empty_attestation.json.ots` | `results/e422/e426_batch1_empty_attestation.json.ots` | `3dc9f3d1c1f919f6b7c3c5c46e1636501bae4386dfe4cc84e1fbc44608fc65ac` | identical | = HEAD |
| `data/derived/results/e422/evals/e422_batch04_evaluation.json.ots` | `results/e422/evals/e422_batch04_evaluation.json.ots` | `6cc5c2287748ec24fc5d36fac6488337835bff1ebf0122ba2c2bf239ceadf3d7` | identical | = HEAD |
| `data/derived/results/e422/evals/e422_batch04_execution_receipt.json` | `results/e422/evals/e422_batch04_execution_receipt.json` | `5ffcf544bd146ecda19bbd2f88980ce7a7a485768f82db3c30f0cda9f84589b8` | identical | = HEAD |
| `data/derived/results/e422/evals/e422_batch04_execution_receipt.json.ots` | `results/e422/evals/e422_batch04_execution_receipt.json.ots` | `06d4c9c9514951f674d3140c07ecaae10e41bed40d49ee2b56ee200772d6e041` | identical | = HEAD |
| `data/derived/results/e422/manifest_checkpoint.json` | `results/e422/manifest_checkpoint.json` | `1513e43ce3f2178da8da72650b92f178390ff10b64f9b518e92c0b8487cff3e9` | identical | = HEAD |
| `data/derived/results/e422/manifest_checkpoint.stamped_440566b4.json` | `results/e422/manifest_checkpoint.json` | `07853dd9d2d247d3d5e6ff71b37622585efc59135930b4bfebedb2f0a6983102` | identical | git `440566b4` (OTS-stamped version) |
| `data/derived/results/e422/manifest_checkpoint.stamped_440566b4.json.ots` | `results/e422/manifest_checkpoint.json.ots` | `2794bf15d745899cf25ff618ac44133a91398d951e89a1458972bc71f8ed6f5e` | identical | = HEAD |
| `data/derived/results/e427/audit/independent_label_sample.json` | `results/e427/audit/independent_label_sample.json` | `151cb01ad6137ca03e779ad910e57cbd769cfee0d0340439e48fa0cfd45ce215` | identical | = HEAD |
| `data/derived/results/e427/audit/phaseB_subagent/recompute_p5.json` | `results/e427/audit/phaseB_subagent/recompute_p5.json` | `9613b8b19af18e77db58d04423c718d23461221cf5a42d57227e2750adf42d4d` | identical | = HEAD |
| `data/derived/results/e427/audit/phaseB_subagent/sensitivity_p5.json` | `results/e427/audit/phaseB_subagent/sensitivity_p5.json` | `875ad12af53c122507c1f3d159e5d919cf91412c55b825937f262141faab52f2` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_batch01_checkpoint.json` | `results/e427/batches/e422_batch01_checkpoint.json` | `526fd4c54603ae128a2b6313a09c10c7a18e5f5699af3a3d7fc46bfb81c5d029` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_batch01_freeze.json` | `results/e427/batches/e422_batch01_freeze.json` | `bf46fee816e41ab68457fc06a84eb96a6abccd87d940a7ca7a40ccec3f028011` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_batch02_checkpoint.json` | `results/e427/batches/e422_batch02_checkpoint.json` | `98b5ccb8af4e7899f665d703d1f55f75cda3d2387391a92199e9817881548935` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_batch02_checkpoint.json.ots` | `results/e427/batches/e422_batch02_checkpoint.json.ots` | `fe768b09d8e92e1a138f24372409216c47791b4a9e7314d2026972a7b5390470` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427/batches/e422_batch02_freeze.json` | `results/e427/batches/e422_batch02_freeze.json` | `c96a2b8b88bc2fde46375227b71a02ddb344bf5797242abcc83f59429d52c03b` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_batch03_checkpoint.json` | `results/e427/batches/e422_batch03_checkpoint.json` | `7153a7f69b7cbd3d4c894b9a9638ef3e28d9c7dea3979960fc0b4bd189f8c8d3` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_batch03_checkpoint.json.ots` | `results/e427/batches/e422_batch03_checkpoint.json.ots` | `d043a907b64be669ed977c188aa6d99946f02638f9d1da582853617cdb557a90` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427/batches/e422_batch03_freeze.json` | `results/e427/batches/e422_batch03_freeze.json` | `654e5b78a5d66290239cef2957fa7b02dca915db25acad3a0a42c28ea50abf48` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_batch04_checkpoint.json` | `results/e427/batches/e422_batch04_checkpoint.json` | `2a2810a02d9b5a14ac70bebe5a92462353da81b229b875971ff68ec7720330d8` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_batch04_checkpoint.json.ots` | `results/e427/batches/e422_batch04_checkpoint.json.ots` | `31e1637d9a4b1e12d150f43eae074bb09a1e687f8f0cc3d029c933c466ff5a9b` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427/batches/e422_batch04_freeze.json` | `results/e427/batches/e422_batch04_freeze.json` | `5dbe67e1a15678e22ffb5545e2f52251c140c7e022de6002f88d1628c9f5775b` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_binding_state.json` | `results/e427/batches/e422_binding_state.json` | `f023cdbd2c5cf351cdde05f9e9e690a9115b0806f362338d209e34db99189c6f` | identical | = HEAD |
| `data/derived/results/e427/batches/e422_binding_state.json.ots` | `results/e427/batches/e422_binding_state.json.ots` | `75740a8616f8cb329739660809bc2b31ae96444818fff2bc875a0cfe22ab161d` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427/e421_relabel/census_results.json` | `results/e427/e421_relabel/census_results.json` | `7cec7709b0fd97d44193835e72ae6abe4d36a6886dd9bf0d619889c21c1d47db` | identical | = HEAD |
| `data/derived/results/e427/e421_relabel/census_results.json.ots` | `results/e427/e421_relabel/census_results.json.ots` | `a3fa82eb27f28741ca57c0ae3da66f446acbe03f1b4887c788faafe60670a917` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427/e421_relabel/receipts.json` | `results/e427/e421_relabel/receipts.json` | `253545ebf225190f16c8eddb188be58ebc028f32c3e569a2593942a03f9e2990` | identical | = HEAD |
| `data/derived/results/e427/e421_relabel/receipts.json.ots` | `results/e427/e421_relabel/receipts.json.ots` | `39c47fbc6e18fd3260c1d761ca69893d7f57296b30fe94b736b229d1b272bab5` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427/e421_relabel/summary.json` | `results/e427/e421_relabel/summary.json` | `4aac0678f21057391597319675a6950111f4c1bf89eb60c87546154e4a06d511` | identical | = HEAD |
| `data/derived/results/e427/e421_relabel/summary.json.ots` | `results/e427/e421_relabel/summary.json.ots` | `2721f8ce54b80589e3a63ef7b675297b22703428079242c0fc630fcc8ae48297` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427/evals/e422_batch04_evaluation.json.ots` | `results/e427/evals/e422_batch04_evaluation.json.ots` | `1e1470b3bb2e5223ed8c9fbd3c1f72a14006c690adaa397be0af6d8dbf678c0d` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427/evals/e422_batch04_execution_receipt.json` | `results/e427/evals/e422_batch04_execution_receipt.json` | `766f4c1f8a72ed2e89d275aea17ef22e49ec4419e3c137c3db4bb1e54587d198` | identical | = HEAD |
| `data/derived/results/e427/evals/e422_batch04_execution_receipt.json.ots` | `results/e427/evals/e422_batch04_execution_receipt.json.ots` | `7eccd0a6a7aca34dcdd40ee1d29029f6962c214a16223edb9e8fc2dd51481b0b` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427/replay_report.json` | `results/e427/replay_report.json` | `f9b540026b35abc9d271b7121d8017259f58f0d44148d44aae9c36dd15f779bc` | identical | = HEAD |
| `data/derived/results/e427/replay_report.json.ots` | `results/e427/replay_report.json.ots` | `3dc53768674d73ae604d1b98f93602287ae9a162d55450233c9081f62a18a23e` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `data/derived/results/e427_repair/offline_relabel_comparison.json` | `results/e427_repair/offline_relabel_comparison.json` | `ff969c4f37115f6d5325cfe539d1a2f59551d979ee900622f296b1886e24186a` | identical | = HEAD |
| `data/derived/results/e427_repair/repaired_rows_batches01-04.json.gz` | `results/e427_repair/repaired_rows_batches01-04.json.gz` | `712dedc48e49ddd2bb27c28231ff107e84938d69febb96bdf761e1801034f1fd` | identical | = HEAD |
| `paper/figures/rev3/f1_mechanism.pdf` | `paper/arxiv/figures/rev3/f1_mechanism.pdf` | `e3d67e95069094ba1463b1b587ada0459d3183181f11f9ec4ed7ac47e9092139` | identical | = HEAD |
| `paper/figures/rev3/f2_scope.pdf` | `paper/arxiv/figures/rev3/f2_scope.pdf` | `5c9cf2552ecb46a42beb670c60a27c721fe70e897fbfab35a2a22083ea4ea305` | identical | = HEAD |
| `paper/figures/rev3/f3_exposure.pdf` | `paper/arxiv/figures/rev3/f3_exposure.pdf` | `75d3a20a57a16c36f87d49019068b3b71b86c0fbd245f2e2eb105479f87dd403` | identical | = HEAD |
| `paper/figures/rev4/f1_mechanism.pdf` | `paper/arxiv/figures/rev4/f1_mechanism.pdf` | `c6d6d49a1d02e2a67198a25afe31c9e807ea5fff7b45937e4c110d44730ffbac` | identical | git `52345ba7` |
| `paper/figures/rev4/f2_scope.pdf` | `paper/arxiv/figures/rev4/f2_scope.pdf` | `f12520c1a2a87cfcdddbf6375c8bbf61cca6648d26a688420592c5f80fef95a0` | identical | not tracked |
| `paper/references.bib` | `paper/arxiv/references.bib` | `3cf2d521d1ddfa00b3258892bd0fc88220cf9ca9984332faac9db647dd3ee940` | identical | git `52345ba7` |
| `paper/temporal_leakage.tex` | `paper/arxiv/temporal_leakage.tex` | `e2f86213677daf3b9b152962e88d299c28cb3c5e6f9968954dd34509a6b03c85` | identical | git `6c899fd9` |
| `paper/temporal_leakage_paper.pdf` | `paper/arxiv/temporal_leakage_paper.pdf` | `efdd6b21cc796d88e8ba7a1768d93a34dece78720907c0927fe824593aaf8404` | identical (built with tectonic) | git `6c899fd9` |
| `registrations/e420_registration.md` | `docs/e420_registration.md` | `ea6d8725e23f504607e75a869ae33c3e189bcd28f6533a2f739b388462789d99` | identical | = HEAD |
| `registrations/e421_post_run_checkpoint.json` | `docs/e421_post_run_checkpoint.json` | `a0ebfa59719a19c0af32a6a40eb628435b37ec44f2809a31f99d7249243a13a1` | identical | = HEAD |
| `registrations/e421_post_run_checkpoint.json.ots` | `docs/e421_post_run_checkpoint.json.ots` | `7fac3db8adaff048052281a4f938d3d483dd9210db5ddd96bb2a7d5cc1679764` | identical | = HEAD |
| `registrations/e421_post_run_checkpoint_v2.json` | `docs/e421_post_run_checkpoint_v2.json` | `2201ecf6f02c39b7c75ff8984d1e02771c0cf859c24d02123053d23f0c924fb4` | identical | = HEAD |
| `registrations/e421_post_run_checkpoint_v2.json.ots` | `docs/e421_post_run_checkpoint_v2.json.ots` | `9b78ef03b21ed019fe578c20f843a8d4434186aa85df917873a8fa5ba4476f01` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `registrations/e421_post_run_checkpoint_v3.json` | `docs/e421_post_run_checkpoint_v3.json` | `0d34ba08bdbba075d367e77b876fd0c69823f7d3dc7d5cf5c22f3bf0214157b7` | identical | = HEAD |
| `registrations/e421_post_run_checkpoint_v3.json.ots` | `docs/e421_post_run_checkpoint_v3.json.ots` | `d4870af2dfb4f40f52dfd8a1e6ee5e691d926401170225fa963f9d96ce12cc5a` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `registrations/e421_pre_source_checkpoint.json` | `docs/e421_pre_source_checkpoint.json` | `afc0cae93d118a8616d205720d5d9cd0c54e8685dc1be743bd9ed37525687f8e` | identical | = HEAD |
| `registrations/e421_pre_source_checkpoint.json.ots` | `docs/e421_pre_source_checkpoint.json.ots` | `281d3399fc1244ed87873aaf77365160e235b7aacc9bc4cfce7321410cae5e30` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `registrations/e421_registration.md` | `docs/e421_registration.md` | `c298bef55fe21c4d73bce93266de1166540f5c77e0788fe5fc9cfb7111a1cd7c` | identical | = HEAD |
| `registrations/e422_anchor_receipt_gate_a_20260909.json` | `docs/e422_anchor_receipt_gate_a_20260909.json` | `0174948fbe40f273d5a98a895289a2476a678407122352fed7fa4b9bd23d44cc` | identical | = HEAD |
| `registrations/e422_batch_runbook.md` | `docs/e422_batch_runbook.md` | `c0cdb9b7dff3ee7ed567d656a41606f4914b10cd0ca144dffcd8d3777b7adee3` | identical | = HEAD |
| `registrations/e422_correction_of_record_20261007.md` | `docs/e422_correction_of_record_20261007.md` | `e8fe3fa520738cd9929147a7b7cd7ae6be943da9699c83dc7bd99167363cac7b` | identical | = HEAD |
| `registrations/e422_correction_of_record_20261007.md.ots` | `docs/e422_correction_of_record_20261007.md.ots` | `7bef47c63746d8be5439cb654688ea6c300774e9e720c16e93925da018b1c09c` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `registrations/e422_freeze_gate_a_20260909.md` | `docs/e422_freeze_gate_a_20260909.md` | `6220f7bc8f199b32ea1bb484693bacba5dc69824944cfec54c850f829e05b383` | identical | = HEAD |
| `registrations/e422_freeze_gate_a_20260909.md.ots` | `docs/e422_freeze_gate_a_20260909.md.ots` | `388bc60fe41bcf814ed00e8df8bc4b1eb0e4822f9f7bc606e771f64a0314315a` | identical | = HEAD |
| `registrations/e422_registration.md` | `docs/e422_registration.md` | `dc7a362773fbf3a03e759aa7d178125963a14d88d9bbd61b66389a04ba8d7b3b` | identical | = HEAD |
| `registrations/e422_registration.md.ots` | `docs/e422_registration.md.ots` | `35f0cdc61b3cbbbf2616a6c208308c9a36dd7dcf4bffcc4bfd39027aba017dbe` | identical | = HEAD |
| `registrations/e423_dispatch_repair_record_20260910.md` | `docs/e423_dispatch_repair_record_20260910.md` | `6b7acb4b879de6f3a0b01fd9d6b89d6d2a6ce09112d28797d426b0039231ce76` | identical | = HEAD |
| `registrations/e424_binding_gate_record_20260910.md` | `docs/e424_binding_gate_record_20260910.md` | `6b9696dd58a71f9146cfac9ae7a6b5de471611cb37b73c40ad4e175edbf4b3bd` | identical | = HEAD |
| `registrations/e425_execution_receipts_record_20260910.md` | `docs/e425_execution_receipts_record_20260910.md` | `9dceab94cb80afec03a95bc385b53b4a09a6933b486166c66b0dad96a6061544` | identical | = HEAD |
| `registrations/e425_execution_receipts_record_20260910.md.ots` | `docs/e425_execution_receipts_record_20260910.md.ots` | `5d1b6866daeaf814e6406b7d902b65df03a701eef2e1aa5781242a95564fb038` | identical | = HEAD |
| `registrations/e427_registration.md` | `docs/e427_registration.md` | `dbec42fbe1db1dc18592d651b988755fc6e02a7b883c61ea0f5f9030ab2f4fd5` | identical | = HEAD |
| `registrations/e427_registration.md.ots` | `docs/e427_registration.md.ots` | `90aa7c8c9f531eb8410cf222d3f84a82fa93f9836f6483084082a23ba86a1b5c` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `registrations/stamped_versions/e422_correction_of_record_20261007_add20261008.md` | `docs/e422_correction_of_record_20261007.md` | `faf68a06e4fc13988e0861898b20c7ac2dd7670ed42cf29846e41ee3f76dc863` | identical | git `a950b9fc` (OTS-stamped version) |
| `registrations/stamped_versions/e422_correction_of_record_20261007_add20261008.md.ots` | `docs/e422_correction_of_record_20261007_add20261008.md.ots` | `7771ee05d7f8745030feb9513e3b2e2bb144131d9ab70ee06bf1586aa8e29f50` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `registrations/stamped_versions/e422_correction_of_record_20261007_add20261009a.md` | `docs/e422_correction_of_record_20261007.md` | `f3675e3dc9504c7290be16b2f50869b63735b7a80e59479bc651127e6b1e0de1` | identical | git `bf55c601` (OTS-stamped version) |
| `registrations/stamped_versions/e422_correction_of_record_20261007_add20261009a.md.ots` | `docs/e422_correction_of_record_20261007_add20261009a.md.ots` | `088d1d2e08dd40a1ba7bac980a2105ceaa1b195fc3b1fc3e07cff9929d26bbb3` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `registrations/stamped_versions/e423_dispatch_repair_record_20260910.stamped_4e31778c.md` | `docs/e423_dispatch_repair_record_20260910.md` | `e4967215dd60c53db30e7d2612be7c3e21e5d9f7bfea65d8f758bff0eb416d4c` | identical | git `4e31778c` (OTS-stamped version) |
| `registrations/stamped_versions/e423_dispatch_repair_record_20260910.stamped_4e31778c.md.ots` | `docs/e423_dispatch_repair_record_20260910.md.ots` | `ad8aa09269273430c582b2aabfdbac7f69cd00d4f73963cc3e4c50ba5dabd1d6` | identical | = HEAD |
| `registrations/stamped_versions/e424_binding_gate_record_20260910.stamped_88310235.md` | `docs/e424_binding_gate_record_20260910.md` | `b17189f431fe99a799da95d34e8317cbb8737dacc1ba0d41f7b5d43a7bf31bf5` | identical | git `88310235` (OTS-stamped version) |
| `registrations/stamped_versions/e424_binding_gate_record_20260910.stamped_88310235.md.ots` | `docs/e424_binding_gate_record_20260910.md.ots` | `354139c448c25f1fee3f70b5a5513c08ad0afc8c5ebbeead5a4844bcd85f9934` | identical | = HEAD |
| `registrations/verification/REGISTRATION.md` | `paper/arxiv/verification/REGISTRATION.md` | `9c88e1c14c88fa8a2f52596dab5bec9217101814a1649f832f738c52fc889857` | identical | = HEAD |
| `registrations/verification/REGISTRATION_A3_24a51a94.md` | `paper/arxiv/verification/REGISTRATION.md` | `b0b917a9638dd88c7ef41b9322ff1d120c44ac0fd640b7e1ade60f46050e6008` | identical | git `24a51a94` (OTS-stamped version) |
| `registrations/verification/REGISTRATION_A3_24a51a94.md.ots` | `paper/arxiv/verification/REGISTRATION_A3_24a51a94.md.ots` | `a6d556ed01ce23c60e4e9c01afbb7436d7ededa383bdb1757359f8e39a979c92` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |
| `registrations/verification/REGISTRATION_A4_de0e7a32.md` | `paper/arxiv/verification/REGISTRATION.md` | `74cc5059eeff4340143f24ee9a5e21565004a523cd36bb06237b62017a5f3969` | identical | git `de0e7a32` (OTS-stamped version) |
| `registrations/verification/REGISTRATION_A4_de0e7a32.md.ots` | `paper/arxiv/verification/REGISTRATION.md.ots` | `4d67f7ba90b2c26b469140f42bd8c3c96ba6538e326811f2633580a43fd03b17` | identical | = `c261d0fe` (upgraded proof; first-release SHA-256 below) |

| `registrations/e428_registration.md` | `docs/e428_registration.md` | `2a2cc1e7b2c59efefc6a57a4cfdea29ea47a812f44474c5fe6fdd67110bf8661` | identical | = HEAD |
| `registrations/e428_registration.md.ots` | `docs/e428_registration.md.ots` | `8ba02bdb5ed7353813f1fa816dd606cd536c5ffd34fe7231718a29060372bd3f` | identical | = HEAD |
| `code/experiments/e428_exposure_audit.py` | `experiments/e428_exposure_audit.py` | `82965c631a63928a5f7cc4b7e9a48e8a1d57ba2745dd28f8a744acd27409bb72` | identical | = HEAD |
| `audits/e428/audit_e428.md` | `docs/audit_e428.md` | `ab2e2955c34b7df7766da968a9f7d4fb313f18b40d7565e22739f590505ed265` | identical | = HEAD |
| `audits/e428/audit_e428.py` | `experiments/audit_e428.py` | `d1e9f13b5b400b32f0e820abdf51d79c7f1c6a5cbf66d5e4756fd4cf55851a0e` | identical | = HEAD |
| `audits/e428/e428_results.md` | `docs/e428_results.md` | `0f2c3754e1e119842d6f93b9eb408d2b838856eb38635d7082edb9508f615d99` | identical | = HEAD |
| `data/derived/results/e428/a_chains.json` | `results/e428/a_chains.json` | `82706cd06f72aded34472f8c6249c00c7dfb5d68f47c3ab7f6aaf009def82289` | identical | = HEAD |
| `data/derived/results/e428/stats.json` | `results/e428/stats.json` | `1f1c3fbfb6e2e26d9ce5fa050352851f6250bd623da6e6d974d4002fce1f5b7e` | identical | = HEAD |
| `data/derived/results/e428/receipts.json` | `results/e428/receipts.json` | `9240acc21ed158489463445cd803fe68cf17e4a4dbaac077355e45b39b734a99` | identical | = HEAD |
| `data/derived/results/e428/t_structures.json` | `results/e428/t_structures.json` | `ebf29c8502086d52af7b3ad27c37c3aff1a3b6b7cc794b4a8ad67c6ba3d8157e` | differs (see Revision 7) | derived (rmsd, plddt removed) |
| `data/derived/data/e428/alphaflow/atlas_train.csv` | `data/e428/alphaflow/atlas_train.csv` | `f9618d3195a6ed48a42eeece38efcd6d4485ef07231a2a0b9886e7b5f488a4f4` | identical | = HEAD |
| `data/derived/data/e428/alphaflow/atlas_val.csv` | `data/e428/alphaflow/atlas_val.csv` | `a4399b4004c694c7760fdffcf3e2db983278bbfd40c2c72f3357ddc91c42383c` | identical | = HEAD |
| `data/derived/data/e428/alphaflow/atlas_test.csv` | `data/e428/alphaflow/atlas_test.csv` | `c5029649db0b7ff6ba41be6094341dde1ecd3036f07b011df4fe05ba2b305a9c` | identical | = HEAD |
| `data/derived/data/e428/alphaflow/SHA256SUMS` | `data/e428/alphaflow/SHA256SUMS` | `f7cfdbef1c670f1d4421bfe5f58860a0d6f83e54ad6e164911ca89a1331f405e` | identical | = HEAD |
| `data/derived/data/e428/alphaflow/repo_commit.txt` | `data/e428/alphaflow/repo_commit.txt` | `2a8e9d76b30fe98cd59707aa84021821dcc837de7d01a6e8fe382d5db34415ae` | identical | not tracked |
| `data/derived/data/e428/holdings/7DRH_removed.json` | `data/e428/holdings/7DRH_removed.json` | `631f4664bcc26fdef447977cee2bd72dbfcce71ec245cefa603a09a47469d652` | identical | = HEAD |

Plus 4,280 cached API response files under `data/derived/data/paper_verification/http/` (1,096), `data/derived/paper/arxiv/verification/cache/v2c/http/` (136), `data/derived/paper/arxiv/verification/cache/v5/http/` (815) and `data/derived/data/e428/http/` (2,233), each byte-identical to the file of the same name in the source repository (each file's name is the SHA-256 of its request; its `sha256` field is the SHA-256 of the response body).
