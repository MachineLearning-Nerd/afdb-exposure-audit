# Independent recount of E422 corrected-label F-scope counts

Date: 2026-10-08 (Asia/Kolkata)

## Scope and method

This is an arithmetic audit of the registered F-scope before/after statistic,
not validation of the label repair, external source responses, licensing, or
the corrected-label P5 evaluation. It uses the original E422 and corrected
E427 batch 1–4 checkpoint JSON files already in the repository. No raw source
structures, network access, pLDDT values as predictors, or P5 evaluation output
were used.

The independent stdlib script
[`recount_e422_fscope.py`](recount_e422_fscope.py) was run in a fresh Python
process with `-I` and working directory `/tmp`. It checked each file's raw
SHA-256, schema, batch index, `n_rows`, and canonical checkpoint digest
(SHA-256 of canonical JSON after removing the top-level `sha256`). It then
aggregated finite residue-level lDDT by the registered key
`(entry, chain, accession)`, took each key's median, and compared corrected
minus original medians. “Changed” means exact nonzero difference; the registered
large-change count uses the strict threshold `abs(delta) > 0.2`.

Command:

```text
/home/dineshai/Drives/Code/AllCode/moluq/.venv/bin/python -I \
  /home/dineshai/Drives/Code/AllCode/moluq/paper/arxiv/verification/audit/recount_e422_fscope.py \
  --original-dir /home/dineshai/Drives/Code/AllCode/moluq/results/e422/batches \
  --corrected-dir /home/dineshai/Drives/Code/AllCode/moluq/results/e427/batches
```

## Recomputed result

All eight checkpoints passed their row-count and canonical self-digest checks.
The result reproduces the table/figure summary totals: original 117,905 rows,
488 groups and 138 accessions; corrected 123,454 rows, 495 groups and 141
accessions. There are 488 groups in common, zero original-only groups, and
seven corrected-only groups. Among the 488 common groups:

- 188 have any nonzero median-lDDT change;
- 151 have `abs(delta) > 0.2`, across 40 accessions.

The `151` and `40` results match both the existing registered F-scope
summaries. The `188` exact-change count matches
`figures_rev3_numbers.json` (`segments_changed_any`).

## Finding: manuscript's “217 changed” is not reproduced

The current manuscript says “217 of 488 chain segments” changed in median
lDDT at [the abstract/summary text](../../temporal_leakage.tex#L42) and repeats
217 in the E422 correction discussion (currently around lines 550–552).
Under the registration's `(entry, chain, accession)` grouping and exact
nonzero-change definition, the independently recomputed count is 188, not
217. The existing figure-number artifact also records 188. The separate
large-change statement—151 groups across 40 accessions with `abs(delta) >
0.2`—does reproduce.

I did not edit the manuscript: the source of 217 has not been identified, and
it could reflect a different grouping or earlier calculation. Before paper
release, reconcile that number with its intended definition; if no distinct
registered definition and reproducible source are found, replace it with the
registered 188 count and retain the 151/>0.2 result. Do not present 217 as
verified on this audit.

## Exact inputs and cross-check artifacts

| Input | Raw SHA-256 |
|---|---|
| Original batch 1 | `526fd4c54603ae128a2b6313a09c10c7a18e5f5699af3a3d7fc46bfb81c5d029` |
| Original batch 2 | `08cc71203d3830813f5b096e27a73e6d1e53cfc1776b20a0a04cb95c772a4aea` |
| Original batch 3 | `b31e0ad5a00ea6e12c8c809d51f59e5eba75966f3dd8f38f42b61f51b6d61c3c` |
| Original batch 4 | `a08e538db7dd54f2a046f2bbc7ed3c9f37e6c2e1ba2c30f42403aa66f6fa934b` |
| Corrected batch 1 | `526fd4c54603ae128a2b6313a09c10c7a18e5f5699af3a3d7fc46bfb81c5d029` |
| Corrected batch 2 | `98b5ccb8af4e7899f665d703d1f55f75cda3d2387391a92199e9817881548935` |
| Corrected batch 3 | `7153a7f69b7cbd3d4c894b9a9638ef3e28d9c7dea3979960fc0b4bd189f8c8d3` |
| Corrected batch 4 | `2a2810a02d9b5a14ac70bebe5a92462353da81b229b875971ff68ec7720330d8` |
| Existing V7 linked F-scope output | `37efe7281f13173493432aaf10cba595aa8c7b255c5d519ae856d9b265bf1e16` |
| Existing figure numbers | `43ee0bc7f3056a2448e321d3fbc8f3c02af38b1793274ba4779ee574cb0e8246` |
| Manuscript bytes reviewed | `60de181923e29e75fd0523a08112c39d7a01f5823d16b01518ce9f9a3e42f2ae` |

The original and corrected E422/E427 evidence remain distinct; this recount
does not transfer E421 source or audit credit. E421 remains exploratory, the
corrected-label P5 rerun still requires its registered output and audit, and
all publication gates remain HOLD / NOT_GRANTED / NON-FINAL NO-GO.
