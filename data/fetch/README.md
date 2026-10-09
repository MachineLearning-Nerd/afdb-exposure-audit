# Re-fetching the raw inputs

Raw third-party files are not redistributed here: PDB/PDBe mmCIF files, PDBe
SIFTS mapping responses, AFDB model files and metadata, UniProt sequences and
the ATLAS archives. Every request the pipelines made was receipted with its
URL, method, HTTP status, UTC time and the SHA-256 of the response body. Those
receipts are part of the derived data, so the exact bytes used can be
identified and, where the source still serves them, re-fetched and checked.

## How digests verify identical bytes

For each receipt the original pipeline stored the response body under a file
named by its SHA-256 (`raw_name` = `sha256` in the receipt). `fetch_raw.py`
downloads the same URL again and hashes the body:

- the same SHA-256 means the bytes are identical to those analysed, and the
  file is stored under the recorded name, where the registered replay code
  (`code/experiments/e427_replay_relabel.py`, `e427_offline_relabel.py`,
  `code/verification/v5_templates.py`, `v7_trigger_prevalence.py`) expects it;
- a different SHA-256 means the source has changed (PDB entries are
  re-versioned, AFDB re-issues files). The body is kept under
  `<raw_dir>/_changed/<recorded>.<new>` and never under the recorded name, so a
  replay that needs it stops with a missing-file error instead of silently
  analysing different bytes.

Spot check (2026-10-09): the three files behind Figure 1 (9qj6 mmCIF, 9qj6
SIFTS mapping, AF-K7PQ54-F1-model_v6) re-fetched with identical SHA-256.

## Receipt sets

```sh
python -I code/tools/make_workroot.py _work            # scratch copy of the layout the code expects
python -I data/fetch/fetch_raw.py --set 9qj6 --work _work         # 3 files, Figure 1 / test_e427_label_repair
python -I data/fetch/fetch_raw.py --set e422 --work _work         # ~1,070 files: preregistered census, batches 1-4
python -I data/fetch/fetch_raw.py --set e421relabel --work _work  # ~410 files: post-2022 relabel (October 2026 bytes)
python -I data/fetch/fetch_raw.py --set e420census --work _work   # ~3,600 files: ATLAS availability census (AFDB metadata, PDBe)
python -I data/fetch/fetch_raw.py --set e420labels --work _work   # ~1,240 files: ATLAS reference / AFDB coordinate files
```

Add `--dry-run` to list what would be fetched, `--limit N` to sample.

| Set | Receipts (in `data/derived/`) | Stored under (in the work root) | Used by |
|---|---|---|---|
| `e422`, `9qj6` | `results/e422/batches/e422_batch0{1..4}_checkpoint.json` → `receipts` | `results/e422/raw/<sha256>` | e427 replay/offline relabel, V5, V7, 9qj6 audit, `test_9qj6_real_bytes_repaired` |
| `e421relabel` | `results/e427/e421_relabel/receipts.json` | `data/e427/e421_raw/<sha256>` | `e427_e421_relabel.py` |
| `e420census` | `results/e420/source_availability_manifest.json` → `requests` | `data/e420/census_raw/<sha256>` | `e420_overflow_check.py`, `e420_label_run.py` |
| `e420labels` | `results/e420/label_run_receipt.json` → `requests` | `data/e420/labels_raw/<sha256>` | `e420_label_run.py` |

Not covered by the script:

- **RCSB Search API POST queries** (census selection, homology search). The
  census queries are rebuilt by `e423_runner.build_census_query`; a later
  query returns more entries because the archive grows. The homology-search
  responses used for Table 2 are included verbatim (RCSB data are CC0) in
  `data/derived/data/paper_verification/http/` and
  `data/derived/paper/arxiv/verification/cache/v2c/http/`, keyed by
  SHA-256(method, URL, body), so V2/V2c replay them without network; only the
  UniProt sequences are fetched live.
- **UniProt FASTA** (`https://rest.uniprot.org/uniprotkb/<acc>.fasta`): fetched
  by `v2_prior_exposure.py`; receipts with SHA-256 are in
  `out/v2_prior_exposure.json` and `out/v2c_extension.json`.
- **V5 template lookups** (PDBe SIFTS `api/mappings/uniprot/<pdb>` and RCSB
  `rest/v1/core/entry/<pdb>`): re-fetched by `v5_templates.py`; receipts in
  `out/v5_templates.json`.
- **ATLAS census universe**: the ATLAS parsable archive
  (`https://www.dsimb.inserm.fr/ATLAS/api/parsable`, retained snapshot
  `atlas_parsable_full_20260830.zip`, Last-Modified 2024-11-10,
  SHA-256 `c9cba2b7190676814bb83daf11fcfe92ad0c3a3cdaca4eb4782ca6b98c67b75e`;
  member `ATLAS_parsable_latest/2023_03_09_ATLAS_pdb.txt`,
  SHA-256 `667a5ebc28cd3f160a4fa313fe713d337dc1148aca877216dd9153f44a68cd0b`).
  The e420 code also validates against retained endpoint-contract snapshots
  (`data/source_contracts/*`, listed with digests in
  `registrations/e420_registration.md`); these are not redistributed.

## Provenance caveat

The preregistered census (e422) was replayed from its original stored bytes,
which are identified by the receipts above. The post-2022 relabel (e427 §4.5)
used files fetched again on 2026-10-07, because the September 2026 bytes of the
exploratory e421 census had not been kept; its receipts identify the October
bytes, not the September ones.
