# Label cross-check with biotite (corrected census ledger)

An independent recomputation of every chain segment of the corrected
preregistered-census ledger (`data/derived/results/e427/batches`) with a
different mmCIF parser and lDDT implementation (biotite 1.7.1) and a
different mapping source: the residue-level SIFTS annotation of the PDBe
*updated* mmCIF files (`pdbx_sifts_xref_db_acc` / `pdbx_sifts_xref_db_num`)
instead of the segment-level PDBe mapping API the study used. No study or
audit code is imported. It backs the sentences in the paper on the second
label check and on the residues the labeller excludes.

## Result (run 2026-10-09)

| Quantity | Value | Source |
|---|---|---|
| corrected chain segments compared | 495 (123,454 rows); 22 further groups the study rejected were recomputed too | `full_ledger_out.json` |
| segments with identical positions and every lDDT within 1e-4 | **494 of 495** (123,390 rows) | same |
| largest \|Δ lDDT\| in those 494 | 5.0e-5 (the ledger stores 4 decimals) | same |
| remaining segment | 21xg chain A (Q4KH59) | `analyse_21xg_out.json` |
| 21xg A: lDDT values that differ | 18 of 63 shared positions, largest 0.0536 | same |
| 21xg A: chain median | 0.9405 (ledger) vs 0.9375 (biotite), difference 0.003 | same |
| 21xg A: residues below 0.60 | 5 vs 5 (one residue crosses 0.60 in each direction) | same |
| observed residues excluded as modified (other than MSE) | **57** (MLY 30, FME 12, TPO 6, LLP 3, TH5 2, A1E7L 1, A1FDA 1, CME 1, CSD 1; 14 entries) | `nonstd_out.json` |
| observed residues excluded as altloc other than A | **16** (10 entries) | same |

**21xg chain A.** Two crosslinked residues (`A1E7L` at label 17, UniProt 48,
and `A1FDA` at label 24, UniProt 56) are not standard amino acids. Next to
the second one, the study's sequence alignment has a tie: it pairs the
reference Thr at label 23 (UniProt 55 in residue-level SIFTS) with AFDB
position 56 (also Thr), so the ledger has position 56 where residue-level
SIFTS has 55. The identity gate cannot see this (Thr = Thr). Because lDDT
is computed over the aligned set, 18 neighbouring values move (by up to
0.054), the chain median moves by 0.003 and the number of residues below
0.60 is unchanged (position 53 drops below 0.60 in the ledger only; the
55/56 residue is below 0.60 in biotite only). `alignment_label_seq_1_40`
in the output lists the alignment around the tie.

The cross-check outputs contain no pLDDT value and no pLDDT-binned or
coverage quantity (the AFDB B-factor column is not read).

## Files

| File | What it is |
|---|---|
| `full_ledger_biotite.py` | per-segment recomputation of the whole corrected ledger → `full_ledger_out.json` (`full_ledger_out.log`: its printed summary) |
| `analyse_21xg.py` | the 21xg A discrepancy in detail → `analyse_21xg_out.json` |
| `nonstd_residues.py` | residues excluded as modified or altloc-not-A → `nonstd_out.json` |
| `fetch_updated_mmcif.py` | downloads the PDBe updated mmCIF files and compares them with `updated_mmcif_sha256.json` |
| `updated_mmcif_sha256.json` | SHA-256 of the 223 updated mmCIF files used for the recorded run |

## Running it

Requirements: Python 3.12 with `biotite==1.7.1` (it brings numpy 2.x; use a
separate virtual environment from the one in `requirements.txt`), and
network access for the two downloads.

```sh
python3.12 -m venv .venv-biotite && .venv-biotite/bin/pip install biotite==1.7.1
python -I code/tools/make_workroot.py _work
python -I data/fetch/fetch_raw.py --set e422 --work _work          # census source bytes (reference mmCIF, SIFTS JSON, AFDB models)
python -I audits/label_crosscheck_biotite/fetch_updated_mmcif.py --derived data/derived --out-dir _updated_mmcif
D=audits/label_crosscheck_biotite
.venv-biotite/bin/python -I $D/full_ledger_biotite.py --derived data/derived --raw _work/results/e422/raw --updated _updated_mmcif --out rerun_full_ledger.json
.venv-biotite/bin/python -I $D/analyse_21xg.py      --derived data/derived --raw _work/results/e422/raw --updated _updated_mmcif --out rerun_21xg.json
.venv-biotite/bin/python -I $D/nonstd_residues.py   --derived data/derived --updated _updated_mmcif --out rerun_nonstd.json
```

PDBe regenerates updated mmCIF files weekly; `fetch_updated_mmcif.py`
reports which downloaded files differ from the recorded ones. The census
source bytes are verified against their recorded SHA-256 on every read.
