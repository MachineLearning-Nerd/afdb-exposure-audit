# Main's bounded Rev-3 table recount: V1, V2, V5

Date: 2026-10-08 (Asia/Kolkata)
Scope: arithmetic and wording audit of the paper's ATLAS exposure, pre-cutoff
sequence exposure, and AFDB template-exposure summaries. This is a narrow
derived-artifact audit, not a source-data audit, full-paper review, or
publication approval.

## Reproduction

The registered V1 parser was copied with its two input JSON artifacts to a
fresh audit directory and run using the repository's supported Python
3.12.13 virtual environment in isolated mode (`-I`). Its output SHA-256 is
`f0e5e255b8b0fc19024fe93eba0820036a75942e886e8484faf1938a6acf5fa0`.
The recount reported an empty discrepancy map against the registered V1
expectations. Input hashes:

| Artifact | SHA-256 |
|---|---|
| `results/e420/label_run_receipt.json` | `b2c36992ec965b77a0d040e44c25fcebbd1ff1fea03a7e3e9a1747d6f593a452` |
| `results/e420/overflow_temporal_check.json` | `6c8cf0a26757a252d44c500ca1619812d1c6349a0bf1b24313ad685bf06f95d7` |
| Registered parser `paper/arxiv/verification/audit/v1_channel_recount.py` | `57ac8a57d5cd01aa4427e6919312b564e49a46b274d707e2830d68b7040caf88` |

V2 and V5 were independently reaggregated in a fresh isolated Python 3.12.13
process from copies of the two existing derived JSON outputs, using
`main_recount_v2v5_20261008.py` (SHA-256
`cc5aac99181e69931b6a263602e842d74e0e190c8cc61d1f4e7291f5eca622db`). The
recount output SHA-256 is
`acf2152d5879ddb2729fb341d9511d4af02304ef59654fb642505d96a40e3bb8`.
All seven checks passed, including row-level V5 reaggregation against its
registered summary and the documented template-date cutoff. The inputs were:

| Artifact | SHA-256 |
|---|---|
| `paper/arxiv/verification/out/v2_prior_exposure.json` | `19fcabc7a36bc73f02fa5be0f963e53cc84cdc702c80285c7ad7fe88878ce720` |
| `paper/arxiv/verification/out/v5_templates.json` | `d8b184c902819f2f7969ea19778b5bf3496438eb9de915b88f294f1c6ef7e69a` |

The reaggregation script is a working-copy file at
`paper/arxiv/verification/audit/main_recount_v2v5_20261008.py`. The fresh
outputs `v2v5_recount.json` and `v1_fresh_py312.json` are under the isolated
sibling directory `/home/dineshai/Drives/Code/moluq-paper-audit.1OdlGJ/`, not
inside the repository. No network request, raw source structure, label row,
pLDDT value, or P5 output was read by these recounts.

## Reproduced values

- V1 admitted pool: 553/599 accessions (92.3%) have all dated references on
  or before the AF2 training cutoff (2018-04-30); the registered
  model-created-date rule classifies 599/599 before that model date. At the
  accession-entry level, 595/646 references satisfy the training cutoff and
  646/646 satisfy the model-date rule.
- V1 overflow: 238/269 accessions (88.5%) satisfy the training cutoff;
  264/269 satisfy the model-date rule. At the entry level, the counts are
  256/289 and 284/289, respectively.
- V1 union: 791/868 accessions (91.1%) satisfy the training cutoff and
  863/868 satisfy the model-date rule. At the entry level, the counts are
  851/935 and 930/935. The accession-level channels disagree for 72 cases;
  the median dated reference release year is 2008.
- V2's 138 E422 accessions classify as SAME_ACCESSION=78, SEQ95=9,
  SEQ30=41, NOVEL=10. Thus 87/138 (63.0%) are SAME_ACCESSION or SEQ95, and
  128/138 have a reported hit at 30% identity or above.
- V5's registered entry-level template classes reproduce as SAME_ACC=45,
  OTHER=88, UNRESOLVED=5 across 138 accessions. All 138 records list four
  templates. The latest recorded template release is 2021-01-20, earlier
  than the registered/documented 2021-02-15 limit.

## Finding: qualify the V5 “same-protein template” prose

The registered V5 rule assigns SAME_ACC when any listed template PDB entry
has a SIFTS UniProt mapping to the target accession. It does not require
that the particular template chain named in AFDB's model details be the chain
carrying that accession mapping. The paper currently describes 45 models as
having been “given a template of the same protein” and says that a model has
such a template if “any template maps to its accession.” That wording can be
read as chain-specific evidence, which the registered entry-level test does
not establish.

The V5 artifact contains a chain-specific count explicitly marked
`NOT_PRESPECIFIED`: SAME_ACC=44, OTHER=89, UNRESOLVED=5. Reaggregation of
the per-template `template_chain_maps_to_accession` flags reproduces 44; the
one registered SAME_ACC accession without a matching named template-chain
flag is Q5UE59. This is a sensitivity analysis, not the registered estimate,
so it must not replace the registered 45/138 result. For release, either
qualify the claim as “45 models had at least one listed template PDB entry
with a SIFTS mapping to the target accession under the registered entry-level
rule,” or independently establish the corresponding chains before retaining
the stronger chain-specific wording. Any report of 44 must be labelled
post-hoc/not prespecified.

No LaTeX manuscript source was edited in this audit because the author
requested coordination before changes to `temporal_leakage.tex`. This finding is
nonblocking for the arithmetic recount but is an unresolved precision issue
for the manuscript's interpretation of V5.

## Limits and disposition

The recount verifies arithmetic on the listed derived artifacts; it does not
validate the original external responses, their licensing or scientific
meaning, the complete E421/E422 label pipeline, P5, or the manuscript's other
claims. The reaggregation code and this report are working-copy artifacts and
are not committed or independently audited. The corrected E421 census remains
exploratory; the corrected-label E422 P5 result remains unverified here; E422
remains paused. G1-G4 remain unmet. No claim was promoted, and HOLD /
NOT_GRANTED / NON-FINAL NO-GO is unchanged.
