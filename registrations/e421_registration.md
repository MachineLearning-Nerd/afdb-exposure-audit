# e421 registration — B-factor F-channel, post-AFDB candidate universe

> **Status: EXPLORATORY_EXECUTED — exploratory census completed with incomplete provenance.**
> Results: 24,457 residue labels across 96 entries (A=3411, B=1, OC=176, O=223).
> Verdict: E421_UNDERPOWERED_POOL. NOT admissible scientific evidence for the AD-CP claim
> (provenance envelopes incomplete). Registered positive claim requires a new prospective
> registration with full provenance.
> **E-number:** e421 (allocated by main, 2026-09-05)
> **Predecessor:** e420 (ATLAS-route G1, verdict E420_INFEASIBLE_CURRENT_POOL — temporal rule voided retrospective E-channel)
> **This registration changes:** candidate universe (post-AFDB PDB entries, not ATLAS census) + F-channel (B-factors from the same structure, not ATLAS RMSF). ATLAS is REMOVED from scope.
> **Everything else carries forward:** E-channel (AF vs reference, Kabsch + lDDT), temporal rule, retention, non-commercial, public endpoints, governance infrastructure.

## 1. Hypothesis

pLDDT can predict the flexibility–error relationship in a leakage-free protein structure benchmark:
where the protein is crystallographically flexible (high B-factor), AlphaFold's model is more likely to
deviate from the experimental structure; where rigid (low B-factor), AF is more likely to be correct.

## 2. Candidate universe

X-ray PDB entries whose initial release date strictly post-dates the modelCreatedDate of their AFDB model.
Assessment: 33,517 X-ray entries released after 2022-06-02 at resolution ≤ 2.5 Å (RCSB search API, 2026-09-05).

## 3. Data sources

| Source | Role | Authorization |
|---|---|---|
| AFDB | AF model (cifUrl child from /prediction/{accession}) | EBI terms; CC-BY-4.0 |
| PDBe entry-files | Reference structure + B-factors ({entry}.cif) | EBI terms; wwPDB CC0 |
| PDBe mappings | SIFTS UniProt crosswalk | EBI terms |
| ATLAS | **NOT USED** — removed from e421 scope | N/A |

## 4. Label channels

- **E:** AF model vs experimental structure (Kabsch + 4-threshold lDDT-Cα, frozen thresholds from e420)
- **F:** normalised B-factor from the SAME experimental structure (per-structure z-score or percentile; normalisation method frozen before POOL_FREEZE)
- **Combined:** OC (correct+ordered), A (correct+flexible), B (error+ordered), O (error+flexible), U (any endpoint unknown)

## 5. Temporal rule

modelCreatedDate STRICTLY earlier than the reference structure's initial release date
(in-file `_pdbx_audit_revision_history` ordinal 1). Pre-AFDB structures → typed
E420_AF_TEMPORAL_FAILURE → excluded (no repair).

## 6. Support floors (unchanged from e420)

20 proteins + 20 groups per primary class; per-protein ≥10 unique non-U residues;
per-fold ≥4 A + ≥4 B groups; ≥200 unique residues per cell; OC and O presence.

## 7. Fences

All e420 fences carry forward: no raw redistribution, public endpoints only,
no credentials, typed failures only, no silent imputation, no mutable latest,
predictor-blindness until POOL_FREEZE, retention 90-day event + 2027 backstop.

## 8. Lifecycle gates

(a) checkpoint freeze → (b) five-role reviews → (c) detached audit → (d) anchor
(OTS) → (e) bounded dispatch → (f) POOL_FREEZE → (g) label construction →
(h) G1 evaluation + independent audit.

### Amendment (2026-09-05) — minor clarifications from gate-(b)/(c) review

- Normalisation: frozen to **robust z-score (median/IQR)** per structure; the
  registration §4 "z-score or percentile" phrasing is superseded by the
  label_rules.json freeze.
- Analysis plan §5 "entry-size < 5000 atoms" pre-census filter is a pure
  **cost bound** (bounded download size), not a biological selection criterion.
  It does not affect the registration universe predicate.
- Gate (e) phase (ii) round 5 receipt = the binding receipt for grammar
  validation (10/10 PASS, 0 typed failures).
- cap_bound in the census artifact **must derive from the bound pool_freeze**
  (selection.cap_bound), never hardcoded (registration 15.3).

### Amendment (2026-09-08) — verdict-label correction (detached final audit)

- The header verdict string previously read `E420_UNDERPOWERED_POOL` — the
  wrong run's prefix (e420's cap-induced-shortage verdict; e421 records
  `cap_bound=false`, so that verdict class is mechanically inapplicable).
  Corrected to **E421_UNDERPOWERED_POOL**, matching the census artifact
  (`results/e421/four_cell_census.json`), the analysis plan §2, and the
  registered e421 vocabulary (`docs/e421_no_source_config.json`).
- Note: the **post-run checkpoint payload**
  (`docs/e421_post_run_checkpoint.json`, aggregate digest `c103eb06...`)
  carries the same wrong-prefix verdict plus `temporal_cause` /
  `outcome_equivalent` fields; the checkpoint is digest-bound anchor evidence
  and is deliberately NOT mutated — this amendment is the correction of
  record (see `docs/audit_e421_final_audit_20260908.md`, finding P-2).
