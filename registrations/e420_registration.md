# e420 registration — G1 successor ground truth for the AD-CP line (fresh predictor-blind cohort)

> **Status: DRAFT — registration only; no source access yet or permitted.**
> e-number **e420** was allocated by `main` on 2026-09-04 (sole active session;
> single-session mode per `AGENTS.md`) to unblock gate **G1**. This document
> instantiates the successor-registration requirements of
> `docs/G1_UNBLOCK_PACKET.md` §7 under a main-allocated e-number. It does not
> authorize any fetch, source read, cache write, candidate retrieval, predictor
> invocation, embedding, label, model, or conformal operation. Every lifecycle
> gate in §10 must complete in order before any current provider byte is
> requested. G1–G4 remain **NOT MET**; the package remains
> **HOLD / NOT_GRANTED / NON-FINAL NO-GO**.

- **e-number:** e420 (allocated by `main`, 2026-09-04)
- **Date:** 2026-09-04
- **Owner:** `main` (single-session mode: main sequentially performs the
  manager, builder, reviewer, audit, and scribe functions while preserving the
  registered lifecycle, fresh-context review, and independent-audit separation)
- **Task:** G1 unblock — `docs/G1_UNBLOCK_PACKET.md` §10 blocker 3
  ("Registration: no current main-allocated successor registration has passed
  exact-byte review and detached audit").
- **Predecessor lineage:** e413 terminal infeasible -> e418/e419 contract
  captures -> this successor. The e413 pool is terminally infeasible
  (`E413_INFEASIBLE_CURRENT_POOL`) and is **not reusable, rescuable, or
  extendable**; no e418/e419 tuple, fixture receipt, review, checkpoint, or
  audit transfers as authority.
- **Authoring basis:** local repository records only. No network, source, API,
  cache, candidate, sequence, coordinate, label, embedding, or model access of
  any kind occurred while drafting this registration.

## 0. Binding context and lineage

### 0.1 Authoritative requirement sources

| Source | Role for e420 | SHA-256 |
|---|---|---|
| `docs/G1_UNBLOCK_PACKET.md` | Authoritative statement of what G1 requires; §7 successor draft is instantiated here | `019162989c339fdd3b0d1c41c60188d0fa9c4b52a7fa2885e2bba1e2aea51e5c` (committed 2026-09-04, commit 7f04c3e6; living document amended after draft — supersedes draft-time `1e4e2868…`; §4–7 instantiation re-verified against committed bytes by the theory/audit review) |
| `docs/adcp_roadmap.md` | G1–G4 definitions; Phase-1 ground-truth requirements | `2ef939aabec0156918c65c4fed13c7d4c98bc91220fc8d347d278f5a1c7beebb` |
| `docs/human_gateC_authorization_20260904.md` | Binding human policy layer: scope, pins, permissions, retention, rights, anchor | `64bdaa4c6496db31c43840374841f4decd8fe1aba28abbaa4c0aa2bcf7ae96fa` |
| `docs/e413_registration.md` | Prior ground-truth design; support floors and guards carried forward (pool NOT reusable) | `95c2c4d4dcca7823808629f140414c15ed093f0da78886a4b438b07a97fecc86` |
| `docs/e413_results.md` | Terminal feasibility result of the predecessor | `646eed99bf6d414af2a3bc0241a14a7e62be17e02ed159fcf7da81daa152b351` |
| `docs/e418_registration.md` | Contract/wire-design context (superseded legacy field names; historical design input only) | `931640ced65f2f664658c11e664879372e64ce21757086da942f1f7c157db71a` |
| `docs/e418_registration_amendment_2026-08-28.md` | A1.3 explicit wire roles (RCSB search pages, AFDB child, HEAD availability, UniProt FASTA grammar) | `cb3872db3446d4c447679f424e7b95ef3f9193aa19855c26ef076bee8e997640` |
| `docs/e418_registration_amendment_2026-08-29_analysis.md` | Phase-1 label / held-out split / AD-CP-vs-pooled analysis contract (later-gate context) | `b918e4834856e54f1e1e72ba82056c24d63dd6c51b6a1dc53e321e009b2bb008` |
| `docs/e418_registration_amendment_2026-08-30.md` | v2 contract rebind candidate; ATLAS gap semantics; production authority | `8bd574aaf3eb0f1e52cb8830a27ec91c98f77d638795b06debf8afecc23df757` |
| `docs/e419_registration.md` | Retrospective diagnostic design context (no authority transfer) | `b717c36a03ffc56662545ce486329553e0ae80674a9965e433024939433da83b` |
| `docs/e419_registration_amendment_2026-08-29.md` | Strict all-pair diagnostic scope; Kabsch/lDDT restatement | `52db12e1287d2c4884ddcb224ef24a06b305b457dd428bd55be83e3cf574fb51` |
| `docs/e418_registration_amendment_2026-09-04_gateC.md` | Gate-C policy-binding pattern (PROPOSED/NOT_EFFECTIVE; pattern only) | `b9c720e79219cf1f3d28ebfd590f124bcc73285a0fbcc94dc69c7603d46cfb4c` (committed 7f04c3e6; supersedes draft-time `6f59fafb…`) |
| `docs/e419_registration_amendment_2026-09-04_gateC.md` | Gate-C policy-binding pattern (PROPOSED/NOT_EFFECTIVE; pattern only) | `5148478a147a54e8a46690318c681f2457227904dcee9553ae0cedbb869bc283` (committed 7f04c3e6; supersedes draft-time `d77bdb43…`) |
| `docs/data_source_contracts_20260830.md` | Captured contract notes; receipt lineage; fences observed | `e05790fdb95b7dcba416d540ad6867f9e972bb3e048f127b3256557010d4cfa2` |
| `docs/theory_e418_amended_contract_fields_spec.md` | Amended provider-contract field/snapshot spec; provenance-schema lineage | `c074348157b6b83a93285d716e1bd62db1c7c4ffe3f72bea79f0423fff442996` (working-tree value) |

SHA-256 values above were recomputed from the working tree at draft time and
are now COMMITTED (requirement-layer commit `7f04c3e6`, 2026-09-04 — review
repair M1/M2/TM1/TM2): the checkpoint binds exact committed bytes for every
file listed here that is in its `files[]`. Draft-time packet/gateC pin values
(`1e4e2868…`, `6f59fafb…`, `d77bdb43…`) are superseded by the committed values
above (living documents amended after draft); the contract-spec bytes
(`c0743481…`) are pinned here and adopted by e420 as the binding semantics
reference for "contract spec §1/§3/§4" citations — authority for the e418 lane
is unchanged (O10).

### 0.2 Predecessor lineage statement

1. **e413 (terminal):** the frozen census contained 172 target rows, 3
   accession-bearing proteins in 3 groups, class-A/B support 0/0 proteins and
   0/0 groups, all four 2x2 availability cells zero; exact verdict
   `E413_INFEASIBLE_CURRENT_POOL` (`docs/e413_results.md`). Its pool, CASP
   releases, and releases are not reusable authorization.
2. **e418/e419 (contract captures):** design-only registrations plus
   amendments captured and pinned the current PDBe v2 / AFDB / ATLAS contracts
   (receipt lineage now v54; provenance schema v27) and the amended AFDB
   temporal/provenance disposition. Both lanes remain
   PROPOSED/NOT_EFFECTIVE/HOLD; their gate-C amendments
   (`docs/e418_registration_amendment_2026-09-04_gateC.md`,
   `docs/e419_registration_amendment_2026-09-04_gateC.md`) bind policy values
   only and are non-authorizing.
3. **e420 (this registration):** fresh predictor-blind pool, fresh release
   pins under the human gate-C record, fresh groups/folds, fresh checkpoint
   and external anchor. Nothing is inherited except explicitly cited contract
   context and the carried support floors/guards of §6.

### 0.3 e-number label note (recorded, unresolved disposition)

Several 2026-08-29/30 `docs/DECISIONS.md` entries informally used the label
"e420/T-038" for a **future Phase-3 quotient-space conformal method**
registration ("Register/execute E420 only after G1"). Main's 2026-09-04
allocation assigns **e420 to this G1 ground-truth successor** instead. No
Phase-3 authority, protocol, or artifact transfers to e420. An append-only
DECISIONS entry recording this label reassignment (and that the Phase-3 method
will receive a fresh allocation only after G1) is an OPEN item (§13, O8).

## 1. Purpose and non-goals

**Purpose.** Establish the AD-CP Phase-1 real ground truth required by gate G1:
a fresh, predictor-blind, source-backed protein cohort with independent
coordinate-error (E) and dynamics/flexibility (F) label channels, complete
raw-byte provenance, group-atomic held-out folds, and adequate per-class
support — under the human gate-C policy layer
(`docs/human_gateC_authorization_20260904.md`).

**Non-goals.** e420 does not: construct or promote any label-free
flexible-vs-error predictor; run ESM inference, embeddings, absorption scores,
model fitting, or conformal calibration; compare AD-CP with a CalPro-style
marginal interval; register or execute Phase 2/Phase 3; or promote any paper
claim. Those require separate registrations after a G1 pass (packet §7.8).
A source-availability or support count is a feasibility eligibility count,
never a biological state or a label.

**Label-channel precedence note.** The program-level Phase-1 sketch in
`AGENTS.md` §4 names "SIFTS↔UniProt + PDB B-factors (missing-residue =
disordered, high-B = flexible, low-B = error), cross-checked against ATLAS
MD-RMSF". The authoritative G1 requirement supersedes that sketch for this
gate (`docs/G1_UNBLOCK_PACKET.md` §7.3–7.4; `docs/adcp_roadmap.md` Phase 1):
**E** is AF-versus-released-reference coordinate disagreement; B-factors,
pLDDT, absorption, and any AF-derived feature are barred from E; the
same-reference missing/construct-span signal is a **sensitivity role only**
and cannot satisfy primary flexible support; **F** is the independent ATLAS
three-replicate MD-RMSF channel; SIFTS/UniProt mappings provide
identity/numbering joins only. This registration follows the packet.

## 2. Hypotheses and pre-registered G1 gate

### 2.1 H1 — hypothesis (ground truth established on a predictor-blind pool)

After the frozen predictor-blind availability census, the declared source-
opening chronology (§10 gates f–g), and conservative label construction, the
cohort contains:

- at least **20 unique accession-bearing proteins and 20 independent homology
  groups** in **A** (`E=correct, F=flexible`);
- at least **20 unique accession-bearing proteins and 20 independent homology
  groups** in **B** (`E=error, F=ordered`);
- at least **one protein and one group in each of OC, A, B, and O**;
- at least **200 unique non-`U` residue identities in each cell**;
- at least **10 valid residue positions per protein** supporting A or B;
- **A and B support in all five held-out group folds**, with at least four A
  groups and four B groups in every fold; and
- complete raw-byte provenance for every source row and an independent audit
  of the support census.

The four cells are retained and never pooled:

| Cell | Coordinate-error state E | Dynamics state F | Interpretation |
|---|---|---|---|
| OC | low/error-correct | low/ordered | correct and ordered |
| A | low/error-correct | high/flexible | primary flexible-correct |
| B | high/error | low/ordered | primary AF-wrong/ordered-error |
| O | high/error | high/flexible | retained overlap |

`U` is used for missing, ambiguous, conflicting, unmapped, unavailable, or
otherwise ineligible evidence. Empty cells remain visible and fail the
registered support gate; they are not hidden by relabeling or pooling.

An H1 pass is a cohort feasibility result. It enables main's separate G1
admission decision after independent audit (§10 gate h); it is not a method
claim and not a publication claim.

### 2.2 H0 — no-go

Any release, license, provenance, template-overlap, chronology,
identity/mapping, coordinate, dynamics/archive, group, fold, raw-byte,
support, mutation-test, or terminal-manifest failure is a reported no-go.
Missing values remain `U`. No row is imputed, silently dropped, favorably
rounded, relabeled as error, or rescued by changing the candidate set, the
release, the thresholds, or the folds.

### 2.3 Source separation

- **E** uses only the exact frozen AFDB model coordinate artifact for the
  exact accession/model entity and only a separately released experimental
  coordinate source approved per §3.5/O2. A predicted or simulated structure
  is never truth.
- **F** uses only the fresh authorized ATLAS per-chain analysis archive and
  its three registered RMSF replicates, joined through the declared
  ATLAS PDB/chain -> selected-chain SIFTS -> canonical UniProt correspondence.
- Neither channel may use pLDDT, AF confidence, AF templates as evidence, AF
  provenance as evidence, B-factors, sequence-only features, embeddings, or a
  downstream model score. The two channels never read one another's endpoint
  values; shared identity selection is disclosed, not relabeled as statistical
  independence.

## 3. Data: authorized sources, pins, and release policy

### 3.1 Source scope (human-authorized, non-commercial research line)

Exactly three sources are authorized
(`docs/human_gateC_authorization_20260904.md` §1): **PDBe (SIFTS mappings
served via PDBe API v2)**, **AFDB**, and **ATLAS**. Public, read-only GETs
only; no credentials; no paid APIs. ATLAS is CC-BY-NC-4.0, admitted for the
non-commercial research line and academic publication only; ATLAS-derived data
is excluded from any commercial Samma Labs artifact (fence restated in §11).
Access executes only through the registered workflow: this registration ->
five-role pre-access reviews -> detached audit -> external anchor
(OpenTimestamps, §10(d)) -> bounded gate-C capture dispatch (§10).

### 3.2 Pinned contract evidence (verified on disk 2026-09-04)

The binding pins are the receipt-v30-retained artifacts of the human policy
record §2. All eight were verified present under `data/source_contracts/`
with exactly these SHA-256 values at draft time:

| Source | Pin (file) | SHA-256 (verified) |
|---|---|---|
| PDBe v2 OpenAPI (info.version 2.10.5) | `pdbe_v2_openapi_20260830.json` | `e86ab8d3c2c9b6865f158526c93b742db344fa62ec7277dea1e4d17e6852899d` |
| AFDB API OpenAPI | `afdb_api_openapi_20260830.json` | `714607265fd8edc581baf28df038ea804d96d871baa6034ff60d22d0cf893163` |
| AFDB v6 release evidence (2025-09-15; UniProt 2025_03 sync) | `afdb_changelog_20260830.txt` | `b6cc3567ff0cd9232c94bfdc080c2871cc27142a39ec983d97ecc26b21b59b86` |
| AFDB canonical optional-types table | `afdb_optional_types_generated_20260830.json` | `060de272d3a0de2d0d9346f67cb03efab49f00db15363e90f26811a791176eb3` |
| ATLAS API OpenAPI (v2.0.0) | `atlas_api_openapi_20260830.json` | `824150cc10bd365c82176bc11884a168d02a64029da65b475620a44054d3ee44` |
| ATLAS parsable census snapshot (Last-Modified 2024-11-10) | `atlas_parsable_full_20260830.zip` | `c9cba2b7190676814bb83daf11fcfe92ad0c3a3cdaca4eb4782ca6b98c67b75e` |
| — census member `2023_03_09_ATLAS_pdb.txt` | CRC-validated member | `667a5ebc28cd3f160a4fa313fe713d337dc1148aca877216dd9153f44a68cd0b` |
| ATLAS MD parameters (Last-Modified 2024-11-10) | `atlas_md_parameters_20260830.bin` | `0fb6e4ffbbdc0f624d071b63bf89839d83e543acb1d21e75973b76cec850a951` |

**Receipt/schema generation note (OPEN O4).** The human record names contract
receipt v30 (`4b9c4fe3e3a2db0c0030c3f96325b96c9a89e68469d8cf5ff9b6ff192f54445d`).
The on-disk canonical pair has since advanced to **receipt v54**
(`data/source_contracts/capture_receipt_20260830.json`, schema
`e418-t037-source-contract-capture-receipt-v54`, 87 items, working-tree file
SHA `bfcca135256038265450cb8c045d1959b8abc29b349bafea0316fcdb501475f3`,
self-excluded final hash announced outside the file) paired with **provenance
schema v27** (`data/e418_retained_provenance_schema_20260830.json`, schema
`e418-retained-provenance-schema-v27`, SHA
`a223f271162fa921cffc23d4b040a603403b76df288cbf439021984df6b14efe`, 24,213 B,
marker `e418-provenance-schema-generator-v27`). The eight artifact pins above
are byte-identical across generations. The exact receipt/schema generation
pair bound by the e420 checkpoint must be fixed by main before checkpoint
freeze (§13, O4).

### 3.3 SIFTS release: typed unknown

The retained PDBe v2 OpenAPI exposes **no SIFTS release identifier**; the
assumed `X-PDBe-SIFTS-Release` header is `NOT_DOCUMENTED_IN_V2_OPENAPI` and is
**forbidden** as a guard (contract spec §1). The SIFTS release identifier is
recorded as a **typed unknown, never invented**. The honest pin is the PDBe
OpenAPI snapshot SHA (§3.2) plus per-response capture timestamps, retained raw
bytes, and byte hashes for every mapping response. Retrieval time alone is not
a release pin.

### 3.4 Candidate universe: retained ATLAS parsable census

The candidate universe is enumerated **only** from the retained census member
`2023_03_09_ATLAS_pdb.txt` (SHA `667a5ebc…`) inside the retained snapshot
`atlas_parsable_full_20260830.zip` (SHA `c9cba2b7…`). The registered ATLAS
release for e420 is the recovered `ATLAS_parsable_2023-03-09` member set
within that aggregate snapshot; e420 records explicitly that the source is
the recovered 2023-03-09 census inside a Last-Modified 2024-11-10 archive,
not a floating latest and not the download-page 2024-11-18 release-history
files. Membership is an ATLAS availability fact, not a predictor output, so
the universe is predictor-blind by construction. Per-chain profile
availability under the current API is captured separately at gate (e) with
typed success/failure rows; the census alone establishes nothing about
current endpoint behavior.

### 3.5 Registered source-role map (templates; exact bytes bound at checkpoint)

| source_role (literal) | method | endpoint / member rule | permitted use | license/terms role |
|---|---|---|---|---|
| `E420_ATLAS_CENSUS_UNIVERSE_V1` | retained bytes (no fetch) | member `2023_03_09_ATLAS_pdb.txt` of `atlas_parsable_full_20260830.zip` | candidate universe only | ATLAS CC-BY-NC-4.0 |
| `E420_PDBE_MAPPING_PRIMARY_V1` | GET | `https://www.ebi.ac.uk/pdbe/api/v2/pdb/entry/uniprot_mapping/{pdb_id}/{entity_id}` — raw envelope retained and hashed before normalization; per-residue strict join with typed rejections per contract spec §1 | selected-chain PDB↔UniProt crosswalk only | EBI terms (PDBe public data CC0 + EBI ToU; upstream third-party rights retained) |
| `E420_PDBE_MAPPING_SEGMENT_V1` | GET | `https://www.ebi.ac.uk/pdbe/api/v2/mappings/uniprot/{pdb_id}` — entity-id resolution source + guarded interpolation-only fallback; resolution use: selected entity's mappings must meet identity >= 0.90 / coverage >= 0.80; fallback use: identity == 100.0 exact; never flat-parser input | entity-id resolution + mapping fallback only | EBI terms (as above) |
| `E420_AFDB_METADATA_V1` | GET | `https://alphafold.ebi.ac.uk/api/prediction/{accession}` — one-object array per accession; per-row pins `modelEntityId`, `allVersions`/`latestVersion`, `modelCreatedDate`, `sequenceVersionDate` recorded at capture | model identity, temporal anchor, sequence/offset fields, returned coordinate URL | AFDB CC-BY-4.0 + EMBL-EBI terms |
| `E420_AFDB_COORDINATE_CHILD_V1` | GET | exactly the metadata row's **`cifUrl`** field (ModelCIF same-bytes path; field bound 2026-09-04 per review m12 — `pdbUrl`/`bcifUrl` are forbidden substitutes); parent hash/URL/model-tuple bound | E predicted-coordinate input only | AFDB CC-BY-4.0 |
| `E420_ATLAS_ANALYSIS_V1` | GET | `https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/{pdb_chain}` — whole authorized ZIP; member names/headers verified at capture against the frozen grammar | F/RMSF construction only | ATLAS CC-BY-NC-4.0 |
| `E420_ATLAS_METADATA_V1` | GET | `https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/metadata/{pdb_chain}` | chain/context metadata only (never a label) | ATLAS CC-BY-NC-4.0 |
| `E420_REFERENCE_COORDINATE_V1` | GET | `https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif` family (compressed `.cif.gz` + legacy variants; `.cif.gz` existence typed-checked at first capture) — **RESOLVED via Path A** (addendum `169032e8…`; capture `ddfa3e23…`); fetched at gate (g) for frozen-pool accessions only | E reference coordinates + in-file initial-release metadata (§5.1 derivation) | EBI terms + wwPDB upstream public-domain data-owner terms (captured) |

Wire-contract requirements carried from the retained contracts: PDBe and
AFDB adapters hash the exact raw envelope before normalization and reject
truncation, member/header change, duplicate/ambiguous objects, and hash
mismatch; a normalized JSON object never substitutes for the raw response.
ATLAS endpoints are **GET-only** (the retained contract capture verified
`HEAD /parsable` -> 405 `Allow: GET`). The ATLAS sequence/structure search
POST routes are **forbidden** (unregistered source/search step). No
open-ended search, latest-release fallback, URL synthesis, alternate member,
transport substitution, or source-role substitution is permitted.

**Retention policy (restated from the human record §4):** raw cached source
bytes are deleted once, with an exact path+SHA-256 deletion receipt, 90 days
after the earlier of (a) the final publication decision on the primary AD-CP
paper at its submission venue or (b) written project closure by the founder;
hard backstop **2027-12-31**. Kept indefinitely: all SHA-256 checksums, fetch
scripts, derived data, and capture receipts. Cache location is under `data/`;
large raw files are never git-committed.

## 4. Predictor-blind pool construction

### 4.1 Availability census and predictor blindness

The first post-dispatch run is a source-availability census only (packet
§7.2). It may inspect, for frozen availability queries: source existence,
release/version metadata where exposed, accession/chain identity, and typed
failure responses. It must not inspect pLDDT, AF coordinates, AF templates as
values, AF provenance as values, RMSF numeric values, B-factors, coordinate
discrepancies, labels, embeddings, or model outputs. Every attempted row
retains: opaque candidate/group key and source role; request/archive/member
identifier, release/version if present, response status; raw-byte existence
and hash status without reading scientific values; mapping/identity
eligibility status; exclusion reason or inclusion status; and
parser/config/registration/checkpoint hashes. Failed rows are retained; there
is no favorable re-query, release switch, pool expansion, one-site rescue, or
post-hoc group merge.

### 4.2 Deterministic selection preimage

Pool eligibility is a deterministic function of the frozen availability
manifest. The pool-freeze preimage is exactly (packet §6.3):

```text
schema = g1-successor-pool-v1 (AMENDED 2026-09-04 — supersedes the degenerate
          group-key draft order found by review; amended before gate (f), no
          source access has occurred)
candidate_rows = sorted canonical rows with typed availability outcomes
exclusion_rules = exact hash of the rules below; exclusion inputs (e413
          candidate inventory, e394/T035 sets) hash-bound in the checkpoint
selection_seed = SHA256('g1-successor-pool-v1:seed:' || checkpoint_sha256)
          (external-anchor-bound, not a secret)
unit_of_selection = canonical UniProt base accession (a protein); census rows
          travel with their accession
order_key(accession) = SHA256(selection_seed || '|' || base_accession)
selection = lexical order_key over all eligible accessions; admit in order up
          to the 600-accession cap; ALL census rows of an admitted accession
          are retained (group-atomicity is enforced by folds, not selection)
shortage = typed E420_UNDERPOWERED_POOL
```

The exclusion-rule text, the quota/stratification parameters, and the
deterministic candidate cap are frozen in `docs/g1_successor_analysis_plan.md`
and `docs/g1_successor_source_allowlist.json` before the census runs; their
exact values are OPEN parameters resolved at freeze (§13, O5). A
cap-induced shortage is a typed `E420_UNDERPOWERED_POOL` stop; no post-hoc
enlargement is permitted. Availability responses used before the label phases
are metadata-only and cannot become mapping or label inputs.

### 4.3 Prior-pool exclusions

The pool excludes, by a preregistered comparison manifest frozen in the
checkpoint: all exact **e413/CASP target identifiers** (e413 census retained at
`results/e413/candidate_inventory.json`), all **e394** experimental attempt
identifiers, and all **T035 fixture identifiers**. The normalization
precedent is `docs/e418_prior_pool_sentinel.json`
(SHA `d4968c16f24242dff3add50747bf0dd36bfe870b62899ae0165b44d7b5b62541`);
e420 produces its **own** fresh projection/comparison artifact and does not
reuse e418's. Exclusion is not a reason to substitute a favorable candidate.
Any duplicate accession, ambiguous chain, unresolved group, sequence mismatch,
or provider-release mismatch is a typed exclusion.

### 4.4 Identity and mapping guards (carried from e413, binding)

- Unified mapping-guard regime (amended 2026-09-04; resolves review M3):
  **primary-route rows** — structural guards only (the primary schema carries
  no identity/coverage fields): strict one-to-one ranges, nonempty
  sequence/data, row entity = requested entity; **segment-resolution
  response** (the entity-id resolution fetch for the same entry) — the
  selected entity's mappings must satisfy SIFTS identity **>= 0.90** and
  coverage **>= 0.80**; **segment-as-interpolation-fallback rows** (only when
  the primary route is unavailable) — the strictest regime, SIFTS identity
  **== 100.0** exactly; observed C-alpha count **>= 30** and observed
  amino-acid agreement **>= 0.97** always apply, over retained bytes only.
- The AFDB model sequence must equal the expected canonical UniProt-derived
  sequence after uppercase normalization (no gap or substitution); the AFDB
  `sequence` field obeys the offset fence
  `sequence == canonical[sequenceStart-1:sequenceEnd]` and
  `sequenceChecksum` is the required MD5 `^[0-9a-f]{32}$` (offset context, not
  identity); internal SHA-256 is computed over the parsed CIF sequence
  (dual-checksum rule, contract spec §3).
- Residue mapping is the declared one-to-one path
  `PDB author residue number -> UniProt position` via the primary PDBe
  per-residue route (segment fallback guarded, interpolation-only). A segment
  with numbering mismatch, ambiguous chain, non-one-to-one map, or AA guard
  failure is rejected; no construct-local ATLAS `UnP_num`, pLDDT, B-factor, or
  absorption value repairs a mapping failure.
- The exact route, entity_id resolution route, and literal response schemas
  are bound in `docs/g1_successor_source_allowlist.json` before capture
  (packet §4.2: the registration must bind one exact route and its literal
  response schema). The earlier draft's O3-open wording is historical lineage;
  the current entity-id resolution route is the typed, uniquely matched route
  recorded in the append-only O3 correction at the end of this registration
  and bound in `docs/g1_successor_source_allowlist.json`.

### 4.5 Typed exclusions (no repair)

Excluded or retained as typed `U` without repair: missing/ambiguous accession,
chain, construct, or residue identity; sequence mismatch, non-one-to-one
mapping, failed guard, unresolved author/label numbering; missing experimental
coordinates at a needed residue; a predicted/simulated coordinate used as
reference truth; dynamic evidence that is AF-derived or shares the E
calculation; duplicate accession/protein, repeated construct, or homology-group
collision without a predeclared arbitration rule; source bytes without a
release/version, exact path, raw hash, or license scope; time/template
overlap the registration cannot resolve under the amended AFDB disposition
(§5.1); incomplete, mutated, truncated, or synthetically reconstructed archive
members; values in an ambiguity interval or with conflicting evidence. No
value is imputed, rounded favorably, silently dropped, or relabeled.
Missingness is reported separately from a negative scientific result.

## 5. Independent label channels

### 5.1 E — coordinate error: AF model vs released experimental reference

Predictor side: only the exact frozen AFDB coordinate child for the exact
accession/model entity. Reference side: only the separately released
experimental coordinate source resolved under O2; a predicted or simulated
structure is never truth.

**AF temporal/provenance disposition (amended contract, binding).** The AFDB
API does not document `training_data_cutoff`, `provenance_entry_ids`,
`derived_reference_ids`, or an `api_release` field; the retained OpenAPI
confirms the absence (contract spec §3). e420 therefore binds option (A):
`modelCreatedDate` (strict `strptime("%Y-%m-%dT%H:%M:%SZ")`, not-after the
capture anchor) must be **strictly earlier** than the selected experimental
reference's initial release date; that date is read IN-FILE from the
same-bytes reference mmCIF: the `_pdbx_audit_revision_history` row with
`ordinal_id` 1 (`revision_date`) — missing/ambiguous history is a typed
temporal failure, never an out-of-band fetch (bound 2026-09-04, review
M4/tm4). **Reference selection (bound 2026-09-04, review M5):** the reference
for a pool candidate is the candidate's own census PDB entry (the ATLAS
chain's entry); no alternative structure may be substituted post hoc; `training_data_cutoff` is a **typed
unknown** and is never invented; `af_model_identity =
(modelEntityId, latestVersion/allVersions integers, modelCreatedDate,
digest(template-ID set), digest(protocol method_type set))` with the
template/provenance sets derived only from the same-bytes ModelCIF
dictionary-derived REQUIRED categories (`_ma_model_list`,
`_ma_template_ref_db_details`, `_ma_template_details`, `_ma_protocol_step`).
A floating `latest` response is ineligible; every AFDB row binds the per-entry
pins of the human record §2 at capture. A reference whose initial release does
not strictly post-date `modelCreatedDate` is a typed temporal failure.

**E metric (frozen before any threshold interpretation).** After the predictor
freeze and identity guards: registered global rigid alignment; then, per
`docs/e418_registration.md` §6.5 (restated in
`docs/e419_registration_amendment_2026-08-29.md` §3): binary64 row-vector
Kabsch fit with `math.fsum` centroids, pinned
`numpy.linalg.svd(H, full_matrices=False, compute_uv=True, hermitian=False)`
backend, deterministic sign rule, reflection
`R = U diag(1,1,sign(det(UV.T))) V.T`, no per-residue refit; complete
four-threshold C-alpha lDDT with neighbor set `N_i` = eligible mapped
reference C-alpha residues at reference squared distance strictly < 225.0 A^2
and thresholds t in {0.5, 1.0, 2.0, 4.0} A; a GDT-style chain-level score is
QC only, never a residue metric substitute.

- `E=correct` iff local lDDT-C-alpha `ell_i >= 0.60` **and** aligned C-alpha
  distance `d_i <= 4.0` Angstrom;
- `E=error` iff `ell_i < 0.60` **and** `d_i > 4.0` Angstrom (complementary
  strict inequalities);
- missing, disagreeing, non-aligned, nonfinite, rank/tie-degenerate,
  incomplete-neighbor, or threshold-boundary evidence yields `E=U` with a
  retained typed reason.

pLDDT, B-factors, absorption proxies, sequence-only scores, and AF-derived
features cannot become E evidence. Thresholds, atom selection, alignment, and
missingness rules freeze in `docs/g1_successor_label_rules.json` before the
census; no metric or threshold is tuned to a result.

### 5.2 F — dynamics/flexibility: ATLAS three-replicate RMSF

F uses only the fresh authorized ATLAS analysis archive per selected
`{PDB}_{CHAIN}`, joined through the declared path
`ATLAS PDB_num -> selected-chain SIFTS pdb_num -> SIFTS UniProt_num ->
canonical UniProt sequence position`. The expected ZIP member grammar is
`{pdb_chain}_RMSF.tsv` + `{pdb_chain}_corresp.tsv` with exact headers
(`seq,RMSF_R1,RMSF_R2,RMSF_R3` and
`UnP_seq,PDB_seq,PDB_num,CIF_num,UnP_num,gap,hetatm,corrected,no_fullbb,no_ca`),
equal row counts, and original ordinal identity (no sorting or inferred
pairing); `CIF_num == PDB_num` is a per-row validated invariant. These names
are **project requirements verified against the fresh authorized archive at
capture** (packet §4.4), not assumed contract evidence; a mismatch is a typed
failure, never a silent adaptation.

For numeric RMSF, the registered rule is:

- the mean of **exactly three** finite, nonnegative replicate values, computed
  only as `math.fsum(values) / 3`;
- **high/flexible** if the mean is **>= 1.5** Angstrom;
- **low/ordered** if the mean is **<= 1.0** Angstrom;
- **U** for the open interval (1.0, 1.5), a missing member/value, a mapping
  failure, a nonfinite/negative replicate, a nonzero flag, or an ambiguous
  correspondence.

All five flags (`gap`, `hetatm`, `corrected`, `no_fullbb`, `no_ca`) must be
literal ASCII `0`/`1`. A two-replicate or partial mean is forbidden. The
same-reference missing/construct-span signal is reported only as a sensitivity
role and cannot satisfy primary flexible support. DisProt, PED, NMR, or any
other dynamics source is disabled; enabling one requires a fresh registration
fixing its release, license, endpoint/member schema, and mapping before any
census.

### 5.3 Combined labels

Each residue receives exactly one state recomputed by the constructor from
E and F: `(correct,ordered)=OC`, `(correct,flexible)=A`,
`(error,ordered)=B`, `(error,flexible)=O`; any U endpoint yields `U` with a
retained reason. Callers cannot provide a state, class count, residue count,
or favorable digest.

### 5.4 Gap policy (typed integer rule)

The ATLAS correspondence is rejected as a whole for the profile
(`atlas_correspondence_not_pdb_sifts_uniprot` semantics, contract spec §4)
only for: malformed rows; a PDB number outside the chain range; **zero**
joined rows; or an unmapped fraction violating the strict typed integer rule

```text
5 * unmapped_count > total_count     # integer arithmetic; strictly > 20%
```

SIFTS-absent ATLAS PDB numbers become per-row `U` with
`f_reason=missing_atlas_join` and full no-drop accounting. Gap rows carry
`derivation_source=SIFTS_GAP`, `uniprot_num=null`,
`identity_source=UNRESOLVED_SIFTS_GAP`, and a gap AA cross-check against the
reference atom/canonical AA. Two exact-match declaration channels
(`declared_sifts_absent_pnums`, `declared_atlas_profile_absent_pnums`) are
mandatory; mismatches are typed failures
(`sifts_absent_declaration_mismatch` / `profile_absent_declaration_mismatch`).

## 6. Groups, folds, and support floors

Sequence/group construction and the fold salt are frozen before support is
read (packet §7.6):

```text
pairwise alignment: Smith-Waterman
match=2, mismatch=-1, gap=-2, leftmost tie break
edge: identity > 0.25 and aligned columns >= 100
groups: connected components, ordered by lexical component minimum
fold order: SHA256("g1-successor-fold-v1:" || group_id), lexical hex order
fold assignment (amended 2026-09-04, review m9): WITHIN EACH PRIMARY CLASS
(A, then B): sort that class's groups by fold hex, assign sorted position
i modulo 5 — guarantees the 4+4 per-fold floors whenever the cohort floor
holds; OC and O groups: global hex order, i mod 5 (presence floors only)
```

The exact group-map preimage records all rows, source/documentation hashes,
alignment parameters, rejected/unmapped rows, and the resulting group-map
hash. No group may appear in more than one train/calibration/test partition.
The fold manifest records every candidate, excluded row, group, fold, and
reason.

**Support floors (carried from e413, adjusted per the G1 packet):**

- per-protein primary floor: **10 independently mapped residues for A** and
  **10 for B** (`N_RES_PRIMARY=10`); one-site support is not allowed; counts
  are never pooled across proteins or groups;
- cohort floor: **20 proteins and 20 groups for each primary class** (A and
  B);
- each of the **five folds** must contain at least **four A groups** and
  **four B groups**;
- **OC and O presence:** at least one protein and one group in each, with the
  same per-unit residue floor;
- **cell residue floor:** at least **200 unique non-`U` residue identities in
  each cell**, keyed by `(base_accession, uniprot_num)`.

These floors are the registered minimum sample sizes. A protein is a canonical
UniProt base accession; a group is counted once by its deterministic
`group_id`; repeated structures never inflate a protein and repeated proteins
never inflate a group. Group-atomic member rules (mandatory membership,
`group_unresolved` on any unexplained member) follow
`docs/e418_registration.md` §5.1 as carried context.

## 7. Provenance schema binding and labeled-row construction

Every raw source row carries: role, phase, source/release ID, endpoint or
archive/member path, retrieval UTC, response metadata (status, headers when
available, content type, length), raw byte length and SHA-256,
compressed/decompressed/member hashes as applicable, parser/version/commit,
sequence/mapping hashes, independence ID, and typed status. Raw bytes are the
immutable inputs; canonical parsing is a separate derived artifact.

Each normalized labeled row binds the per-row typed exact-key envelope of the
retained provenance schema (contract spec, per-row envelope update):
`raw_row_sha256`, `normalized_row_sha256`, `source_row_sha256`,
`modelcif_projection_sha256` where applicable, the complete non-empty
`source_member_hashes` map (`pdbe_member`, `atlas_member`, `afdb_child`,
`reference_member`, `sifts_member`), gap-row additions when applicable, and
context `emitted_by`, `release_ids`, `capture_anchor`. The exact schema bytes
(provenance schema generation, currently v27 `a223f271…`) are bound by the
checkpoint; no pointer indirection and no root-only context may replace the
per-row fields. A global hash or synthetic canonical-row digest never
substitutes for a raw member hash.

The preflight/audit must run and record at least these negative checks
(packet §7.7): (1) flip one byte in a copied raw file and verify
hash/manifest rejection; (2) truncate a copied archive/member and verify
length/hash rejection; (3) change a member path or header and verify
schema/manifest rejection; (4) change a release/version scalar and verify
lock rejection; (5) alter the fold/group preimage and verify aggregate-hash
rejection.

The source-free receipt explicitly records
`labels_written=false`, `embeddings_read=false`, `models_run=0`,
`e412_started=false` until the registered label phases; later source receipts
preserve the same registration/checkpoint/authorization hashes and never
overwrite a failed receipt.

## 8. Success criteria, falsifiers, and registered verdicts

**Success (G1 definition; `docs/adcp_roadmap.md` Phase 1 + packet §7.1).** All
of: (i) real ground truth established — the H1 floors of §2.1 met by the
four-cell census built from the two independent channels; (ii) predictor-blind
pool with adequate support — census frozen as `POOL_FREEZE` before any
predictor or label access, floors met without pooling rescues; (iii)
release-pinned provenance audited — every row bound to §3.2 pins (or a
captured, pinned release under the pins-only policy), raw-byte hashes complete,
mutation tests recorded, and an independent audit reproduces the support
census. Only then may main record the G1 admission decision; Phase 2/3 remain
separately registered.

**Falsifiers (any one is terminal for the gate).** Fewer than 20 A or B
proteins or groups; a per-protein/per-group/cell/fold floor failure; absent
OC/A/B/O presence; any predictor-blindness breach (predictor/pLDDT/coordinate/
label value read before `POOL_FREEZE`); E/F channel sharing (dynamic evidence
derived from the E calculation or AF-derived); a SIFTS-release invention or a
mutable `latest` pointer in any receipt; release/pin/license mismatch or
missing raw hash; AFDB temporal failure (`modelCreatedDate` not strictly
before the reference release) or a fabricated cutoff; incomplete gap
declaration or a `5*unmapped > total` violation treated as a repairable
condition; group/fold leak or post-freeze reassignment; provenance-envelope
failure; a failed mutation test; prior-pool overlap; any fetch outside the
registered allowlist.

**Registered verdict vocabulary (e420-owned; declared here, not inherited).**

- `E420_H1_SUPPORT_PASS` — census and provenance satisfy §2.1; enables main's
  separate G1 admission decision and audit; no Phase-2/3 claim.
- `E420_INFEASIBLE_CURRENT_POOL` — lawful census completed but support floors
  unmet; write the negative result, preserve the census, stop.
- `E420_UNDERPOWERED_POOL` — deterministic cap-induced shortage.
- `E420_SOURCE_CONTRACT_FAILURE` — release/pin/license/schema/raw-byte or
  mutation-test failure.
- `E420_INDEPENDENCE_FAILURE` — predictor-blindness or channel-separation
  breach.

Any gate HOLD state (§10) is typed and blocks all later gates. Negative
results are logged, never buried.

## 9. Compute and access budget

- **Hardware envelope (AGENTS §5):** GTX 1050 4GB, 4-core i5, 23GB RAM +
  4GB zram.
- **Memory cap (mandatory):** any single Python job with expected RSS > 2GB
  runs under `systemd-run --user --scope -p MemoryMax=8G` so a runaway dies
  alone. Parsers/grouping/census jobs are CPU-only, single-process or bounded
  workers, resumable, and logged.
- **GPU = ESM-2 inference only, fp32,** models <= t36/150M; iterate with
  t6/8M and t12/35M; reserve t36/150M for final confirmation; pinned **cu126
  wheels — never upgrade to cu130+**; CPU fallback always acceptable; no full
  fine-tuning; LoRA only on t6/t12 if ever. **No GPU and no ESM inference is
  used anywhere in e420's G1 scope** (census, labels, audit are CPU-only);
  any later predictor evaluation is a separately registered phase after a G1
  pass.
- **Access discipline:** public endpoints only; bounded, logged request
  schedule; honor HTTP backoff/status signals; record actual request counts;
  no numeric provider rate/retention terms are documented, so conservative
  low concurrency is mandatory. Probe budgets are frozen in the allowlist
  before dispatch (O5). Cache every download under `data/`; never re-download
  retained bytes; never git-commit large files.
- Statistics, conformal calibration, data pipelines, and audit reruns run on
  CPU.

## 10. Lifecycle gates in order (typed entry/exit; each fails closed)

**(a) Pre-source immutable checkpoint freeze.**
Entry: this registration committed with its exact SHA recorded; companion
checkpoint files drafted and committed — `docs/g1_successor_analysis_plan.md`,
`docs/g1_successor_label_rules.json`, `docs/g1_successor_fold_rules.json`,
`docs/g1_successor_source_allowlist.json`,
`docs/g1_successor_no_source_config.json`,
`docs/g1_successor_dependency_lock.txt`,
`docs/g1_successor_source_authorization.json`,
`docs/g1_successor_external_anchor.json` (packet §6.1 set; the authorization
JSON must materialize the human record's seven policy values into the packet
§5 field set). Exit: `docs/g1_successor_pre_source_checkpoint.json` created
fresh — canonical UTF-8 JSON, sorted keys, compact separators, no floating
timestamps, no self-referential hash — binding exact relative path, byte
length, and SHA-256 for every file above plus the repository commit, the
e-number, the documentation snapshot hashes (§3.2), the legal-authorization
hash, the pool-selection preimage, the label-rule hash, and the fold-manifest
hash; digest computed as
`SHA256(UTF8(json.dumps(payload_without_self_hash, sort_keys=True, separators=(",", ":"))))`;
produced in a detached clean checkout and read back byte-for-byte. The packet
itself is not part of the manifest and never enters its own hash preimage.
Failure state: typed `CHECKPOINT_HOLD`; no source access.

**(b) Five-role pre-access reviews on exact committed bytes.**
Entry: the exact committed tuple (registration + companion files + checkpoint
candidates) exists. Exit: five same-SHA reviews — lab, data, theory, audit,
lead — recorded against those exact bytes; in single-session mode main
performs these as sequential, fresh-task-context role passes with independent
recomputation where possible; a self-read or a green command is not an
independent audit. Failure: typed `REVIEW_HOLD`.

**(c) Detached audit from clean checkout.**
Entry: reviews complete. Exit: a detached audit from a clean checkout
independently reproduces the no-source self-tests, exact hashes, typed
boundaries, and checkpoint bindings. Failure: typed `AUDIT_FAIL`; repair and
re-freeze; no source access.

**(d) External anchor (OpenTimestamps; main-executed, no human step).**
Entry: checkpoint digest, registration/analysis-plan SHAs, pool-selection
algorithm + seed/preimage, E/F thresholds, four cells, U handling, support
floors, and five-fold design prepared as the registration snapshot file.
Exit: the snapshot bundle is SHA-256-stamped via OpenTimestamps and the
`.ots` proof committed (mechanism per DECISIONS 2026-09-04T13:05Z, commit
`b1c12d62`; the gate-C authorization record is already stamped, proof
`cb6210c97f990c0553df07d23771a56e733db3d65b30086868490a9da5b6e2d0`). The
proof is PENDING until Bitcoin confirmation (~hours; automatic), then
verified with `ots verify` and the confirmation recorded in
`docs/DECISIONS.md`. No account is required and document content never
leaves the machine — only the hash is anchored; this also removes the
self-registered-registrar limitation. OSF/Zenodo remain optional later
upgrades if a reviewer-facing registry DOI is ever wanted. Typed state:
`ANCHOR_PENDING` until BTC confirmation, then `ANCHOR_BOUND`.

**(e) Bounded gate-C capture dispatch.**
Entry: (a)–(d) complete and main records the explicit dispatch. Order:
(i) **grammar validation first** — offline validation of every parser/adapter
against the retained OpenAPIs, census ZIP, provenance schema, and ModelCIF
allowlists, with the mutation matrix of §7; (ii) **small samples** — bounded
per-source probes (ATLAS per-chain metadata/analysis, PDBe mapping, AFDB
metadata; exact probe budgets frozen in the allowlist per O5) with full
receipts (URL, method, status, headers, bytes, hashes, UTC); (iii) full
census acquisition. Exit: complete capture receipts; every success and typed
failure retained; request counts logged. Failure: typed
`DISPATCH_HOLD`/`E420_SOURCE_CONTRACT_FAILURE`.

**(f) Pool construction.**
Entry: dispatch receipts complete. Exit: predictor-blind availability census
completed per §4.1; inventory, accession/sequence identity guards, homology
groups, and folds frozen and hash-sealed as `POOL_FREEZE` **before any
predictor value, coordinate, reference byte, or label is opened**; prior-pool
comparison manifest empty-intersection receipt. Failure: typed no-go per §8.

**(g) Label construction.**
Entry: `POOL_FREEZE` sealed. Chronology (packet §6.3): freeze the AFDB release
and exact model IDs/versions for the frozen candidates (no floating latest);
open AF predicted coordinates/metadata under the §5.1 temporal/provenance
disposition; then open the separately released experimental reference/mapping
artifacts under the declared chronology; then open the independent ATLAS
dynamics archives; construct the four-cell ground-truth census. Exit:
labeled-row ledger with §7 provenance binding; every E/F/U reason retained.
Failure: typed `E420_SOURCE_CONTRACT_FAILURE` / `E420_INDEPENDENCE_FAILURE`.

**(h) G1 evaluation and independent audit.**
Entry: label ledger complete. Exit: H1 census per §2.1/§6 floors; mutation
receipts; source-free flags verified; independent audit (fresh task
context/process/path, independent recomputation of load-bearing numbers,
starting from recorded artifacts) returns PASS / PASS-with-notes / FAIL; only
after audit PASS may main record the G1 admission decision in DECISIONS and
board. On PASS: Phase 2/3 require separate registrations. On floor failure:
`E420_INFEASIBLE_CURRENT_POOL`, negative result logged, stop before any
embedding, model, conformal, or Phase 2/3 analysis.

## 11. Fences (hard rules)

1. **No raw redistribution** of source bytes by any channel (repo, paper
   supplements, OSF, GitHub, or otherwise).
2. **No commercial use of ATLAS-derived data.** ATLAS dataset is CC-BY-NC-4.0;
   use is limited to the non-commercial research line and academic
   publication; any commercial use requires a separate ATLAS license from
   DSIMB or exclusion of ATLAS-derived assets. License inheritance:
   AFDB-derived data CC-BY-4.0; attribution and required citations in any
   derived release.
3. **No capture outside registered lanes;** no capture before gates (a)–(e)
   complete; no endpoint, member, or source role not in the frozen allowlist;
   ATLAS POST search routes forbidden; no credentials, ever.
4. **Typed failures only.** No silent imputation, favorable rounding, silent
   dropping, relabeled `U`, or retry-until-favorable; failed rows and failed
   receipts are retained; negative results are logged, never buried.
5. **No mutable `latest` pointer** in any ledger, receipt, result, or claim;
   every "latest" resolves once at capture into concrete bytes + SHA-256.
   The SIFTS release identifier stays a typed unknown; inventing a release ID
   is a falsifier.
6. **Predictor blindness is inviolable:** no pLDDT, coordinate, template
   value, provenance value, RMSF value, B-factor, embedding, or model output
   influences pool membership, groups, folds, thresholds, or the census before
   `POOL_FREEZE`.
7. **Write-scope:** in single-session mode main is the sole writer/committer;
   large raw files are never committed; raw caches live under `data/`.
8. **Honesty:** missing data is reported as missing; a noted limitation beats
   a fabricated capability; no claim enters a paper before a reviewer audit
   PASS.

## 12. Required artifacts (all under registered paths)

- `results/e420/candidate_inventory.json` — frozen universe rows from the
  census member, typed availability outcomes, inclusion/exclusion reasons,
  pool preimage and SHA-256, predictor-blind lock timestamp.
- `results/e420/source_availability_manifest.json` — per-row role/phase/
  endpoint/member, raw bytes + hashes, release pins, typed failures.
- `results/e420/prior_pool_comparison.json` — e413/e394/T035 exclusion
  projection with empty-intersection receipt.
- `results/e420/pool_freeze.json` — sealed inventory/group/fold hashes
  (`POOL_FREEZE`).
- `results/e420/group_map.json`, `results/e420/fold_manifest.json` — preimage,
  parameters, per-fold membership.
- `results/e420/label_ledger.json` — per-residue E/F/U with the §7 envelope;
  `results/e420/four_cell_census.json` — cells, floors, fold support, verdict.
- `results/e420/mutation_receipt.json`, `results/e420/run_receipt.json` —
  negative checks; registration/commit, commands, dependency versions,
  request counts, memory policy, `labels_written` / `embeddings_read` /
  `models_run` / `e412_started` flags.
- Checkpoint set per §10(a); five review records; detached audit outputs;
  external-anchor receipt (`docs/g1_successor_anchor_receipt.json`, binding
  the checkpoint `.ots` sha and the transitive-binding claim); capture
  receipts.

Every artifact self-hash uses the canonical compact sorted-key UTF-8 JSON
rule excluding its own hash field.

## 13. OPEN items (explicit; each blocks its named gate)

- **O1 — RESOLVED 2026-09-04 (main review): committed as
  `docs/e420_registration.md`;** the checkpoint binds this exact path+SHA and
  records the packet §6.1 slot-name difference
  (`g1_successor_registration.md`) as a superseding naming decision
  (DECISIONS 2026-09-04).
- **O2 — Experimental reference-coordinate endpoint family (gates a/e/g;
  blocks E entirely).** The three authorized sources, as literally scoped in
  the human record, name PDBe SIFTS mappings, AFDB prediction data, and ATLAS
  analysis data — none is a documented experimental-coordinate download
  contract, and the retained evidence records the PDBe coordinate-download
  contract as INCOMPLETE (`docs/data_source_contracts_20260830.md` addendum:
  `/pdbe/entry-files/` 404 negative evidence). Resolution path A: a recorded
  scope clarification that PDBe-served entry coordinate files are within the
  authorized PDBe scope, plus completion of that endpoint contract capture
  before checkpoint freeze. Resolution path B: a new human authorization
  adding a coordinate source (e.g., RCSB files/entry endpoints) with license
  binding (WWPDB policy). Until resolved and frozen in
  `docs/g1_successor_source_allowlist.json`, E labels cannot lawfully be
  constructed and gate (g) cannot run.

  **RESOLVED 2026-09-04 via Path A:** the human confirmed PDBe-served
  experimental coordinate entry files are within the authorized PDBe scope
  (answer "yes", recorded in
  `docs/human_gateC_authorization_20260904_addendum_o2.md`, sha256
  `169032e8b35f7f1fcb0c63b50c08d5bead074931afae65e53a9d78607cfebcf0`, OTS
  proof `31261940…`). Documentation-class contract capture of the
  entry-files endpoint family retained at
  `data/source_contracts/pdbe_entry_files_contract_capture_20260904.json`
  (`ddfa3e2338b2b5ed096bc5568133b985a4bdc0ff978a0b820fcd310f7c8d4783`):
  `.cif` 200 OK (Range ignored on this path), ranged reads 206 on legacy
  `.ent`, directory listing 404, wwPDB public-domain data-owner statement
  captured. Final endpoint
  grammar verification against fresh authorized bytes remains a gate-(e)
  typed requirement.
- **O3 — RESOLVED 2026-09-04:** identity guards close over retained
  PDBe-mapping + AFDB-sequence (`uniprotSequence`) + ATLAS-corresp bytes; no
  canonical FASTA fetch is needed or permitted. The entity_id resolution route
  is frozen in the allowlist (segment mapping response; typed ambiguity
  rejection; segment never feeds the flat parser). The historical draft-open
  wording is retained in the review records as provenance.
- **O4 — RESOLVED 2026-09-04: main bound receipt v54**
  (`bfcca135256038265450cb8c045d1959b8abc29b349bafea0316fcdb501475f3`,
  69,564 B) + provenance schema v27
  (`a223f271162fa921cffc23d4b040a603403b76df288cbf439021984df6b14efe`,
  24,213 B) as the e420 checkpoint pair (DECISIONS 2026-09-04; both files
  committed). The human record's v30 reference remains historically accurate
  for the authorization record; the eight artifact pins are byte-identical
  across generations. The e418 theory spec's v52/v24 binding is historical to
  that lane and does not transfer (O10).
- **O5 — Quota/stratification, candidate cap, and probe budgets (gates
  e/f).** Mechanism fixed (§4.2); exact values must be frozen in
  `docs/g1_successor_analysis_plan.md` / `docs/g1_successor_source_allowlist.json`
  before the census; not invented in this draft.
- **O6 — RESOLVED 2026-09-04: anchor mechanism = OpenTimestamps** (DECISIONS
  2026-09-04T13:05Z, commit `b1c12d62`); no human step; the checkpoint
  snapshot is stamped at gate (d); BTC confirmation pending and verified
  later. OSF/Zenodo = optional upgrades only.
- **O7 — RESOLVED 2026-09-04: `docs/g1_successor_source_authorization.json`
  materialized from the human record into the packet §5 field set**
  (reference-coordinate scope field marked
  `UNRESOLVED_PENDING_HUMAN_SCOPE_CONFIRMATION`); the independent auditor
  reads it back at gate (c).
- **O8 — RESOLVED 2026-09-04: DECISIONS entry appended** recording main's
  assignment of e420 to this G1 successor and the retirement of the informal
  "e420/T-038" Phase-3 label (§0.3); the Phase-3 method receives a fresh
  allocation only after G1.
- **O9 — ATLAS analysis-ZIP member/headers verification (gate e).** The
  OpenAPI does not document member names/headers; the frozen grammar
  (§5.2) is a project requirement to verify against the fresh authorized
  archive at first capture; a mismatch is a typed failure requiring
  amendment, never adaptation.
- **O10 — Non-transferability confirmations (standing).** The e418/e419
  gate-C amendments remain PROPOSED_NOT_EFFECTIVE; e420 inherits no tuple,
  review, checkpoint, or audit from any prior lane; each confirmation is
  re-recorded in the e420 review records.

## 14. Registration and approval log

- **2026-09-04 —** e420 allocated by `main` (single-session mode) to unblock
  G1 per `docs/G1_UNBLOCK_PACKET.md` §10 blocker 3. This DRAFT was authored
  from local repository records only; **no network, source, API, cache,
  candidate, sequence, coordinate, label, embedding, or model access
  occurred.** All eight authorization pins were verified against on-disk
  bytes (`data/source_contracts/`); the receipt v54 / provenance-schema v27
  on-disk state versus the human record's v30 reference is recorded as O4.
  Status remains DRAFT until §10 gates (a)–(c) complete on the committed
  bytes.

- **2026-09-04 (main review) —** Draft reviewed in full by main. Resolved at
  review: O1 (committed as `docs/e420_registration.md`), O4 (checkpoint pair
  receipt v54 `bfcca135…` + provenance schema v27 `a223f271…`, both files
  committed), O6 (anchor = OpenTimestamps, commit `b1c12d62`; no human step),
  O7 (`docs/g1_successor_source_authorization.json` materialized;
  reference-coordinate scope field UNRESOLVED pending human confirmation),
  O8 (supersession recorded in DECISIONS). Still open: O2 (one-line human
  scope confirmation on PDBe entry coordinate files), O3/O5 (exact routes,
  quotas, caps, probe budgets frozen in
  `docs/g1_successor_source_allowlist.json` / `docs/g1_successor_analysis_plan.md`
  at checkpoint freeze), O9/O10 (verified at their own gates).

- **2026-09-04 (O2 resolution) —** Human confirmed PDBe coordinate-file scope
  ("yes"); addendum recorded and stamped; contract capture retained; O2
  closed via Path A. No other scope change.

```text
E420_STATUS=FINALIZED_REGISTRATION_COMMITTED
SOURCE_OPENED=NONE
PREDICTOR_BLIND_POOL=NOT_CONSTRUCTED
CHECKPOINT=ABSENT
ANCHOR=OPENTIMESTAMPS_PENDING_BTC_CONFIRMATION
G1=NOT_MET
HOLD / NOT_GRANTED / NON-FINAL NO-GO
```

## 15. Append-only O3 resolution correction (2026-09-04)

This section is the current O3 disposition and supersedes the historical
open-question wording retained in §4.4 and §13 above. It is a governance
correction only; it does not open a source, create a candidate pool, or grant
capture authority.

- `rest.uniprot.org` remains outside the three-source scope. Identity guards
  use only the retained PDBe-mapping, AFDB-sequence, and ATLAS-correspondence
  bytes named by this registration.
- For the primary route
  `GET /pdbe/api/v2/pdb/entry/uniprot_mapping/{pdb_id}/{entity_id}`, resolve
  `entity_id` from the retained segment envelope
  `GET /pdbe/api/v2/mappings/uniprot/{pdb_id}` by grouping returned mappings
  per entity and selecting the unique entity whose accession set contains the
  candidate's canonical UniProt accession.
- Zero matching entities, multiple matching entities, a missing accession, or
  any segment envelope/schema failure is a typed rejection with no arbitration
  or fallback. The segment response is used only for entity resolution and the
  guarded interpolation-only fallback; it never feeds the flat per-residue
  parser. The selected primary mapping response remains the binding label-join
  route.
- This resolution is bound by the `entity_id_resolution_route` object and the
  current `open_items.O3` disposition in
  `docs/g1_successor_source_allowlist.json`. A fresh checkpoint and fresh
  review are required after this correction; the prior checkpoint is not
  capture-authorizing.

```text
O3=RESOLVED_2026-09-04_ENTITY_ROUTE_BOUND
SOURCE_OPENED=NONE
AUTHORITY=NOT_GRANTED
G1=NOT_MET
HOLD / NOT_GRANTED / NON-FINAL NO-GO
```


## 15. Amendment 2026-09-04 — phase-(ii) round-1 findings: ATLAS member-superset + metadata-superset policies

**Trigger (typed findings, phase (ii) round 1, 10/16 probes executed, HALT per
§10(e) stop conditions; receipts `results/e420/probes_receipt.json`, raw bytes
`data/e420/probes/`, request log in the receipt):**

1. `atlas_archive_member_set_invalid` on BOTH analysis probes: fresh
   `/ATLAS/analysis/{pdb_chain}` archives (16pk_A 93,614,576 B; 1b2s_E
   21,671,996 B) contain the two registered F-channel members **plus 14
   additional members** (pLDDT/Bfactor/Neq/RMSD/gyrate/contacts tsv, R1-R3
   xtc/tpr, .pdb, README.txt). The registered exact-2-member set does not
   match current provider reality; the registered members themselves
   **validated perfectly** (headers/row-shapes §5.2, equal counts 415/415 and
   90/90) when extracted from the retained probe bytes.
2. `atlas_metadata_keys_invalid` on all 3 metadata probes: real responses
   carry a large superset (CATH_*, ECOD_F, ...) beyond the retained OpenAPI
   contract's documented fields; the wrapper shape `{pdb_chain: {...}}` and
   identity discipline matched exactly.
3. `reference_cif_token_invalid: loop arity` on 1crn.cif: the mmCIF loop
   reader did not stop loop-value collection at scalar tags — a phase-(i)
   code defect (fixed; real-bytes fixture added).
4. Probe-runner shape bugs (segment accession extraction) — runner-side only,
   fixed; not a provider issue.

**Amended policy (binding for e420):**

- **ATLAS analysis member superset:** the archive MUST contain the two
  registered F-channel members (`{pdb_chain}_RMSF.tsv`,
  `{pdb_chain}_corresp.tsv`), which are validated per §5.2. Additional
  members are permitted: they are hash-recorded in the outcome
  (`extra_member_hashes`) and **never opened** (trajectories are the bulk);
  a MISSING required member remains
  `ATLAS_ARCHIVE_MEMBER_SET_INVALID`. RMSF values are still never read
  before POOL_FREEZE.
- **ATLAS metadata superset:** required contract fields must be present and
  type-checked; unknown fields are recorded (`unknown_fields_recorded`) and
  are non-fatal. Wrapper `{pdb_chain: {...}}` and identity discipline
  unchanged.
- The validator implements both policies (commit with this amendment);
  tests encode them (mutation (e) rewritten; missing-required still
  rejects; real-bytes 1crn fixture added to tests/fixtures/).

**Effect on gates:** phase (i) re-validated (159+23 green); phase (ii)
re-probe follows (round 1 consumed 10 of 16; round 2 re-runs the probe set
within a fresh <=16 budget; cumulative request count disclosed in receipts).
No labels, no predictor values, no pool effects — probes are identity/grammar
only. This amendment is pre-census, therefore lawful under §2.2.

```text
AMENDMENT=PHASE_II_ROUND1_2026-09-04
HALT_RESOLVED_BY=AMENDMENT_NOT_ADAPTATION


### 15.1 Amendment extension 2026-09-04 (phase-(ii) rounds 2-5): PDBe identifier semantics + probe disclosure

- PDBe segment per-accession object's `"identifier"` field is the **UniProt
  entry name** (e.g. `PGKC_TRYBB`), NOT the accession; the accession is the
  map key under `"UniProt"`. Validator enforces presence/nonemptiness only.
- Phase-(ii) probe disclosure: round 1 (10 requests, pre-amendment typed
  findings — preserved in lead.md 19:22-era notes), rounds 2-5 (10 each,
  cumulative 50 requests total across amendment rounds) — every request
  receipted; final round 5: **10/10 PASS, 0 typed failures**. The <=16-per-run
  budget holds; cumulative counts are disclosed here for honesty.
- Round 3 finding retained: ATLAS corresp author-numbering offset policy
  (constant PDB_num - CIF_num offset per chain; recorded per row-set).


### 15.2 Amendment extension 2026-09-04 (census round-1 findings, 1,531/1,938 rows executed before policy halt):

Census round 1 failure distribution drove four validator policy corrections
(all pre-POOL_FREEZE, no labels affected; raw bytes retained; partial manifest
preserved as results/e420/source_availability_manifest_round1_partial.json):

1. **Non-standard AA3 codes (28 rows: pdbe_primary_code_invalid, e.g.
   startCode='MSE'):** SeMet entries are overrepresented in the ATLAS census.
   Policy: unknown 3-letter codes are RECORDED
   (`nonstandard_code_segments_recorded`), mapped to null 1-letter, non-fatal.
2. **AFDB sunset fields (156 rows: afdb_field_sunset, e.g. entryId,
   uniprotStart/End, paeImageUrl):** sunset fields are documented-but-
   deprecated in the retained contract table. Policy: RECORDED
   (`sunset_fields_recorded`), non-fatal; contract-UNKNOWN fields remain fatal.
3. **AFDB multi-object responses (7 rows: afdb_metadata_cardinality_invalid,
   e.g. objects=6):** real /prediction responses may return one object per
   UniProt sequence. Policy: every object validated; census records
   `objects_summary`; label-phase selection is deterministic — earliest
   modelCreatedDate, then lexicographic modelEntityId.
4. **PDBe segment real nesting + identifier semantics:** response nests under
   `"UniProt"`; per-accession objects carry `{name, mappings}` (no
   identifier); `identifier` is the UniProt ENTRY NAME. Policy: unwrap the
   "UniProt" level; required keys {name, mappings}; extras recorded;
   identifier presence-only. Corresp numbering: constant per-chain offset
   `PDB_num - CIF_num` (e.g. 16pk_A +4) is valid; varying offsets remain
   `ATLAS_CORRESP_CIF_PDB_MISMATCH`.

Entity-ambiguity typed exclusions (5 rows) are legitimate availability
findings and remain exclusions.


### 15.3 Amendment extension 2026-09-04 (independent-review reconciliation: anchor policy v2, cap_bound rule, overflow measurement)

Trigger: independent review (codex-c) PASS-with-notes on the G1 negative
verdict chain with three authorization/evidence issues. Repairs:

1. **Anchor policy v2 (resolves issue 1):** the gate anchor for a checkpoint
   re-freeze lineage is satisfied by (a) the LINEAGE anchors BTC-confirmed and
   block-merkle-verified (record proof @965448, addendum proof @965447,
   a0af0425 rebind proof @965451) PLUS (b) the current freeze's own fresh
   stamp in ANCHOR_PENDING state, passively confirmed and block-verified when
   the calendar lands (recorded in the anchor receipt). A re-freeze resets the
   confirmation clock by design; it does not reopen the gate. The current v6
   stamp (770 B, pending at amendment time) is covered by (b).

2. **cap_bound rule (resolves issue 2):** `four_cell_census.cap_bound` MUST
   derive from the bound pool_freeze (`selection.cap_bound`, true when
   admitted == cap), never hardcoded. The label runner's hardcoded False is a
   bookkeeping bug repaired in code. Verdict-mapping refinement: when floors
   fail AND cap_bound=true, the mechanical verdict is
   `E420_UNDERPOWERED_POOL`; when the recorded per-accession temporal-failure
   rate is 100% across all measured rows, the artifact records
   `temporal_cause=AFDB_TEMPORAL` alongside, and INFEASIBLE remains
   defensible as outcome-equivalent (independent review concurs). Both
   labels are recorded in the census artifact to remove ambiguity.

3. **Overflow measurement (resolves issue 3):** a BOUNDED temporal-check-only
   capture over the 269 overflow accessions (canonical sequences not admitted
   by the seeded selection) is AUTHORIZED: per overflow CENSUS ENTRY, exactly
   one reference cif GET (entry-files lane, ~70-500 KB) to extract the
   initial-release date; the modelCreatedDate comes from the ALREADY-CAPTURED
   census AFDB raw bytes (no new AFDB fetch). Schedule/receipts per analysis
   plan. Purpose: replace the cap-not-causal inference with measured
   evidence, mechanically resolving INFEASIBLE vs UNDERPOWERED.
