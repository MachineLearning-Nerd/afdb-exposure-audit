# e422 registration (FROZEN) — per-residue conformal model-quality estimation on leakage-free post-snapshot PDB releases

> **Status: FROZEN (gate (a), 2026-09-09; T0 = 2026-09-11T00:00:00Z; see
> the §6 gate-(a) FREEZE block).** Draft prepared 2026-09-09 by
> the e421 execution session (w1:pG) following main's strategic decision to
> draft in parallel with e421-P (DECISIONS 2026-09-09); iterated to
> PASS-with-notes through the codex-main review cycle, rounds 4–11
> (docs/audit_e422_review_cycle_20260909.md). **E-number "e422": ALLOCATED
> 2026-09-09 under delegated authority** (AGENTS.md §4 places e-number
> allocation with `main` as sole allocator; here allocation was exercised
> under main's own recorded delegation — main delegated E-number allocation and
> every gate/dispatch decision for THIS registration to w1:pG — DECISIONS
> 2026-09-09, "e422 dispatch authorization and decision delegation"). Under that
> delegation the lifecycle proceeds: five-role review (b) → detached audit
> (c) → registration checkpoint freeze (a) → OTS anchor (d) → bounded
> dispatch (e), each gate committed and ledgered before the next.
> **Predecessors:** e420 (E420_INFEASIBLE_CURRENT_POOL — temporal-leakage
> audit), e421 exploratory (E421_UNDERPOWERED_POOL — binary four-cell estimand
> starves the B cell), e421-P (four-cell prospective draft — remains a
> separate, independently-decided registration).
> **Everything in §2 carries e421-P §§2–3 forward AS AMENDED BY §2 ITSELF
> (T0 definition and presence-manifest mechanism, gate-(b) disposition
> 2026-09-09):** frozen AFDB v6 T0-snapshot semantics, rolling weekly
> batches, X-ray ≤2.5 Å release>T0 predicate, fail-closed temporal rule,
> sources, fences, typed failures.
>
> **Provenance caveat:** all e421 figures cited here are from the
> EXPLORATORY_NOT_FULL_REGISTERED run — motivation and design inputs only,
> no evidential weight; every claim rests solely on the prospective pool.

## 1. Motivation and estimand split

Two registered findings motivate a new estimand:

1. **Scarcity is an artifact of the binary endpoint.** e421's four-cell
   design starves the binding B cell (~0.9 residues/week → ≈4.3-year floor
   horizon). The exploratory data show the scarcity is in the
   *discretization*, not the data: lDDT-Cα exists for 24,456 of 24,457
   residue evaluations, and the A cell alone accumulates ~3.1k residues/week.
2. **The live phenomenon is conditional, not marginal.** On the exploratory
   ledger (development-only), the committed classwise AD-CP construction
   on the BINARY correct/error endpoint (pLDDT-derived scores;
   results/e421/adcp_conformal_analysis.json) achieves near-nominal
   MARGINAL coverage (0.896 / 0.946 at 90/95% — accession-level
   protein-disjoint split, R = 200 means)
   while per-pLDDT-bin CONDITIONAL coverage collapses mid-range (five of
   six bins, 0.08–0.66 vs 0.90 — canonical RESIDUE-level split, seed 0;
   the two number families come from different exploratory protocols and
   are tagged here so neither is read as one study's result). The
   registered question TRANSFERS this marginal-OK / conditional-collapsed
   structure to the primary continuous estimand (split-conformal
   intervals for lDDT-Cα on the same feature): whether conditional
   validity fails there too — and which conditioning repairs it — is
   what e422 tests.

**Estimand split (E422-2 resolution).** The PRIMARY estimand is a
**deployable per-residue quality estimator**: distribution-free intervals
for continuous lDDT-Cα built from **AF-observable inputs only** — within
this program's registered channels, pLDDT (the per-residue AlphaFold
confidence output) is the sole AF-observable feature; reference-derived
quantities (B-factor robust-z, |AF−reference| distance, lDDT itself) are
label/evaluation channels and are NEVER features. This is what makes the
result consumable by downstream users of AlphaFold models, who have no
reference structure. The SECONDARY estimand is a **reference-conditioned
diagnostic analysis**: conditional validity evaluated on strata that use
reference data (pLDDT × B-factor-z), reported as analysis, never as a
deployable predictor. CQR/Mondrian below are specified so no
reference-derived quantity enters any interval's construction.

## 2. Design (by reference to e421-P, as amended here)

Rolling prospective census as e421-P §§2–3 **as amended by the T0 /
presence-manifest block below**: T0-snapshot AFDB v6 source (incremental
presence manifest per accession, below; absent accessions →
`ABSENT_FROM_SNAPSHOT`); UTC batch windows (T0 + 7(k−1) days, T0 + 7k
days] — batch 1 = (T0, T0+7 days]; identical
admissibility predicate (including e421-P §2.1's date/instant
normalization — a date-only release date D parses to the UTC instant
`D T00:00:00Z`; a release dated T0's UTC calendar day fails
`release > T0` for any T0 time-of-day and is `PRE_SNAPSHOT_RELEASE`,
fail-closed); `PRE_SNAPSHOT_RELEASE` typed exclusions; identical
sources, rate limits (e421-P §7: ≤2 workers, ≥250 ms inter-request delay
+ jitter, ×2 backoff, max 3 attempts), and receipts (allowlist binding:
docs/e421_source_allowlist.json, SHA-256 recorded at the gate-(a)
freeze).

**T0 and the snapshot presence-manifest (e422 amendment; gate-(b)
disposition, 2026-09-09 — supersedes the carried e421-P §2.1 T0
definition).** For e422, **T0 is the exact second-level UTC instant
bound in the frozen text at gate (a)** — this AMENDS e421-P §2.1's
definition of T0 as the manifest-freeze instant, and correspondingly
supersedes e421-P §2.2 point 1's "observed at T0" observability
argument. The complete AFDB-v6 accession snapshot is not enumerable over
the registered allowlist lanes (no bulk index lane is allowlisted; the
per-accession metadata lane cannot enumerate the ~2×10⁸-entry release
within the registered rate limits), so the presence manifest is built
**incrementally at each accession's FIRST appearance** in a batch, with
the temporal rule enforced **fail-closed per accession at build time:**
the metadata response's date-only `modelCreatedDate` D maps to the UTC
instant `D T00:00:00Z`; `D ≥ T0.date()` (UTC) — **including T0's own UTC
calendar day** — is the typed failure `PRE_SNAPSHOT_RELEASE` (exclusion;
no substitution), because a model created any time during T0's UTC day
may postdate the T0 instant; an accession whose metadata fetch fails or
carries no v6 model is `ABSENT_FROM_SNAPSHOT` (fail-closed). The
manifest records `{accession, v6_cif_url, model_created_date,
afdb_version, model_entity_id}` per accession (multiple v6 records bind
the registered selection rule: earliest `modelCreatedDate`, then
lexicographic `modelEntityId`); manifest content is checkpointed per
batch as canonical-JSON (`sort_keys=True, separators=(',',':')`)
SHA-256 and OTS-anchored with the batch checkpoint (§5/§6). At the
first child fetch, `modelCreatedDate` and `afdb_version` are verified
against the manifest — disagreement or a missing/unparseable created
date is a typed failure (`MODEL_CREATED_DATE_MISMATCH`,
`MODEL_CREATED_DATE_MISSING`; version disagreement →
`AFDB_VERSION_MISMATCH`, source-contract class) — and a child sha256
that changes across fetches is `SOURCE_CONTRACT_FAILURE` (e421-P §2.2
point 4 fence, unchanged). **Code-prefix note (gate-(c) disposition,
2026-09-09):** failure codes emitted by the committed implementation
carry the `E422_` prefix; names in this text are unprefixed
(`PRE_SNAPSHOT_RELEASE` ≡ `E422_PRE_SNAPSHOT_RELEASE`, etc.).
**Registered residual (disclosed, not
assumed away):** first-observation pinning bounds what is observable —
the operating guarantee is "**no model whose creation date falls on or
after T0's UTC day ever enters the ledger**" (fail-closed,
per-accession verified), not "byte-identical to the T0 release state";
an AFDB-side replacement of a model between T0 and its first fetch
without a creation-date bump is outside registered observability and is
disclosed as source-contract risk.

**Label payload per residue:** each
admitted residue records `(lDDT-Cα, d_Kabsch, pLDDT, B-factor robust-z)` —
four continuous channels. Channel construction is e421-P §4 carried
forward unchanged and bound at the gate-(a) freeze to
docs/e421_label_rules.json (its SHA-256 recorded in the freeze): Kabsch
d (fsum centroids, pinned SVD, sign rule), lDDT-Cα (15 Å neighbor
cutoff, squared distance strictly < 225 Å², thresholds 0.5/1.0/2.0/4.0
Å), B-factor robust-z (median/IQR), pLDDT (per-residue B-factor column
of the AF v6 child). **pLDDT parsing is text-registered here, not only
via the label-rules binding (gate-(c) disposition, 2026-09-09):** pLDDT
is the value in the residue's B-factor column of the AFDB v6 CIF, taken
as-is on the 0–100 scale (AFDB stores pLDDT×100 there), matched to the
aligned reference residue by the mapping; no transformation, clipping,
or rescoring is applied.
Channel roles are frozen: **pLDDT = feature
(AF-observable)**; **lDDT-Cα = primary label**; d_Kabsch = descriptive;
**B-factor robust-z = evaluation-only stratification** (never a feature —
E422-2). **No four-cell discretization enters analyses A, B, or D, or
any interval construction; the sole registered exception is secondary
analysis C, which computes e421-P §11's binary E-channel endpoint
(correct = lDDT-Cα ≥ 0.60 ∧ d ≤ 4.0; error = lDDT-Cα < 0.60 ∧ d > 4.0;
class-evaluable rows only, per e421-P §11's own frozen design;
e-ambiguous rows are typed and counted, never silently dropped).**

**Row provenance and chain selection (e422 amendment; gate-(c)
disposition, 2026-09-09 — the registration was previously silent on
chain selection; flagged by the runner implementation).** Each batch
enumerates entries via the census lane, fetches the whole-entry PDBe
UniProt mapping, and processes **EVERY protein chain mapping to a
manifest-admitted accession** — the unit of label construction is the
(entry, chain, accession) triple, matching the exploratory yield basis
(~250 rows/entry over all mapped chains). A triple contributes rows
when ≥ 30 residues align under the mapping (the e421 threshold); each
UniProt position contributes at most one row per triple (an author
residue colliding with an already-matched position is skipped, first
occurrence in reference-file order). All rows of all triples enter the
cumulative label ledger; per-triple typed failures
(`REF_INSUFFICIENT_CA`, `INSUFFICIENT_ALIGNED`, `KABSCH_FAILED`, …)
are tallied in the batch checkpoint, never silently dropped. A UniProt
accession reached through several entries or chains contributes rows
from EACH route — the accession, not the entry, is the analysis unit
(§3's permutation splits accessions). Committed implementation:
experiments/e422_runner.py.

## 3. Pre-registered analyses (executable precision)

**Frozen protocol (analyses A, B, D; analysis C executes e421-P §11's own
frozen design unchanged — e422 does not restate or reuse its 50/50
draw):** R = 200 repeats, seeds 0–199 via `np.random.default_rng(seed)`;
bins = the committed schema-v2 SIX-bin grid, half-open [lo, hi):
[0,50), [50,60), [60,70), [70,80), [80,90), [90,100.01) — the 100.01
top edge is inherited verbatim from the committed e421 schema-v2 grid
(experiments/e421_conformal.py `PLDDT_BINS`) so the half-open convention
still admits exact-100 values; the VALUES LISTED HERE are authoritative
for this registration (the code file is provenance — the detached audit
binds the grid source blob current at gate (c) freeze); levels
0.90/0.95 (analysis A runs both and always reports both; analyses B
and D are 0.90-ONLY — level scope restated at analysis B below), with
the SHARED conformal-level constant
α_c = 1 − level — α_c = 0.10 at level 0.90 and α_c = 0.05 at level
0.95 — used by EVERY interval method's order statistics and tails
below (pooled, Mondrian, CQR, isotonic-D); small-bin rule — bins with test n < 30 in a repeat are
excluded from that repeat's error means, per-bin n always reported;
material margin M = 0.01.

**Exact E422 split algorithm (every repeat; the repeat's only random
draw):** first the analysis POOL is bound: accessions =
`sorted(set(ids))` over the ids of rows that are ANALYSIS-ELIGIBLE in
the cumulative label ledger as of the evaluation cut — **an
analysis-eligible row has all four channels present and finite**
(lDDT-Cα computable, d_Kabsch computable, pLDDT within the registered
grid [0, 100.01), B-factor robust-z computable); **every non-eligible
row is typed and tallied** — in registered precedence order:
`null_label` (lDDT-Cα absent), `nonfinite_label` (lDDT-Cα
non-finite), `null_feature` (pLDDT absent), `nonfinite_feature`
(pLDDT non-finite), `feature_out_of_grid` (pLDDT outside [0, 100.01)),
`null_descriptive_channel` (d_Kabsch or B-factor robust-z absent or
non-finite; checked after the feature and grid tests) — first failure
in that order names the row, and **the six-key tally is reported with
every evaluation** (committed implementation:
experiments/e422_protocol.py `_prep_rows`). Accessions admitted to the census that yielded ZERO
eligible rows do NOT enter the permutation — including them would
change n and every split assignment — and their exclusion is tallied in
the batch checkpoint. The canonical accession list is Python
`sorted()` order (PDB ids are ASCII, so codepoint order = byte
order); it is recorded in the batch checkpoint. Then
let n = len(accessions), rng = `np.random.default_rng(seed)`, perm =
`rng.permutation(n)` — no other rng call precedes perm, and entry i of
the permutation denotes accessions[perm[i]]. Slice the permuted
indices with integer arithmetic: t = n//2; c = (n − t)//2; train =
perm[:t]; calibration = perm[t:t+c]; test = perm[t+c:]. Accession-level
(all rows of an accession stay in one split), so the three splits are
mutually accession-disjoint; all methods run on the SAME three splits
within a repeat (paired). STRICT role separation: train fits each
method's center/quantile functions (pooled median; per-bin medians; CQR
quantile functions), calibration fits the conformal quantiles, test
evaluates coverage — no method touches another split's rows. "No
tuning" below means no hyperparameter search; the fitted quantities
named in §3A are the entire learning done by any method.
**Minimum-pool precondition (LAB-06):** an evaluation executes only if
n ≥ 4 accessions AND every level's calibration slice satisfies
k = ⌈(n_cal+1)(1−α_c)⌉ ≤ n_cal for every method's conformal step; a
pool violating either is typed `DEFERRED_INSUFFICIENT_SUPPORT` for the
whole evaluation (a §6 non-outcome; the next scheduled evaluation
supersedes it). **The k ≤ n_cal check is a PER-SEED, PER-LEVEL
pre-screen on CALIBRATION ROW counts, not accession counts (gate-(c)
disposition, 2026-09-09):** for each registered seed, n_cal is the
number of LEDGER ROWS landing in that seed's calibration slice
(accessions carry unequal row counts, so row counts are the binding
basis); k is computed against α_c of that level; if k > n_cal for ANY
(seed, level), the whole evaluation is `DEFERRED` — no repeat runs, no
partial results are emitted. Committed implementation:
`run_evaluation` returns status `DEFERRED` below the accession floor,
runs the per-seed calibration-row pre-screen before any repeat, and
`conformal_order_stat` raises a specification violation when k > n_cal
(the pre-screen makes this unreachable in the registered path — a
belt-and-braces assertion, not a control-flow branch).

**Effective repeats and MCSE.** R_eff(method, level) = the number of
repeats in which that method produced intervals at that level (typed
`FIT_FAILURE` drops, §3A.3). Paired statistics use R_eff(pair) =
repeats where BOTH compared methods contribute AND ≥1 bin is eligible.
Bin-level aggregates use R_bin(bin, method, level) = the number of
repeats in which that bin is ELIGIBLE for that method AND the method
produced intervals — a bin's per-bin coverage exists only in those
repeats, so small-bin eligibility varying across repeats is handled PER
BIN, never collapsed to R. Every reported aggregate carries its count
(R_eff, R_eff(pair), or R_bin). MCSE = sample sd (ddof=1) over the
statistic's contributing repeats / √count, against the standard-normal
2-sigma critical value. If the count is < 2, sd is undefined: the
affected statistic is typed `DEFERRED_INSUFFICIENT_SUPPORT` (the
e421-P §2.1 status discipline) — point estimate reported where one
exists, no 2-sigma claim; for P7a/P7b, whole-comparison deferral
applies ONLY when the comparison's own count R_eff(pair) < 2 (§7
deferral scope) — a per-bin R_bin < 2 defers only that bin's statistic.
A repeat with no eligible bins contributes no bin-level
statistics (recorded, with per-bin n).

- **A (primary).** Split-conformal intervals for lDDT-Cα at 0.90/0.95:
  1. **Pooled** (comparator): center ŷ_pool = training-pool median lDDT
     (fitted on the train split — the method's ONLY fitted quantity, a
     single constant); conformity score = |y − ŷ_pool|; q = ⌈(n_cal+1)(1−α_c)⌉
     order statistic of calibration scores; interval = ŷ_pool ± q.
  2. **Mondrian** (pLDDT strata — AF-observable): per-bin center
     c_k = training-pool median lDDT within pLDDT bin k (train split);
     calibration residuals r = y − c_k(bin(x)); per-bin
     q_k = ⌈(n_k+1)(1−α_c)⌉ order statistic of |r| over calibration rows in
     bin k; interval = c_k(bin(x)) ± q_k. A bin whose train split has
     n_k < 30 OR whose calibration split has n_k < 30 defers entirely to
     the pooled interval for that bin (recorded per bin per repeat).
  3. **CQR** (pLDDT feature only — E422-2): linear quantile regression of
     lDDT on the design [pLDDT, pLDDT²] (column order fixed) with the
     single constant column supplied by the estimator:
     `sklearn.linear_model.QuantileRegressor` with fit_intercept=True —
     registered this way so the design contains EXACTLY ONE constant
     column ("[1, pLDDT, pLDDT²] plus fit_intercept=True" would
     duplicate the intercept and is NOT the registered design);
     estimator penalty alpha = 0.0 (sklearn's L1 regularization weight —
     UNRELATED to the conformal level; 0.0 = unregularized pinball
     loss); solver='highs'; quantile=τ; all constants registered here,
     no tuning. **Level naming: α_c is the SHARED protocol constant
     above (1 − level; 0.10 at 0.90, 0.05 at 0.95) — distinct from the
     estimator's penalty alpha = 0.0.** The CQR-specific mapping, bound
     so an implementation cannot misread it: level 0.90 → τ_lo = 0.10,
     τ_hi = 0.90, score quantile 0.90; level 0.95 → τ_lo = 0.05,
     τ_hi = 0.95, score quantile 0.95. Tail assignment
     FIXED: q̂_lo ≡ fit at τ_lo = α_c (lower-bound function), q̂_hi ≡
     fit at τ_hi = 1−α_c (upper-bound function) — no implementer
     discretion, and reversing the pair is a specification violation.
     The α_c/1−α_c pair (NOT the textbook α_c/2, 1−α_c/2) is the
     INTENTIONAL registered convention: conformalization recalibrates
     marginal coverage to 1−α_c for any fitted pair — the α_c-tailed
     fits start narrower and the conformal quantile Q compensates; only
     the starting point differs, the guarantee is identical. Quantile
     crossing between the fitted curves is LEFT AS FITTED (CQR's
     conformal validity holds for arbitrary fitted functions; any
     post-hoc crossing repair would be unregistered tuning). A fit that
     fails to converge is a typed `FIT_FAILURE` for that (repeat,
     level): the repeat is dropped and R_eff reported (§7 P7b) — no
     silent substitution. **Registered detection mechanism (committed in
     experiments/e422_protocol.py `_fit_cqr_quantile`, mandated for any
     reimplementation):** every fit executes under
     `warnings.catch_warnings()`; a raised exception OR any emitted
     `sklearn.exceptions.ConvergenceWarning` is `FIT_FAILURE` —
     sklearn's QuantileRegressor signals solver failure by EMITTING the
     warning and RETURNING the fit, so an un-guarded fit would silently
     use non-converged quantiles in violation of the no-silent-drops
     fence. Scores = max(q̂_lo(x) − y, y − q̂_hi(x)) on
     calibration; Q = the ⌈(n_cal+1)(1−α_c)⌉-th order statistic of the
     scores; interval = [q̂_lo(x) − Q, q̂_hi(x) + Q].
  4. **Isotonic-pLDDT baseline** (analysis D): sklearn IsotonicRegression
     (out-of-bounds="clip", no tuning) fit on train; pooled conformal
     residual about its predictions — conformal step identical to A.1's
     construction at the shared α_c.
  **No hyperparameter tuning exists anywhere in this registration; any
  tuning or model change post-dispatch requires a new E-number.** Typed
  statuses `UNDER_SUPPORTED_STRATUM` (§3B), `FIT_FAILURE` (§3A.3), and
  `DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE` (§6 close status) are
  enumerated in the finalized analysis-status artifact BEFORE anchor (the
  e421-P §2.1 discipline, carried here).

- **B (primary).** Conditional coverage per pLDDT bin (six-bin grid above)
  at 0.90: per-bin |coverage − nominal| with per-repeat paired MCSE.
  **Level scope (T-11):** analyses B and D are registered at 0.90 ONLY
  (consistent with P6's level binding); analysis A runs at both levels
  and always reports both. **Secondary (reference-conditioned
  diagnostic — E422-1/E422-2):** the
  joint grid pLDDT bin × B-factor-z band (≤−1, (−1,+1), ≥+1) is
  evaluated DESCRIPTIVELY at 0.90 only, per repeat across R = 200 on
  that repeat's test split (the same repeats as the primaries), the
  cells conditioning on the POOLED intervals (the analysis-A.1
  comparator construction — the joint grid extends analysis B's
  conditional view, whose coverage statistic is the pooled method's):
  joint
  cells with test n < 100 (rows) or protein n < 20 IN A REPEAT are
  ineligible in that repeat; a joint cell is REPORTED iff it is
  eligible in ≥ 2 contributing repeats, with MCSE over those repeats
  (the same count ≥ 2 discipline as the per-bin tables), otherwise it
  is typed `UNDER_SUPPORTED_STRATUM` — reported with its n (per-repeat
  means over all R repeats), carrying NO
  coverage claim; no joint cell is a primary endpoint and no joint
  floor is claimed at registration. The joint grid's pLDDT axis
  INHERITS the bound conditional grid at each evaluation (original
  six-bin, or the §4 merged grid once bound — S-04).

- **C (secondary).** The e421-P §11 Phase-3 comparison (classwise AD-CP vs
  pooled-CalPro on the binary endpoint) executes here unchanged as a
  secondary analysis, cross-registered by reference to e421-P §11's frozen
  design. **Binding (F-03): the analysis-C design of record is
  docs/e421_prospective_registration.md §11 as pinned by whole-file
  content SHA-256 at the gate-(a) freeze (recorded in the freeze
  artifact); any change to e421-P §11 after that freeze does NOT enter
  e422 — executing the changed design here would require a new
  E-number.** Analysis C evaluates at e422's §6 checkpoints under the
  same schedule; the §2 analysis-C exception to the four-cell ban applies, and
  the §3 grid (half-open intervals) is AUTHORITATIVE for every analysis
  that uses it, including analysis C — e421-P §11's "(0,50)"-style
  notation is non-normative prose (the committed
  experiments/e421_conformal.py `PLDDT_BINS` is the shared ground
  truth; LAB-09). If analysis C's own e421-P §11 minimum-support gate
  has not passed by batch 12, analysis C's status at close is
  `DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE` (the §6 close discipline
  applies to analysis C as to every registered prediction; F-11).

- **D (secondary).** Isotonic-pLDDT calibration quality at 0.90: per-bin
  residual bias (mean signed residual y − method center) and width (mean
  interval width) vs the pooled baseline (descriptive; no win rule;
  MCSE over R_bin with counts reported — the §3 effective-repeats
  discipline).

## 4. Feasibility (registered arithmetic and its limits)

From the exploratory census (development-only): five-bin pLDDT shares were
0.72% / 1.76% / 2.82% / 9.74% / 84.96% (177/431/690/2,381/20,778 of
24,457; the five bands are [0,50), [50,70), [70,80), [80,90),
[90,100.01) — the §4 merge-step-1 grid below). At ~21.9k prospective residues/week, 12 weeks yield ≈ **1.9k /
4.6k / 7.4k / 25.6k / 223k** residues per five-bin band. On the registered
six-bin grid the committed exploratory ledger counts are **177 / 177 /
254 / 690 / 2,381 / 20,778 of 24,457** (shares **0.724% / 0.724% /
1.04% / 2.82% / 9.74% / 84.96%**): **the sparsest six-bin strata are
0–50 AND 50–60, TIED at 0.724% each**, each reaching 1,000 residues at
**≈6.3 weeks** (1,000 / (0.00724 × 21,900) ≈ 6.3) — all pLDDT strata
exceed 1,000 residues within the 12-week primary window.

**Limits of this arithmetic (E422-1):** these are MARGINAL (per
pLDDT bin) projections only. They are NOT evidence about the joint
pLDDT × B-factor grid — the exploratory readback already contains a
0-row/0-protein 0–50 × low-B joint cell, and several other joint cells are
far below 100 rows/20 proteins; joint cells therefore carry no floor and no
claim (§3B secondary). **Primary pLDDT-bin floors: ≥1,000 ledger residues
AND ≥20 proteins per bin**, measured at the **batch-8 (day-56)
evaluation** (floor-checkpoint relocation, T-02/LAB-08 resolution below);
if any bin misses EITHER floor there (protein floor OR residue floor),
the
registered fallback — decided BEFORE unblinding
any coverage table and recorded in the batch-8 checkpoint — is to merge
adjacent bins for the conditional analysis, with the EXACT resulting
grids registered (no interval-reading ambiguity): **step 1** merges
50–60 ∪ 60–70 → the grid {[0,50), [50,70), [70,80), [80,90),
[90,100.01)}; **step 2** (taken only if a floor is still missed after
step 1) merges [0,50) ∪ [50,70) ∪ [70,80) into the single contiguous
[0,80) → the grid {[0,80), [80,90), [90,100.01)}; or, failing that, to
re-scope the
conditional analysis to the three top bins and report the low bins
descriptively. The decision binds P6/P7 to the original or merged grid
from that point on — no post-hoc grid redefinition. Committed
implementation (gate-(c) disposition, 2026-09-09): the decision logic
is `floor_decision` in experiments/e422_evaldriver.py (merged-bin
residues fold as constituent sums; merged-bin protein counts are exact
set unions of per-bin accession sets, so a protein with rows in several
constituents counts once), and the conditional evaluation on a merged
or re-scoped grid runs the SAME registered code path with the analysis
grid rebinded via `use_grid` in experiments/e422_protocol.py — the
analysis pool, accession list, and splits are UNCHANGED by the grid
decision; only the conditional strata change (rows below a re-scoped
grid's low edge are retained in the pool but unstratified).

**Floor-checkpoint relocation (gate-(b) disposition, 2026-09-09;
T-02/DP-09/LAB-08):** the floor checkpoint was moved from the batch-4
(day-28) to the batch-8 (day-56) evaluation. On the registration's own
arithmetic, [0,50) and [50,60) reach the 1,000-residue floor only at
≈6.3 weeks, so a batch-4 floor measurement would fire the merge
fallback for bins whose floors were UNREACHABLE at measurement time —
the merge would be an artifact of checkpoint placement, not evidence.
At day 56 (≈175k cumulative residues) both low bins project ≈1.27k
residues — above the residue floor; the protein floor (≥20 proteins per
bin) remains an empirical question and the merge fallback covers its
failure. **Registered consequence:** P6, P7a, and P7b are EXPECTED to
run on the original SIX-bin grid (their first evaluable evaluation is
the batch-8 checkpoint); the merged grids above are the registered
fail-safe, not the expected path. The batch-4 (day-28) evaluation
remains scheduled: P5-90/P5-95 are evaluable there (no strata needed),
and the batch-4 checkpoint additionally records the DESCRIPTIVE
floor-trajectory table (per-bin cumulative residues/proteins to date,
non-binding). A protein- or
residue-starved low bin does NOT affect the primary
marginal-coverage analysis (P5), which needs no strata.

**Timeline: G2-grade evidence in one 12-week primary window** — versus
≈4.3 years under the four-cell design. The exploratory 24,457-label ledger
is development-only (method development and sanity checks may use it);
**no headline number in any output may be computed from it.**

## 5. Provenance (e421's gaps closed by construction, from residue 1)

v27 per-row envelopes, per-fetch cryptographic receipts, and fold/group
bindings are emitted AT LABEL TIME for every batch (the e421 corrective
actions are entry conditions here, not retrofits). Batch lifecycle as
e421-P §7: census receipt → per-batch freeze receipt → label construction →
batch checkpoint (SHA-256 manifest) → OTS anchor → ledger append; labels
never read pLDDT for cell construction — here pLDDT IS the declared
feature, so the applicable fence is the §1 channel-role split (labels
never enter feature construction; features never enter label
construction).

## 6. Lifecycle gates

(b) five-role reviews → (c) detached audit → (a) registration checkpoint
freeze (exact T0 bound here, as in e421-P) → (d) OTS anchor → (e)
bounded dispatch → (f) weekly batches with per-batch freezes → (g) 12-week
primary window close → (h) G2 evaluation + fresh independent audit.

**Ordering amendment (2026-09-09, on main's dispatch authorization,
DECISIONS same date):** the text gates (b)/(c) now precede the freeze so
the anchored object is the FINAL text. Exact T0 (second-level UTC
instant) is bound IN the frozen text at gate (a); the snapshot
presence-manifest mechanism (§2 amendment: schema, selection rule,
fail-closed temporal check, canonical-JSON SHA-256 checkpointing)
freezes at dispatch preparation and must precede batch 1 — the
manifest CONTENT then accumulates incrementally per accession at first
appearance (§2), each per-batch checkpoint digest OTS-anchored with the
batch checkpoint. Gate (e)
dispatch: authorized by main 2026-09-09 (DECISIONS same date) — dispatch
and execution decisions for THIS registration are delegated to w1:pG
under the standing governance rules (receipted fetches, typed failures,
per-batch freezes, ledger discipline, no silent drops), every action
reviewable and auditable.

**Gate-(a) FREEZE (2026-09-09, w1:pG under delegated authority;
DECISIONS same date).** This text as committed at the freeze commit is
the FROZEN artifact of record — the anchored object at gate (d); ANY
change after this commit is a new E-number. **T0 =
2026-09-11T00:00:00Z** (exact second-level UTC; batch 1 = (T0, T0+7d] =
(2026-09-11T00:00:00Z, 2026-09-18T00:00:00Z]; the batch-4 evaluation is
day-28, batch-8 day-56, batch-12 close day-84, all UTC). **Bound
objects at freeze (content SHA-256s and the full implementation-set
enumeration recorded in the freeze artifact,
docs/e422_freeze_gate_a_20260909.md):** the analysis-grid source blob
f55ac83de71ac2b0fb65a64a0a2c660e4e6842a9
(experiments/e421_conformal.py `PLDDT_BINS`, cat-file-verified);
docs/e421_label_rules.json; docs/e421_source_allowlist.json; analysis
C's design of record — docs/e421_prospective_registration.md by
whole-file content SHA-256 at the freeze commit (§3C binding); and the
committed implementation set experiments/e422_protocol.py,
e422_manifest.py, e422_status.py, e422_runner.py, e422_evaldriver.py
(each content-SHA-256-bound in the freeze artifact). **Dependency lock
verified at freeze:** python 3.12.13 / numpy 1.26.4 / scikit-learn
1.9.0 / scipy 1.17.1 (shared `.venv`). The §3B joint-grid
implementation was reviewed pre-freeze by detached review wf_75c31e8b
(2 auditors, READY_WITH_FINDINGS, 13 findings — 0 BLOCKER/MAJOR — all
disposed; record docs/audit_e422_joint_grid_review_20260909.md).

**Evaluation schedule and binding (analyses A–D; P5, P6, P7a, P7b).**
Cumulative evaluations at the batch-4 (day-28), batch-8 (day-56), and
batch-12 (day-84) checkpoints, UTC — each evaluation uses the FULL
cumulative label
ledger to that date, never a sliding window. **Checkpoint roles
(floor-checkpoint relocation, §4):** the batch-4 evaluation evaluates
P5-90/P5-95 (no strata needed) and records the DESCRIPTIVE
floor-trajectory table; **its artifact carries those two things ONLY —
all conditional-coverage outputs (per-bin tables, the §3B joint
grid, P6/P7 predictions) are suppressed at batch-4 until the
batch-8 checkpoint, so no
conditional table exists before the count-based grid decision**
(registered suppression fence; committed implementation
`_scope_batch4` in experiments/e422_evaldriver.py); the §4 floor
measurement and merge decision
run at the batch-8 evaluation (their binding point), which is also
P6/P7a/P7b's first evaluable evaluation; batch-12 is the close. A comparison
in `DEFERRED_INSUFFICIENT_SUPPORT` is a NON-outcome: it does not bind,
and the next scheduled evaluation supersedes it. The first evaluation
at which a prediction is evaluable and its decision rule yields an
outcome — P5-90 / P5-95 hold/fail (each level independently; no
combined P5 status exists); P6 replicate/falsify; P7a/P7b win,
mirrored win, or NO_DECISION — BINDS that prediction. Each prediction binds
independently at its own first evaluable checkpoint; bound outcomes
are never revised by later evaluations; every evaluation (bound or
non-outcome) is reported. Gate (g), the 12-week primary window close
(batch 12), is the final scheduled PRIMARY checkpoint: any prediction
still unbound at close is recorded
`DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE` — a registered possible
outcome, reported as such.

**Compute budget and dependency lock (AGENTS.md §4; versions observed in
the shared `.venv` on 2026-09-09 and pinned at gate (a)):** CPU-only
(AGENTS §5 statistics lane; no GPU): a full registered evaluation
(R = 200 repeats × {pooled, Mondrian, CQR, isotonic} × 2 levels on
≤ ~260k rows) fits in ≤ 4 cores and ≤ 8 GB RAM, with per-evaluation
wall-clock budgets as stated below (≤ 24 h at batches 4 and 8; ≤ 72 h
at the batch-12 close).
**Registered execution detail (LAB-05/F-02 resolution): CQR fits MAY run
in parallel across independent worker PROCESSES bounded by the ≤4-core
fence, and the wall-clock bound is per-evaluation as stated below.** A
full evaluation fits 200 repeats × 2 levels × 2 quantile tails =
**800** QuantileRegressor fits; the fits are independent across
(seed, level, tail) — each consumes only its own repeat's train split
and no randomness enters the estimator (sklearn QuantileRegressor with
solver='highs' is deterministic for a fixed design) — so a
parallel-by-(seed, level, tail) execution computes byte-identical
per-fit results to the serial run; results merge on the
(seed, level, tail) key and the R_eff/R_bin accounting is unchanged.
Measured basis in the pinned
environment (alpha=0.0, solver='highs', design [x, x²], single core):
0.84 s @ 5k rows, 4.6 s @ 12k, 12.9 s @ 20k, 82.5 s @ 50k, 128.9 s @
60k, 570 s @ 130k (superlinear, ≈n²). **Per-evaluation wall-clock
budgets (registered):** batch-4 evaluation ≈ 800 × 66 s ≈ 15
CPU-hours — within ≤ 24 h even serially; batch-8 evaluation (train
≈ 88k rows) ≈ 800 × 254 s ≈ 56 CPU-hours — ≤ 24 h ONLY with the
registered 4-process parallelism (≈ 14 h); batch-12 close evaluation
(train ≈ 131k rows) ≈ 800 × 570 s ≈ 127 CPU-hours — **budgeted ≤ 72
wall-clock hours on the same ≤ 4-process lane (≈ 32 h expected)**. Any
evaluation exceeding its budget is reported as an execution deviation
in the evaluation artifact, never silently truncated (no repeat is
dropped for budget reasons). Locked: python 3.12.13; numpy 1.26.4
(`np.random.default_rng`, PCG64);
scikit-learn 1.9.0 (`QuantileRegressor`, `IsotonicRegression`); scipy
1.17.1 (`highs` linprog backend). Any version change post-freeze = new
E-number. (statsmodels is absent from the environment; the CQR
specification therefore pins sklearn's `QuantileRegressor` — a DIFFERENT
implementation of the same estimator family, linear quantile
regression, unregularized; implementation identity with statsmodels
QuantReg is NOT claimed — this registration pins the sklearn
implementation, installed and verified.)

## 7. Falsifiable predictions (registered before any prospective label exists)

- **P5-90 / P5-95 (independent predictions, one per level):** pooled
  split-conformal lDDT intervals achieve marginal coverage ≥ nominal −
  3 MCSE (one-sided, protein-disjoint) at level 0.90 (P5-90) and at
  level 0.95 (P5-95). Each level binds INDEPENDENTLY under the §6
  schedule; a `DEFERRED_INSUFFICIENT_SUPPORT` at one level defers only
  that level's status and never the other's — there is NO combined P5
  status; both levels are always reported together.
- **P6 (level 0.90 — bound):** pooled conditional coverage at level
  0.90 deviates from nominal by >5 MCSE in ≥2 pLDDT bins (the
  exploratory mid-bin failure replicates on fully-provenanced labels).
  Null replication is a registered possible outcome and would be
  reported as such. **Level binding:** P6 is a 0.90-only prediction —
  §3B registers the conditional analysis at 0.90; no P6-95 exists (the
  0.95 level is represented by P5-95's marginal prediction).
  Eligible bins and the bound grid are exactly P7a/P7b's; each bin's
  MCSE uses that bin's own contributing-repeat count
  R_bin(bin, pooled, 0.90) — pooled incurs no `FIT_FAILURE`, but
  small-bin eligibility can vary by repeat, so R_bin ≤ R per bin; every
  per-bin table reports R_bin beside n. **Counting rule:** only bins
  with R_bin ≥ 2 (defined MCSE) count toward the ≥2-bin threshold; a
  bin with R_bin < 2 is typed `DEFERRED_INSUFFICIENT_SUPPORT` for its
  per-bin statistic and is NOT counted as P6 evidence in either
  direction (it neither advances the ≥2 count nor counts as a
  null-replication bin). If fewer than 2 bins carry R_bin ≥ 2 at an
  evaluation, P6 as a whole is `DEFERRED_INSUFFICIENT_SUPPORT` at that
  evaluation (a §6 non-outcome).
- **P7a (co-primary, two-sided, material margin M = 0.01) — Mondrian vs
  pooled:** for method m and repeat r, maxdev_r(m) = max over eligible
  bins of |coverage_bin,r(m) − 0.90|. Eligible bins in repeat r are
  bins passing the small-bin rule (test n ≥ 30) on the grid AS BOUND at
  the batch-8 checkpoint (original six-bin grid, or the §4 registered
  merge — bound pre-unblinding, never redefined after); eligibility
  depends only on the test split, so the eligible set is COMMON to all
  methods within a repeat by construction. d_r = maxdev_r(Mondrian) −
  maxdev_r(pooled), paired within repeat, over R_eff(pair) contributing
  repeats (§3 effective-repeats rule; R_eff reported). Mondrian wins
  iff mean(d_r) < −M and |mean(d_r)| > 2·sd(d_r)/√R_eff (sd ddof=1);
  mirrored for a pooled win; otherwise NO_DECISION; R_eff < 2 →
  `DEFERRED_INSUFFICIENT_SUPPORT` (§3).
- **P7b (co-primary, two-sided, M = 0.01) — CQR vs pooled:** the
  IDENTICAL maxdev_r/d_r rule and outcomes with CQR in place of
  Mondrian; repeats drop additionally via typed `FIT_FAILURE` (§3A.3).
- **Multiplicity policy (registered pre-data; mirrors e421-P §8):**
  P7a and P7b are SEPARATE, STANDALONE registered comparisons — each is
  evaluated on its own against the same margin; no hierarchy and no
  single-winner selection between them, and no post-hoc choice of which
  comparison to report. Both results are always reported in full (all
  per-bin tables, whichever method wins or loses either comparison).
  **Deferral scope (executable):** P7a/P7b defer or bind as WHOLE
  comparisons on R_eff(pair) alone; a bin with R_bin < 2 defers ONLY
  that bin's per-bin statistic (P6 and the per-bin tables) and NEVER
  defers P7a/P7b — maxdev_r is a per-repeat max over that repeat's
  eligible bins and consumes no per-bin MCSE. First evaluable
  evaluation (§6 evaluation schedule and binding) binds each
  comparison independently.

## 8. Relationship to e421-P

Independent registrations, independently decided: shared infrastructure,
separate E-numbers, ledgers, and checkpoint chains. If both are dispatched,
e422's binary-endpoint secondary (analysis C) and e421-P's primary census
are computed from disjoint batch streams to avoid double-use of labels in
registered claims. **"Disjoint batch streams" is defined (DP-05):** each
dispatched registration fetches, receipts, freezes, and ledgers ITS OWN
copies of the weekly census over the allowlisted lanes — no label row,
receipt, or checkpoint artifact is shared between the two streams; a
given residue's label therefore enters exactly ONE registration's
registered claims. The cost is duplicate fetching (accepted); the
guarantee is that no registered claim in either object consumes a label
that entered another object's claim chain; each stream receives the
FULL released-entry census, so e422's §4 yield arithmetic is unaffected
by whether e421-P is also dispatched (F-09). Dispatch priority was
decided by main on 2026-09-09 (DECISIONS same date: e422 dispatches
first; e421-P remains independently decided) — any change to that
priority is main's alone (F-08).

## 9. Fences (all carried forward, unchanged)

No raw redistribution; public endpoints only; no credentials; typed
failures only; no silent imputation; no mutable "latest" references;
labels never enter feature construction and features never enter label
construction; reference-derived quantities are never deployable features
(§1 channel roles); retention 90-day event + 2027 backstop; negative
results logged, never buried; amendments via dated DECISIONS entries only
(pre-dispatch); post-dispatch changes require a new E-number.
