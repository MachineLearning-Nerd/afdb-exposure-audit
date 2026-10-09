"""e422 frozen analysis protocol: executable form of docs/e422_registration.md.

Implements EXACTLY the registered text (sections 3, 3A, 3B, 6, 7). Any
deviation from the registration is a specification violation; the
registration is the artifact of record and this module must never be
changed post-freeze except by a new E-number.

Estimand: distribution-free split-conformal intervals for continuous
lDDT-Ca from AF-observable inputs (pLDDT sole feature; reference-derived
quantities are labels/stratification only). Methods: pooled (comparator),
Mondrian (pLDDT strata), CQR (sklearn QuantileRegressor), isotonic-pLDDT
baseline. R = 200 repeats, seeds 0-199, np.random.default_rng; the
permutation draw in split_indices is each repeat's ONLY random draw.

Offline: operates on an in-memory cumulative label ledger. No network,
no file writes, deterministic given the ledger.

Typed statuses (registered vocabulary, docs/e422_registration.md):
  FIT_FAILURE                  - CQR fit failed to converge (repeat dropped
                                 for that method/level; R_eff reported)
  DEFERRED_INSUFFICIENT_SUPPORT- statistic/comparison lacks support
                                 (count < 2) at this evaluation; a
                                 NON-outcome under the section 6 schedule
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field

import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import QuantileRegressor

# --- registered constants (docs/e422_registration.md section 3) ------------

PLDDT_BINS = ((0.0, 50.0), (50.0, 60.0), (60.0, 70.0), (70.0, 80.0),
              (80.0, 90.0), (90.0, 100.01))
LEVELS = (0.90, 0.95)
SMALL_BIN_MIN_TEST_N = 30      # small-bin rule: test n >= 30 per bin/repeat
MONDRIAN_MIN_N = 30            # per-bin train AND calibration n >= 30
MARGIN_M = 0.01                # material margin for P7a/P7b
P5_MCSE_K = 3.0                # P5 one-sided threshold: nominal - 3 MCSE
P6_MCSE_K = 5.0                # P6 deviation threshold: >5 MCSE
CRIT_2SIGMA = 2.0              # standard-normal 2-sigma critical value
DEFERRED = "DEFERRED_INSUFFICIENT_SUPPORT"
FIT_FAILURE = "FIT_FAILURE"

# Section 3B secondary (reference-conditioned diagnostic): the joint
# grid pLDDT bin x B-factor-z band, evaluated DESCRIPTIVELY at 0.90
# only, on the same repeats as the primaries.  Bands are the registered
# <=-1 / (-1,+1) / >=+1 bands; a joint cell is INELIGIBLE in a repeat with
# test n < 100 rows or protein n < 20; a cell is REPORTED iff eligible
# in >= 2 contributing repeats (MCSE over those repeats), else typed
# UNDER_SUPPORTED_STRATUM (reported with its n, NO coverage claim).
# The joint grid's pLDDT axis INHERITS the bound conditional grid.
BZ_BANDS = ((-math.inf, -1.0), (-1.0, 1.0), (1.0, math.inf))
JOINT_MIN_TEST_N = 100
JOINT_MIN_PROTEINS = 20
JOINT_MIN_REPEATS = 2
UNDER_SUPPORTED_STRATUM = "UNDER_SUPPORTED_STRATUM"

# Registration text of record cited in every evaluation artifact.  The
# value is updated ONLY at a lifecycle gate; gate (a) (2026-09-09) set the
# FINAL frozen bind.  Both identifiers are labeled: git-blob (cat-file
# resolvable) and content SHA-256 of the identical bytes (finding 3B-4).
# ANY change to the registration text after this bind is a new E-number.
REGISTRATION_CITATION = (
    "docs/e422_registration.md (git-blob "
    "6e1b9ce128ac3a1be4678463bdc777c25b50ca8d, content-sha256 "
    "dc7a362773fbf3a03e759aa7d178125963a14d88d9bbd61b66389a04ba8d7b3b, "
    "commit 9dd2380b; FROZEN gate-(a) object of record, "
    "T0=2026-09-11T00:00:00Z)")

# Statuses enumerated before anchor (registration section 3A tail): the
# analysis-side statuses are FIT_FAILURE and DEFERRED_INSUFFICIENT_SUPPORT;
# DEFERRED_INSUFFICIENT_SUPPORT_AT_CLOSE is applied by the section 6
# schedule driver at the batch-12 close, never inside an evaluation.


class E422SpecViolation(RuntimeError):
    """Raised when inputs fall outside the registered specification."""


class _FitFailure(Exception):
    """Internal: CQR fit did not converge -> typed FIT_FAILURE."""


# --- split (section 3, exact algorithm) ------------------------------------


def canonical_accessions(rows):
    """sorted(set(admitted_accession_ids)) - the registered pool order."""
    return sorted({r["accession"] for r in rows})


def split_indices(n, seed):
    """Registered 50/25/25 accession-level split; the repeat's only draw.

    perm = np.random.default_rng(seed).permutation(n) is the ONLY random
    draw; entry i of perm indexes accessions[perm[i]]. Integer slices:
    t = n//2, c = (n-t)//2; train = perm[:t], cal = perm[t:t+c],
    test = perm[t+c:]. Mutually disjoint by construction.
    """
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    t = n // 2
    c = (n - t) // 2
    return perm[:t], perm[t:t + c], perm[t + c:]


def alpha_c_for(level):
    """SHARED conformal-level constant: alpha_c = 1 - level.

    0.10 at level 0.90, 0.05 at level 0.95. Used by EVERY interval
    method's order statistics and tails. DISTINCT from the CQR estimator
    penalty alpha = 0.0.
    """
    return 1.0 - level


def conformal_order_stat(scores, alpha_c):
    """ceil((n+1)(1-alpha_c))-th order statistic (1-indexed).

    Registered form for the pooled q, Mondrian q_k, and CQR Q. Under the
    registered rules (pooled: full calibration set; Mondrian: only bins
    with calibration n_k >= 30; CQR: full calibration set) k never
    exceeds n; a violation is a specification failure, not an inf.
    """
    s = np.sort(np.asarray(scores, dtype=float))
    n = s.size
    if n == 0:
        raise E422SpecViolation("empty_calibration_scores")
    k = math.ceil((n + 1) * (1.0 - alpha_c))
    if k > n:
        raise E422SpecViolation(f"order_statistic_beyond_n:k={k},n={n}")
    return float(s[k - 1])


def bin_index(plddt):
    """Index into the registered six-bin grid (half-open [lo, hi))."""
    for i, (lo, hi) in enumerate(PLDDT_BINS):
        if lo <= plddt < hi:
            return i
    raise E422SpecViolation(f"plddt_outside_registered_grid:{plddt!r}")


_BIN_UPPERS = np.array([hi for _, hi in PLDDT_BINS])
_BIN_LOWS = np.array([lo for lo, _ in PLDDT_BINS])

# --- active analysis grid (section 4 merge fail-safe) -----------------------
# The REGISTERED grid (PLDDT_BINS) stays the eligibility authority
# (bin_index / _prep_rows validate against it unconditionally).  The
# ACTIVE analysis grid governs conditional strata only (_to_bins,
# _bin_stats, the Mondrian per-bin loop); it is the registered grid
# unless use_grid() rebinds it under a recorded section 4 floor decision.
_ANALYSIS_GRID = PLDDT_BINS


def _refresh_grid_state() -> None:
    global _BIN_UPPERS, _BIN_LOWS
    _BIN_UPPERS = np.array([hi for _, hi in _ANALYSIS_GRID])
    _BIN_LOWS = np.array([lo for lo, _ in _ANALYSIS_GRID])


def use_grid(bins):
    """Rebind the ACTIVE analysis grid (section 4 merge fail-safe).

    Registered callers only: the section 6 schedule driver, ONLY on a
    recorded floor decision.  Validation: sorted, half-open, strictly
    contiguous, top edge exactly 100.01, values finite; the low edge MAY
    exceed 0 (top-3 re-scope) — rows below it are then unstratified
    (index -1), never excluded from the pool.  Returns a restore thunk.
    """
    global _ANALYSIS_GRID
    grid = tuple((float(lo), float(hi)) for lo, hi in bins)
    if not grid:
        raise E422SpecViolation("empty_grid")
    if grid[-1][1] != 100.01:
        raise E422SpecViolation(f"grid_top_edge:{grid[-1][1]!r}")
    for (lo, hi) in grid:
        if not (math.isfinite(lo) and math.isfinite(hi) and lo < hi):
            raise E422SpecViolation(f"grid_bin:{lo!r},{hi!r}")
    for (lo_a, hi_a), (lo_b, _) in zip(grid, grid[1:]):
        if hi_a != lo_b:
            raise E422SpecViolation(f"grid_not_contiguous:{hi_a!r}->{lo_b!r}")
    prior = _ANALYSIS_GRID
    _ANALYSIS_GRID = grid
    _refresh_grid_state()
    def restore():
        global _ANALYSIS_GRID
        _ANALYSIS_GRID = prior
        _refresh_grid_state()
    return restore


def _to_bins(p):
    """Vectorized bin membership on the ACTIVE analysis grid.

    Default (registered six-bin) behavior is unchanged and strict:
    out-of-grid values raise.  On a RESTRICTED conditional grid
    (low edge > 0, the section 4 top-3 re-scope fail-safe) rows below
    the low edge map to index -1 — unstratified but RETAINED for method
    fits, so the accession pool and splits are unchanged (the merged /
    restricted grid changes only the conditional strata, never the
    registered eligibility, pool, or permutation).
    """
    p = np.asarray(p, dtype=float)
    idx = np.searchsorted(_BIN_UPPERS, p, side="right")
    if _BIN_LOWS[0] > 0.0:
        low_mask = p < _BIN_LOWS[0]
        idx = np.where(low_mask, -1, idx)
        p_eff = np.where(low_mask, _BIN_LOWS[0], p)
    else:
        p_eff = p
    if idx.size and (idx.max() >= len(_ANALYSIS_GRID) or
                     (p_eff < _BIN_LOWS[idx.clip(min=0, max=len(_ANALYSIS_GRID) - 1)]).any()):
        raise E422SpecViolation("plddt_outside_registered_grid")
    return idx


# --- method fits (section 3A) -----------------------------------------------


@dataclass
class MethodRun:
    """One method's intervals + per-bin statistics for one repeat/level."""
    method: str
    ok: bool
    fit_failure: bool = False
    detail: str = ""
    marginal_cov: float = float("nan")
    n_test: int = 0
    bins: dict = field(default_factory=dict)   # bin_idx -> {n_test, cov, width, bias}
    deferred_bins: tuple = ()                  # mondrian only: bins deferred to pooled


def _bin_stats(p_test, y_test, lo, hi):
    """Per-bin coverage/width/bias over ALL active-grid test bins (n always
    reported; unstratified rows on a restricted grid carry no entry)."""
    te_b = _to_bins(p_test)
    out = {}
    for bi in range(len(_ANALYSIS_GRID)):
        m = te_b == bi
        n = int(m.sum())
        entry = {"n_test": n}
        if n > 0:
            covered = (lo[m] <= y_test[m]) & (y_test[m] <= hi[m])
            entry["cov"] = float(covered.mean())
            entry["width"] = float((hi[m] - lo[m]).mean())
            # registered section 3D sign: mean signed residual y - center
            entry["bias"] = float((y_test[m] - (lo[m] + hi[m]) / 2.0).mean())
        out[bi] = entry
    return out


def _marginal(plddt_test, y_test, lo, hi):
    covered = (lo <= y_test) & (y_test <= hi)
    return float(covered.mean())


def _bz_bands_vec(z):
    """Vectorized registered B-factor-z bands (section 3B) — the SINGLE
    source of truth for the band edges; _bz_band derives from it."""
    z = np.asarray(z, dtype=float)
    return np.where(z <= -1.0, 0, np.where(z < 1.0, 1, 2))


def _bz_band(z):
    """Index into the registered B-factor-z bands (section 3B):
    <=-1 -> 0, (-1,+1) -> 1, >=+1 -> 2."""
    return int(_bz_bands_vec(np.array([float(z)]))[0])


def _joint_cells(p_test, y_test, z_test, acc_test, lo, hi):
    """Section 3B joint-grid cell statistics for ONE repeat (pooled
    intervals).  Cell key "bi:band" on the ACTIVE analysis grid (the
    pLDDT axis inherits the bound conditional grid).  Every cell reports
    its per-repeat n_rows / n_prot / cov; eligibility (test n >= 100 AND
    proteins >= 20) is applied at aggregation, not here."""
    te_b = _to_bins(p_test)
    z_band = _bz_bands_vec(z_test)
    cells = {}
    for bi in range(len(_ANALYSIS_GRID)):
        for band in range(len(BZ_BANDS)):
            m = (te_b == bi) & (z_band == band)
            n = int(m.sum())
            key = f"{bi}:{band}"
            if n == 0:
                cells[key] = {"n_rows": 0, "n_prot": 0, "cov": None}
                continue
            covered = (lo[m] <= y_test[m]) & (y_test[m] <= hi[m])
            cells[key] = {"n_rows": n,
                          "n_prot": len({acc_test[i] for i in np.flatnonzero(m)}),
                          "cov": float(covered.mean())}
    return cells


def aggregate_joint(repeats):
    """Section 3B aggregation over the SAME repeats as the primaries.

    A cell is eligible in a repeat iff n_rows >= 100 AND n_prot >= 20
    there; it is REPORTED iff eligible in >= 2 contributing repeats
    (MCSE over those repeats); otherwise UNDER_SUPPORTED_STRATUM,
    reported with its n (per-repeat means over ALL repeats), carrying NO
    coverage claim.
    """
    if not repeats:
        raise E422SpecViolation("no_repeats_for_joint_aggregation")
    keys = sorted(repeats[0]["joint_cells"], key=lambda k: tuple(
        int(part) for part in k.split(":")))
    out = {}
    for key in keys:
        rows_all = [rep["joint_cells"][key]["n_rows"] for rep in repeats]
        prot_all = [rep["joint_cells"][key]["n_prot"] for rep in repeats]
        eligible = [rep for rep in repeats
                    if rep["joint_cells"][key]["n_rows"] >= JOINT_MIN_TEST_N
                    and rep["joint_cells"][key]["n_prot"] >= JOINT_MIN_PROTEINS]
        entry = {"n_rows_mean_all_repeats": float(np.mean(rows_all)),
                 "n_prot_mean_all_repeats": float(np.mean(prot_all)),
                 "r_joint": len(eligible)}
        if len(eligible) >= JOINT_MIN_REPEATS:
            covs = [rep["joint_cells"][key]["cov"] for rep in eligible]
            mcse, n = _mcse(covs)
            entry["cov_mean"] = float(np.mean(covs))
            entry["mcse"] = mcse
            if mcse is None:
                entry["status"] = DEFERRED
        else:
            entry["status"] = UNDER_SUPPORTED_STRATUM
        out[key] = entry
    return out


def run_pooled(train_y, cal_y, test_y, alpha_c):
    """A.1 pooled comparator: center = train median; |y - center| scores."""
    center = float(np.median(train_y))
    q = conformal_order_stat(np.abs(cal_y - center), alpha_c)
    lo = np.full(test_y.shape, center - q)
    hi = np.full(test_y.shape, center + q)
    return center, q, lo, hi


def run_mondrian(train_p, train_y, cal_p, cal_y, test_p, test_y, alpha_c):
    """A.2 Mondrian: per-bin train-median center + per-bin conformal q_k.

    A bin whose train n_k < 30 OR whose calibration n_k < 30 defers
    ENTIRELY to the pooled interval for that bin (recorded). Takes RAW
    pLDDT arrays and bins them internally on the registered grid.
    """
    pooled_center, pooled_q, pooled_lo, pooled_hi = run_pooled(
        train_y, cal_y, test_y, alpha_c)
    lo = pooled_lo.copy()
    hi = pooled_hi.copy()
    tr_b = _to_bins(train_p)
    ca_b = _to_bins(cal_p)
    te_b = _to_bins(test_p)
    deferred = []
    for bi in range(len(_ANALYSIS_GRID)):
        tr = tr_b == bi
        ca = ca_b == bi
        if tr.sum() < MONDRIAN_MIN_N or ca.sum() < MONDRIAN_MIN_N:
            deferred.append(bi)
            continue
        c_k = float(np.median(train_y[tr]))
        r = np.abs(cal_y[ca] - c_k)
        q_k = conformal_order_stat(r, alpha_c)
        te = te_b == bi
        lo[te] = c_k - q_k
        hi[te] = c_k + q_k
    return pooled_center, pooled_q, tuple(deferred), lo, hi


def _fit_cqr_quantile(design, y, tau):
    """sklearn QuantileRegressor, all constants registered; FIT_FAILURE on
    non-convergence (exception or ConvergenceWarning). No silent retry."""
    model = QuantileRegressor(alpha=0.0, solver="highs", fit_intercept=True,
                              quantile=float(tau))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            model.fit(design, y)
        except Exception as exc:  # noqa: BLE001 - typed, never silent
            raise _FitFailure(f"quantile_fit_exception:{type(exc).__name__}") from exc
    for w in caught:
        if issubclass(w.category, ConvergenceWarning):
            raise _FitFailure("quantile_fit_convergence_warning")
    return model


def run_cqr(train_p, train_y, cal_p, cal_y, test_p, test_y, alpha_c):
    """A.3 CQR: design [pLDDT, pLDDT^2], intercept via fit_intercept=True.

    Tails FIXED: q_lo at tau_lo = alpha_c (lower-bound function), q_hi at
    tau_hi = 1 - alpha_c (upper). Reversal is a specification violation.
    Crossing left as-fitted. Scores max(q_lo - y, y - q_hi) on cal;
    Q = ceil((n_cal+1)(1-alpha_c)) order statistic; interval
    [q_lo(x) - Q, q_hi(x) + Q].
    """
    d_train = np.column_stack([train_p, train_p ** 2])
    d_cal = np.column_stack([cal_p, cal_p ** 2])
    d_test = np.column_stack([test_p, test_p ** 2])
    qlo_m = _fit_cqr_quantile(d_train, train_y, alpha_c)        # tau_lo = alpha_c
    qhi_m = _fit_cqr_quantile(d_train, train_y, 1.0 - alpha_c)  # tau_hi = 1-alpha_c
    scores = np.maximum(qlo_m.predict(d_cal) - cal_y, cal_y - qhi_m.predict(d_cal))
    q = conformal_order_stat(scores, alpha_c)
    lo = qlo_m.predict(d_test) - q
    hi = qhi_m.predict(d_test) + q
    return lo, hi


def run_isotonic(train_p, train_y, cal_p, cal_y, test_p, test_y, alpha_c):
    """D isotonic-pLDDT baseline: pooled-style conformal step (A.1
    construction at the shared alpha_c) about isotonic predictions."""
    iso = IsotonicRegression(out_of_bounds="clip")
    iso.fit(train_p, train_y)
    scores = np.abs(cal_y - iso.predict(cal_p))
    q = conformal_order_stat(scores, alpha_c)
    pred_test = iso.predict(test_p)
    return pred_test - q, pred_test + q


# --- one repeat (one seed, one level) ---------------------------------------


def _prep_rows(rows):
    """Registered channel roles: accession + plddt (feature) + lddt (label).

    Analysis eligibility (registration section 3): a row is eligible iff
    ALL FOUR channels are present and finite — lddt (label), plddt
    (feature, within the registered grid), d_kabsch and b_factor_z
    (descriptive channels).  Every non-eligible row is typed and tallied
    (never silently dropped), with precedence:

      null_label               lddt absent
      nonfinite_label          lddt present but non-finite (NaN/inf)
      null_feature             plddt absent
      nonfinite_feature        plddt present but non-finite
      feature_out_of_grid      plddt outside the registered grid
      null_descriptive_channel d_kabsch or b_factor_z absent/non-finite
                                   (checked AFTER the analysis channels and
                                    grid membership: feature issues precede
                                    descriptive-channel issues)

    Returns (accessions, per-accession (plddt_list, lddt_list, bz_list))
    plus the exclusion tally.  The B-factor robust-z list is carried for
    the section 3B joint grid (descriptive channel; never a feature).
    """
    excluded = {"null_label": 0, "nonfinite_label": 0, "null_feature": 0,
                "nonfinite_feature": 0, "feature_out_of_grid": 0,
                "null_descriptive_channel": 0}
    by_acc = {}
    for r in rows:
        y = r.get("lddt")
        p = r.get("plddt")
        if y is None:
            excluded["null_label"] += 1
            continue
        if not math.isfinite(float(y)):
            excluded["nonfinite_label"] += 1
            continue
        if p is None:
            excluded["null_feature"] += 1
            continue
        if not math.isfinite(float(p)):
            excluded["nonfinite_feature"] += 1
            continue
        try:
            bin_index(float(p))
        except E422SpecViolation:
            excluded["feature_out_of_grid"] += 1
            continue
        d = r.get("d_kabsch")
        z = r.get("b_factor_z")
        if (d is None or z is None
                or not math.isfinite(float(d)) or not math.isfinite(float(z))):
            excluded["null_descriptive_channel"] += 1
            continue
        by_acc.setdefault(r["accession"], ([], [], []))
        by_acc[r["accession"]][0].append(float(p))
        by_acc[r["accession"]][1].append(float(y))
        by_acc[r["accession"]][2].append(float(z))
    accessions = sorted(by_acc)
    return accessions, by_acc, excluded


def run_repeat(rows, seed, level):
    """One repeat at one level: split once, run all methods on the SAME
    split (paired). Returns per-method MethodRun entries."""
    accessions, by_acc, _ = _prep_rows(rows)
    n = len(accessions)
    if n < 4:
        raise E422SpecViolation(
            f"pool_too_small_for_registered_split:n_accessions={n}")
    tr_idx, ca_idx, te_idx = split_indices(n, seed)
    alpha_c = alpha_c_for(level)

    tr_p, tr_y, ca_p, ca_y, te_p, te_y = [], [], [], [], [], []
    te_z, te_acc = [], []
    for i in tr_idx:
        ps, ys, _zs = by_acc[accessions[i]]
        tr_p.extend(ps)
        tr_y.extend(ys)
    for i in ca_idx:
        ps, ys, _zs = by_acc[accessions[i]]
        ca_p.extend(ps)
        ca_y.extend(ys)
    for i in te_idx:
        ps, ys, zs = by_acc[accessions[i]]
        te_p.extend(ps)
        te_y.extend(ys)
        te_z.extend(zs)
        te_acc.extend([accessions[i]] * len(ps))
    tr_p = np.array(tr_p)
    tr_y = np.array(tr_y)
    ca_p = np.array(ca_p)
    ca_y = np.array(ca_y)
    te_p = np.array(te_p)
    te_y = np.array(te_y)
    te_z = np.array(te_z)

    runs = {}
    # pooled - incurs no FIT_FAILURE (single constant center); the
    # section 3B joint grid rides the SAME split on the pooled intervals
    # (computed at the registered 0.90 level only)
    c, q, lo, hi = run_pooled(tr_y, ca_y, te_y, alpha_c)
    runs["pooled"] = MethodRun("pooled", True, marginal_cov=_marginal(te_p, te_y, lo, hi),
                               n_test=int(te_y.size), bins=_bin_stats(te_p, te_y, lo, hi))
    joint_cells = (_joint_cells(te_p, te_y, te_z, te_acc, lo, hi)
                   if level == 0.90 else {})
    # mondrian
    pc, pq, deferred, lo, hi = run_mondrian(tr_p, tr_y, ca_p, ca_y, te_p, te_y,
                                            alpha_c)
    runs["mondrian"] = MethodRun("mondrian", True, marginal_cov=_marginal(te_p, te_y, lo, hi),
                                 n_test=int(te_y.size), bins=_bin_stats(te_p, te_y, lo, hi),
                                 deferred_bins=deferred)
    # cqr - typed FIT_FAILURE drops this (method, repeat, level)
    try:
        lo, hi = run_cqr(tr_p, tr_y, ca_p, ca_y, te_p, te_y, alpha_c)
    except _FitFailure as exc:
        runs["cqr"] = MethodRun("cqr", False, fit_failure=True, detail=str(exc))
    else:
        runs["cqr"] = MethodRun("cqr", True, marginal_cov=_marginal(te_p, te_y, lo, hi),
                                n_test=int(te_y.size), bins=_bin_stats(te_p, te_y, lo, hi))
    # isotonic - deterministic; registered construction cannot fail to fit
    lo, hi = run_isotonic(tr_p, tr_y, ca_p, ca_y, te_p, te_y, alpha_c)
    runs["isotonic"] = MethodRun("isotonic", True, marginal_cov=_marginal(te_p, te_y, lo, hi),
                                 n_test=int(te_y.size), bins=_bin_stats(te_p, te_y, lo, hi))
    return {"seed": int(seed), "level": level, "n_accessions": n,
            "split_sizes": {"train": int(tr_idx.size), "cal": int(ca_idx.size),
                            "test": int(te_idx.size)},
            "methods": runs, "joint_cells": joint_cells}


# --- aggregation (section 3 effective repeats + section 7 rules) ------------


def _eligible_bins(bin_table):
    """Bins passing the small-bin rule (test n >= 30) in one repeat."""
    return {bi for bi, st in bin_table.items() if st["n_test"] >= SMALL_BIN_MIN_TEST_N}


def _mcse(vals):
    vals = np.asarray(vals, dtype=float)
    n = vals.size
    if n < 2:
        return None, n
    return float(vals.std(ddof=1) / math.sqrt(n)), n


def aggregate_level(repeats, level):
    """Aggregate one level's repeats into the registered statistics.

    Effective counts exactly as registered: R_eff(method) = repeats where
    the method produced intervals; R_bin(bin, method) = repeats where the
    bin is eligible AND the method produced intervals; R_eff(pair) =
    repeats where both compared methods contribute AND >=1 bin eligible.
    A count < 2 types the affected statistic DEFERRED_INSUFFICIENT_SUPPORT.
    """
    alpha_c = alpha_c_for(level)
    out = {"level": level, "alpha_c": alpha_c, "n_repeats_requested": len(repeats)}
    r_eff = {}
    r_bin = {m: {} for m in ("pooled", "mondrian", "cqr", "isotonic")}
    bin_cov = {m: {} for m in ("pooled", "mondrian", "cqr", "isotonic")}
    marginal = {m: [] for m in ("pooled", "mondrian", "cqr", "isotonic")}
    deferred_tally = {}
    for rep in repeats:
        for m, run in rep["methods"].items():
            if not run.ok:
                continue
            r_eff.setdefault(m, 0)
            r_eff[m] += 1
            marginal[m].append(run.marginal_cov)
            elig = _eligible_bins(run.bins)
            for bi in elig:
                r_bin[m].setdefault(bi, 0)
                r_bin[m][bi] += 1
                bin_cov[m].setdefault(bi, []).append(run.bins[bi]["cov"])
            if m == "mondrian":
                for bi in run.deferred_bins:
                    deferred_tally[bi] = deferred_tally.get(bi, 0) + 1
    out["r_eff"] = r_eff
    out["r_bin"] = {m: {str(bi): c for bi, c in d.items()} for m, d in r_bin.items()}
    out["marginal"] = {}
    for m, vals in marginal.items():
        mcse, n = _mcse(vals)
        entry = {"r_eff": n, "mean": float(np.mean(vals)) if n else None,
                 "mcse": mcse}
        if mcse is None and n is not None:
            entry["status"] = DEFERRED
        out["marginal"][m] = entry
    out["per_bin"] = {}
    for m in bin_cov:
        out["per_bin"][m] = {}
        for bi, vals in sorted(bin_cov[m].items()):
            mcse, n = _mcse(vals)
            entry = {"r_bin": n,
                     "cov_mean": float(np.mean(vals)),
                     "mcse": mcse,
                     "n_test_mean": float(np.mean(
                         [rep["methods"][m].bins[bi]["n_test"] for rep in repeats
                          if rep["methods"][m].ok and
                          bi in _eligible_bins(rep["methods"][m].bins)])),
                     "width_mean": float(np.mean(
                         [rep["methods"][m].bins[bi]["width"] for rep in repeats
                          if rep["methods"][m].ok and
                          bi in _eligible_bins(rep["methods"][m].bins)])),
                     "bias_mean": float(np.mean(
                         [rep["methods"][m].bins[bi]["bias"] for rep in repeats
                          if rep["methods"][m].ok and
                          bi in _eligible_bins(rep["methods"][m].bins)]))}
            if mcse is None:
                entry["status"] = DEFERRED
            out["per_bin"][m][str(bi)] = entry
    out["mondrian_deferred_bins_tally"] = {str(k): v for k, v in
                                           sorted(deferred_tally.items())}
    out["fit_failures"] = {
        "cqr": sum(1 for rep in repeats if not rep["methods"]["cqr"].ok),
    }
    return out


def decide_p5(marginal_pooled, level):
    """P5-90 / P5-95 (independent): pooled marginal coverage >= nominal - 3
    MCSE (one-sided). Each level binds independently; no combined status."""
    nominal = level
    r_eff = marginal_pooled.get("r_eff", 0)
    if r_eff < 2 or marginal_pooled.get("mcse") is None:
        return {"prediction": f"P5-{int(round(level * 100))}", "outcome": DEFERRED,
                "r_eff": r_eff, "detail": "marginal MCSE undefined"}
    mean = marginal_pooled["mean"]
    mcse = marginal_pooled["mcse"]
    threshold = nominal - P5_MCSE_K * mcse
    outcome = "HOLD" if mean >= threshold else "FAIL"
    return {"prediction": f"P5-{int(round(level * 100))}", "outcome": outcome,
            "mean_cov": mean, "mcse": mcse, "r_eff": r_eff,
            "threshold": threshold, "nominal": nominal}


def decide_p6(per_bin_pooled_090):
    """P6 (level 0.90, bound): >=2 pLDDT bins with |cov - 0.90| > 5 MCSE.

    Only bins with R_bin >= 2 count as evidence in EITHER direction; a bin
    with R_bin < 2 is typed DEFERRED for its own statistic. Fewer than 2
    bins carrying R_bin >= 2 -> the WHOLE P6 is DEFERRED at this
    evaluation (a section 6 non-outcome).
    """
    bins_supported, bins_deferred, deviations = [], [], {}
    for bi, st in sorted(per_bin_pooled_090.items()):
        if st.get("status") == DEFERRED or st.get("mcse") is None:
            bins_deferred.append(bi)
            continue
        bins_supported.append(bi)
        deviations[bi] = abs(st["cov_mean"] - 0.90)
    if len(bins_supported) < 2:
        return {"prediction": "P6", "outcome": DEFERRED,
                "bins_supported": bins_supported, "bins_deferred": bins_deferred,
                "detail": "fewer than 2 bins with R_bin >= 2"}
    deviating = [bi for bi in bins_supported
                 if deviations[bi] > P6_MCSE_K * per_bin_pooled_090[bi]["mcse"]]
    outcome = "REPLICATE" if len(deviating) >= 2 else "FALSIFY"
    return {"prediction": "P6", "outcome": outcome,
            "bins_supported": bins_supported, "bins_deferred": bins_deferred,
            "deviating_bins": deviating,
            "deviations": {str(bi): deviations[bi] for bi in bins_supported}}


def _pair_series(repeats, level, method):
    """(maxdev_r, eligible) pairs per repeat for one method at one level.

    maxdev_r = max over ELIGIBLE bins (test n >= 30 on the bound grid; a
    per-repeat max consuming no per-bin MCSE) of |coverage_bin,r - level|.
    Returns None for a FIT_FAILURE repeat.
    """
    series = []
    for rep in repeats:
        run = rep["methods"][method]
        if not run.ok:
            series.append(None)
            continue
        elig = _eligible_bins(run.bins)
        if not elig:
            series.append(None)  # no eligible bins: contributes to neither side
            continue
        maxdev = max(abs(run.bins[bi]["cov"] - level) for bi in elig)
        series.append(maxdev)
    return series


def decide_p7(repeats, level, method):
    """P7a/P7b: paired d_r = maxdev_r(method) - maxdev_r(pooled) over
    R_eff(pair) contributing repeats (both contribute AND >=1 bin
    eligible). Win iff mean(d) < -M AND |mean(d)| > 2 sd/sqrt(R_eff);
    mirrored for a pooled win; otherwise NO_DECISION; R_eff < 2 -> DEFERRED.
    """
    tag = f"P7{'a' if method == 'mondrian' else 'b'}"
    d_series = []
    a = _pair_series(repeats, level, method)
    b = _pair_series(repeats, level, "pooled")
    for da, db in zip(a, b):
        if da is None or db is None:
            continue
        d_series.append(da - db)
    r_eff = len(d_series)
    if r_eff < 2:
        return {"prediction": tag, "outcome": DEFERRED, "r_eff_pair": r_eff,
                "method": method, "detail": "R_eff(pair) < 2"}
    d = np.asarray(d_series)
    mean = float(d.mean())
    sd = float(d.std(ddof=1))
    half_width = CRIT_2SIGMA * sd / math.sqrt(r_eff)
    if mean < -MARGIN_M and abs(mean) > half_width:
        outcome = "WIN"
    elif mean > MARGIN_M and abs(mean) > half_width:
        outcome = "MIRRORED_WIN"
    else:
        outcome = "NO_DECISION"
    return {"prediction": tag, "outcome": outcome, "method": method,
            "mean_d": mean, "sd_d": sd, "half_width": half_width,
            "r_eff_pair": r_eff, "margin": MARGIN_M}


# --- whole evaluation (one section 6 checkpoint) ----------------------------


def run_evaluation(rows, seeds=range(200), levels=LEVELS):
    """Full registered evaluation over the cumulative label ledger.

    Returns the evaluation artifact: exclusions accounting, canonical
    accession list, per-level aggregates, and P5-90/P5-95/P6/P7a/P7b
    outcomes (typed statuses where support is insufficient). Binding
    across scheduled evaluations is the section 6 driver's role, not this
    function's.
    """
    accessions, by_acc, excluded = _prep_rows(rows)
    artifact = {
        "schema": "e422-evaluation-v1",
        "registration": REGISTRATION_CITATION,
        "n_rows_input": len(rows),
        "exclusions": excluded,
        "n_accessions": len(accessions),
        "canonical_accessions": accessions,
        "grid": [list(b) for b in _ANALYSIS_GRID],
        "seeds": [int(seeds[0]), int(seeds[-1])] if len(seeds) else [],
        "levels": list(levels),
    }
    if len(accessions) < 4:
        artifact["status"] = DEFERRED
        artifact["detail"] = ("pool_too_small_for_registered_split:"
                              f"n_accessions={len(accessions)}")
        return artifact
    # registered minimum-pool precondition, k arm (LAB-06): the conformal
    # order statistic k = ceil((n_cal+1)(1-alpha_c)) must satisfy
    # k <= n_cal against the repeat's CALIBRATION SCORE count (rows, not
    # accessions) at EVERY registered level and EVERY seed; a violation is
    # the typed whole-evaluation DEFERRED non-outcome (never a crash,
    # never silent).  Deterministic pre-screen: no fitting is needed.
    counts = [len(by_acc[a][0]) for a in accessions]
    n = len(accessions)
    for seed in seeds:
        _tr, ca_idx, _te = split_indices(n, int(seed))
        n_cal = sum(counts[i] for i in ca_idx)
        for lv in levels:
            k = math.ceil((n_cal + 1) * (1.0 - alpha_c_for(lv)))
            if k > n_cal:
                artifact["status"] = DEFERRED
                artifact["detail"] = (f"calibration_support_below_registered_"
                                      f"order_statistic:level={lv},k={k},"
                                      f"n_cal={n_cal},seed={int(seed)}")
                return artifact
    repeats_by_level = {lv: [] for lv in levels}
    for lv in levels:
        for seed in seeds:
            repeats_by_level[lv].append(run_repeat(rows, seed, lv))
    artifact["aggregates"] = {str(lv): aggregate_level(repeats_by_level[lv], lv)
                              for lv in levels}
    if 0.90 in repeats_by_level:
        # section 3B secondary: joint pLDDT x B-factor-z grid, 0.90 only,
        # same repeats as the primaries, pooled intervals (descriptive;
        # no joint cell is a primary endpoint, no joint floor claimed)
        artifact["joint_grid"] = aggregate_joint(repeats_by_level[0.90])
    p5 = {}
    for lv in levels:
        p5[str(lv)] = decide_p5(artifact["aggregates"][str(lv)]["marginal"]["pooled"], lv)
    p6 = decide_p6(artifact["aggregates"][str(0.90)]["per_bin"]["pooled"])
    p7a = decide_p7(repeats_by_level[0.90], 0.90, "mondrian")
    p7b = decide_p7(repeats_by_level[0.90], 0.90, "cqr")
    artifact["predictions"] = {"p5": p5, "p6": p6, "p7a": p7a, "p7b": p7b}
    artifact["status"] = "EVALUATED"
    return artifact
