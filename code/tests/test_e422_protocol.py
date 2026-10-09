"""Offline tests for the e422 frozen analysis protocol.

Binds experiments/e422_protocol.py to docs/e422_registration.md sections
3, 3A, 3B, 6, 7: split algorithm, shared alpha_c, per-method interval
constructions, sklearn pins, effective-repeat counts, and the P5/P6/P7
decision rules including typed deferral. No network, no artifacts.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments import e421_conformal  # noqa: E402
from experiments import e422_protocol as P  # noqa: E402
from sklearn.linear_model import QuantileRegressor  # noqa: E402


# --- registered constants ---------------------------------------------------


def test_grid_matches_committed_e421_source():
    """The listed six-bin values are authoritative and identical to the
    committed schema-v2 grid (registration section 3)."""
    assert P.PLDDT_BINS == e421_conformal.PLDDT_BINS
    assert P.PLDDT_BINS[-1] == (90.0, 100.01)


def test_shared_alpha_c_mapping():
    assert P.alpha_c_for(0.90) == pytest.approx(0.10)
    assert P.alpha_c_for(0.95) == pytest.approx(0.05)
    for lv in P.LEVELS:
        assert 0.0 < P.alpha_c_for(lv) < 1.0


def test_bin_index_grid_membership():
    assert P.bin_index(0.0) == 0
    assert P.bin_index(49.999) == 0
    assert P.bin_index(50.0) == 1
    assert P.bin_index(60.0) == 2
    assert P.bin_index(70.0) == 3
    assert P.bin_index(80.0) == 4
    assert P.bin_index(90.0) == 5
    assert P.bin_index(100.0) == 5  # half-open grid admits exact 100
    with pytest.raises(P.E422SpecViolation):
        P.bin_index(100.01)
    with pytest.raises(P.E422SpecViolation):
        P.bin_index(-0.5)


# --- split algorithm (section 3, exact) -------------------------------------


def test_split_indices_matches_registered_algorithm():
    n, seed = 11, 7
    tr, ca, te = P.split_indices(n, seed)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)  # the registered only draw
    t = n // 2
    c = (n - t) // 2
    assert np.array_equal(tr, perm[:t])
    assert np.array_equal(ca, perm[t:t + c])
    assert np.array_equal(te, perm[t + c:])
    assert tr.size == 5 and ca.size == 3 and te.size == 3  # 11 -> 5/3/3


def test_split_disjoint_and_complete_for_various_n():
    for n in (4, 5, 6, 7, 13, 50, 97):
        parts = P.split_indices(n, 0)
        assert sum(p.size for p in parts) == n
        joined = np.concatenate(parts)
        assert sorted(joined.tolist()) == list(range(n))
        for p in parts:
            assert p.size == len(set(p.tolist()))


def test_split_is_deterministic_per_seed():
    for seed in (0, 3, 199):
        a = P.split_indices(23, seed)
        b = P.split_indices(23, seed)
        assert all(np.array_equal(x, y) for x, y in zip(a, b))


def test_canonical_accessions_is_sorted_unique():
    rows = [{"accession": "B"}, {"accession": "A"}, {"accession": "B"}]
    assert P.canonical_accessions(rows) == ["A", "B"]


# --- conformal order statistic ----------------------------------------------


def test_order_statistic_ceil_n_plus_1_times_coverage():
    scores = np.arange(1.0, 11.0)  # 1..10
    # alpha_c = 0.10 -> k = ceil(11 * 0.90) = 10 -> 10th smallest = 10.0
    assert P.conformal_order_stat(scores, 0.10) == 10.0
    # alpha_c = 0.50 -> k = ceil(11 * 0.50) = 6 -> 6.0
    assert P.conformal_order_stat(scores, 0.50) == 6.0


def test_order_statistic_beyond_n_is_spec_violation():
    # alpha_c is the CONFORMAL LEVEL: k = ceil((n+1)(1-alpha_c)).
    # n=2, alpha_c=0.10 -> k = ceil(3*0.9) = 3 > 2 -> spec violation.
    with pytest.raises(P.E422SpecViolation):
        P.conformal_order_stat(np.array([1.0, 2.0]), 0.10)
    with pytest.raises(P.E422SpecViolation):
        P.conformal_order_stat(np.array([]), 0.10)


# --- method constructions (section 3A) ---------------------------------------


def _mk_arrays(n_per_bin, y_of_plddt, noise_rng=None):
    p, y = [], []
    for plddt_val, count in n_per_bin.items():
        for i in range(count):
            p.append(plddt_val)
            y.append(y_of_plddt(plddt_val, i, noise_rng))
    return (np.array(p, dtype=float), np.array(y, dtype=float))


def test_pooled_interval_construction_exact():
    train_y = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    center = float(np.median(train_y))
    cal_y = np.array([center + 0.5, center - 1.0, center + 2.0, center - 3.0])
    test_y = np.array([0.0, 100.0])
    # alpha_c = 0.5 -> k = ceil(5*0.5) = 3 -> 3rd smallest score = 2.0
    c, q, lo, hi = P.run_pooled(train_y, cal_y, test_y, alpha_c=0.50)
    assert c == 3.0
    assert q == pytest.approx(2.0)
    assert np.allclose(lo, np.full(2, 1.0))
    assert np.allclose(hi, np.full(2, 5.0))


def test_mondrian_defers_under_supported_bins_to_pooled():
    rng = np.random.default_rng(0)
    n = 40
    rich_p = np.full(n, 95.0)  # bin 5 well supported
    rich_y = rng.normal(80.0, 2.0, n)
    poor_p = np.full(5, 10.0)  # bin 0 starved (train/cal n < 30)
    poor_y = rng.normal(20.0, 2.0, 5)
    tr_p = np.concatenate([rich_p[:20], poor_p[:3]])
    tr_y = np.concatenate([rich_y[:20], poor_y[:3]])
    ca_p = np.concatenate([rich_p[20:30], poor_p[3:]])
    ca_y = np.concatenate([rich_y[20:30], poor_y[3:]])
    te_p = rich_p[30:]
    te_y = rich_y[30:]
    pooled_center, pooled_q, deferred, lo, hi = P.run_mondrian(
        tr_p, tr_y, ca_p, ca_y, te_p, te_y, alpha_c=0.10)
    # starved bin 0 defers; MONDRIAN_MIN_N=30 means bin 5 (n_k=20) ALSO
    # defers at these tiny n: both bins fall back to the pooled interval.
    assert 0 in deferred and 5 in deferred
    pc, pq, p_lo, p_hi = P.run_pooled(tr_y, ca_y, te_y, 0.10)
    assert np.allclose(lo, p_lo)
    assert np.allclose(hi, p_hi)


def test_mondrian_uses_per_bin_center_when_supported():
    rng = np.random.default_rng(1)
    def block(plddt_val, base):
        p = np.full(40, plddt_val)
        y = base + rng.normal(0, 1.0, 40)
        return p, y
    p1, y1 = block(95.0, 80.0)   # bin 5
    p2, y2 = block(10.0, 20.0)   # bin 0
    tr_p = np.concatenate([p1[:20], p2[:20]])
    tr_y = np.concatenate([y1[:20], y2[:20]])
    ca_p = np.concatenate([p1[20:], p2[20:]])
    ca_y = np.concatenate([y1[20:], y2[20:]])
    te_p = np.array([95.0, 10.0])
    te_y = np.array([80.0, 20.0])
    # still deferred: n_k = 20 < MONDRIAN_MIN_N = 30 in BOTH splits
    _, _, deferred, lo, hi = P.run_mondrian(tr_p, tr_y, ca_p, ca_y, te_p, te_y, 0.10)
    assert 0 in deferred and 5 in deferred


def test_mondrian_supported_bin_gets_own_interval():
    rng = np.random.default_rng(2)
    def block(plddt_val, base, n):
        return np.full(n, plddt_val), base + rng.normal(0, 0.5, n)
    p1, y1 = block(95.0, 80.0, 80)   # bin 5: 40 train / 40 cal
    p2, y2 = block(10.0, 20.0, 80)   # bin 0: 40 train / 40 cal
    tr_p = np.concatenate([p1[:40], p2[:40]])
    tr_y = np.concatenate([y1[:40], y2[:40]])
    ca_p = np.concatenate([p1[40:], p2[40:]])
    ca_y = np.concatenate([y1[40:], y2[40:]])
    te_p = np.array([95.0, 10.0])
    te_y = np.array([80.0, 20.0])
    _, _, deferred, lo, hi = P.run_mondrian(tr_p, tr_y, ca_p, ca_y, te_p, te_y, 0.10)
    # bins 0 and 5 are supported on both splits; empty bins 1-4 defer
    # (n_k = 0 < 30) but hold no test rows
    assert 0 not in deferred and 5 not in deferred
    assert all(b in deferred for b in (1, 2, 3, 4))
    # bin 5 interval centers near 80, bin 0 near 20 (per-bin centers, not pooled)
    assert 70.0 < lo[0] and hi[0] < 95.0
    assert 10.0 < lo[1] and hi[1] < 35.0


def test_cqr_sklearn_pins_match_registration():
    import inspect
    sig = inspect.signature(QuantileRegressor)
    assert "quantile" in sig.parameters
    assert sig.parameters["fit_intercept"].default is True
    assert sig.parameters["alpha"].default == 1.0  # we explicitly pass 0.0
    # registered construction: alpha = 0.0 (L1 weight), solver highs
    m = QuantileRegressor(alpha=0.0, solver="highs", fit_intercept=True,
                          quantile=0.9)
    assert m.alpha == 0.0 and m.solver == "highs" and m.quantile == 0.9


def test_cqr_design_has_exactly_one_constant_column():
    """[pLDDT, pLDDT^2] with fit_intercept=True: two coefs + one intercept,
    NOT [1, p, p^2] which would duplicate the intercept."""
    rng = np.random.default_rng(3)
    n = 300
    p = np.clip(rng.normal(70.0, 15.0, n), 1.0, 99.0)
    y = 10.0 + 0.9 * p + rng.normal(0, 1.0, n)
    tr = slice(0, 200)
    ca = slice(200, 250)
    te = slice(250, 300)
    lo, hi = P.run_cqr(p[tr], y[tr], p[ca], y[ca], p[te], y[te], alpha_c=0.10)
    assert lo.shape == (50,) and hi.shape == (50,)
    # on near-deterministic linear data the band must be tight and ordered
    assert np.all(lo < hi)
    assert np.all(lo <= y[te] + 5.0) and np.all(hi >= y[te] - 5.0)


def test_cqr_tail_assignment_lower_upper():
    rng = np.random.default_rng(4)
    n = 400
    p = np.linspace(1.0, 99.0, n)
    y = 3.0 * p + 5.0 + rng.normal(0, 4.0, n)  # noise separates the tails
    tr, ca, te = slice(0, 300), slice(300, 350), slice(350, 400)
    lo, hi = P.run_cqr(p[tr], y[tr], p[ca], y[ca], p[te], y[te], alpha_c=0.10)
    # q_lo fitted at tau_lo = alpha_c must sit BELOW q_hi fitted at 1-alpha_c
    # (interval [q_lo - Q, q_hi + Q] brackets the line y = 3p + 5)
    assert np.all(lo <= 3.0 * p[te] + 5.0 + 15.0)
    assert np.all(hi >= 3.0 * p[te] + 5.0 - 15.0)
    assert np.all(lo < hi)


def test_cqr_fit_failure_is_typed(monkeypatch):
    calls = {"n": 0}

    def boom(*a, **k):
        calls["n"] += 1
        raise P._FitFailure("forced_for_test")

    monkeypatch.setattr(P, "_fit_cqr_quantile", boom)
    rng = np.random.default_rng(5)
    n = 60
    p = rng.uniform(1, 99, n)
    y = rng.normal(50, 5, n)
    run = P.run_repeat(
        [{"accession": f"A{i // 10}", "plddt": float(p[i]), "lddt": float(y[i]),
          "d_kabsch": 1.0, "b_factor_z": 0.0} for i in range(n)],
        seed=0, level=0.90)
    cqr = run["methods"]["cqr"]
    assert calls["n"] >= 1
    assert not cqr.ok and cqr.fit_failure
    # other methods unaffected (paired on the same split)
    assert run["methods"]["pooled"].ok and run["methods"]["mondrian"].ok


def test_isotonic_conformal_step_is_pooled_construction():
    rng = np.random.default_rng(6)
    n = 200
    p = np.sort(rng.uniform(1, 99, n))
    y = 0.8 * p + rng.normal(0, 2.0, n)
    tr, ca, te = slice(0, 120), slice(120, 160), slice(160, 200)
    lo, hi = P.run_isotonic(p[tr], y[tr], p[ca], y[ca], p[te], y[te], alpha_c=0.10)
    from sklearn.isotonic import IsotonicRegression
    iso = IsotonicRegression(out_of_bounds="clip").fit(p[tr], y[tr])
    scores = np.abs(y[ca] - iso.predict(p[ca]))
    q = P.conformal_order_stat(scores, 0.10)
    assert np.allclose(lo, iso.predict(p[te]) - q)
    assert np.allclose(hi, iso.predict(p[te]) + q)


def test_isotonic_out_of_bounds_clip():
    from sklearn.isotonic import IsotonicRegression
    iso = IsotonicRegression(out_of_bounds="clip")
    p = np.array([10.0, 20.0, 30.0, 40.0])
    y = np.array([1.0, 2.0, 3.0, 4.0])
    iso.fit(p, y)
    assert iso.predict(np.array([0.0, 50.0]))[0] == pytest.approx(1.0)
    assert iso.predict(np.array([0.0, 50.0]))[1] == pytest.approx(4.0)


# --- effective repeats, MCSE, decision rules ---------------------------------


def test_mcse_deferred_below_two_contributors():
    assert P._mcse([0.5]) == (None, 1)
    assert P._mcse([]) == (None, 0)
    mcse, n = P._mcse([0.1, 0.3, 0.5, 0.7])
    assert n == 4 and mcse == pytest.approx(np.std([0.1, 0.3, 0.5, 0.7], ddof=1) / 2.0)


def test_decide_p5_thresholds():
    held = P.decide_p5({"r_eff": 200, "mean": 0.905, "mcse": 0.005}, 0.90)
    assert held["outcome"] == "HOLD"  # 0.905 >= 0.90 - 0.015
    failed = P.decide_p5({"r_eff": 200, "mean": 0.880, "mcse": 0.005}, 0.90)
    assert failed["outcome"] == "FAIL"
    deferred = P.decide_p5({"r_eff": 1, "mean": 0.9, "mcse": None}, 0.90)
    assert deferred["outcome"] == P.DEFERRED
    # levels are independent: 0.95 evaluated against ITS OWN nominal
    held95 = P.decide_p5({"r_eff": 200, "mean": 0.945, "mcse": 0.004}, 0.95)
    assert held95["outcome"] == "HOLD" and held95["prediction"] == "P5-95"


def test_decide_p6_counting_rule():
    table = {
        "0": {"r_bin": 200, "cov_mean": 0.60, "mcse": 0.01},   # deviates: |0.6-0.9|=0.30 > 0.05
        "1": {"r_bin": 200, "cov_mean": 0.89, "mcse": 0.01},   # no deviation
        "2": {"r_bin": 200, "cov_mean": 0.50, "mcse": 0.02},   # deviates: 0.40 > 0.10
        "3": {"r_bin": 1, "cov_mean": 0.10, "mcse": None, "status": P.DEFERRED},
    }
    res = P.decide_p6(table)
    assert res["outcome"] == "REPLICATE"
    assert sorted(res["deviating_bins"]) == ["0", "2"]
    assert "3" in res["bins_deferred"]  # R_bin < 2 is evidence in NEITHER direction

    one_dev = dict(table)
    one_dev["2"] = {"r_bin": 200, "cov_mean": 0.90, "mcse": 0.01}
    assert P.decide_p6(one_dev)["outcome"] == "FALSIFY"

    under = {k: v for k, v in table.items() if k in ("0", "3")}
    assert P.decide_p6(under)["outcome"] == P.DEFERRED


def test_decide_p7_win_mirror_nodecision_deferred():
    # WIN: mean(d) < -M and |mean| > 2 sd/sqrt(R)
    d_win = np.array([-0.05, -0.06, -0.05, -0.04, -0.05])
    mean, sd = d_win.mean(), d_win.std(ddof=1)
    assert mean < -P.MARGIN_M
    assert abs(mean) > 2.0 * sd / math.sqrt(len(d_win))
    # NO_DECISION: small mean
    d_small = np.array([-0.001, 0.002, -0.001, 0.001])
    assert abs(d_small.mean()) < P.MARGIN_M
    # MIRRORED_WIN: mean(d) > M with the same significance
    d_mirror = -d_win
    assert d_mirror.mean() > P.MARGIN_M

    # via the registered path with paired repeats:
    def mk_rep(m_cov, p_cov):
        bins = {bi: {"n_test": 100, "cov": c, "width": 1.0, "bias": 0.0}
                for bi, c in m_cov.items()}
        pbins = {bi: {"n_test": 100, "cov": c, "width": 1.0, "bias": 0.0}
                 for bi, c in p_cov.items()}
        m = P.MethodRun("mondrian", True, bins=bins)
        p = P.MethodRun("pooled", True, bins=pbins)
        cqr = P.MethodRun("cqr", True, bins=bins)
        iso = P.MethodRun("isotonic", True, bins=pbins)
        return {"seed": 0, "level": 0.90, "n_accessions": 10,
                "split_sizes": {"train": 1, "cal": 1, "test": 1},
                "methods": {"pooled": p, "mondrian": m, "cqr": cqr, "isotonic": iso}}

    # mondrian much better in every eligible bin -> d_r < 0 in every repeat
    reps_win = [mk_rep({0: 0.90}, {0: 0.50}) for _ in range(5)]
    res = P.decide_p7(reps_win, 0.90, "mondrian")
    assert res["prediction"] == "P7a" and res["outcome"] == "WIN"
    # mirrored: pooled better everywhere
    res_m = P.decide_p7(reps_win, 0.90, "cqr")  # cqr tracks mondrian covs here
    assert res_m["prediction"] == "P7b" and res_m["outcome"] == "WIN"
    reps_mirror = [mk_rep({0: 0.50}, {0: 0.90}) for _ in range(5)]
    assert P.decide_p7(reps_mirror, 0.90, "mondrian")["outcome"] == "MIRRORED_WIN"
    # NO_DECISION: tiny paired differences
    reps_tiny = [mk_rep({0: 0.90 + 1e-4 * i}, {0: 0.90}) for i in range(5)]
    assert P.decide_p7(reps_tiny, 0.90, "mondrian")["outcome"] == "NO_DECISION"
    # DEFERRED: R_eff(pair) < 2
    assert P.decide_p7(reps_win[:1], 0.90, "mondrian")["outcome"] == P.DEFERRED


def test_pair_series_drops_fit_failures_and_empty_eligibility():
    def rep(ok, n_test, cov):
        bins = {0: {"n_test": n_test, "cov": cov, "width": 1.0, "bias": 0.0}}
        m = P.MethodRun("mondrian", ok, fit_failure=not ok, bins=bins)
        p = P.MethodRun("pooled", True, bins=bins)
        c = P.MethodRun("cqr", ok, fit_failure=not ok, bins=bins)
        i = P.MethodRun("isotonic", True, bins=bins)
        return {"seed": 0, "level": 0.90, "n_accessions": 10,
                "split_sizes": {"train": 1, "cal": 1, "test": 1},
                "methods": {"pooled": p, "mondrian": m, "cqr": c, "isotonic": i}}

    reps = [rep(True, 100, 0.8), rep(False, 100, 0.5), rep(True, 10, 0.1)]
    a = P._pair_series(reps, 0.90, "mondrian")
    assert a[0] == pytest.approx(0.10)
    assert a[1] is None            # FIT_FAILURE repeat: dropped
    assert a[2] is None            # no eligible bin (test n < 30): dropped


# --- end-to-end on synthetic ledgers -----------------------------------------


def _synthetic_ledger(n_acc=40, rows_per_acc=40, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    plddt_pool = np.concatenate([
        np.full(6, rng.uniform(5, 45)),      # bin 0
        np.full(6, rng.uniform(50, 59)),     # bin 1
        np.full(6, rng.uniform(60, 69)),     # bin 2
        np.full(8, rng.uniform(70, 79)),     # bin 3
        np.full(8, rng.uniform(80, 89)),     # bin 4
        np.full(6, rng.uniform(90, 99)),     # bin 5
    ])
    for a in range(n_acc):
        acc = f"AAAA{a:02d}"
        base = 0.15 + 0.008 * rng.uniform(0, 1)
        for _ in range(rows_per_acc):
            p = float(rng.choice(plddt_pool))
            y = float(np.clip(0.2 + 0.007 * p + rng.normal(0, 0.05), 0.0, 1.0))
            rows.append({"accession": acc, "plddt": p, "lddt": y,
                         "d_kabsch": abs(rng.normal(1.0, 0.5)),
                         "b_factor_z": rng.normal(0, 1)})
    return rows


def test_run_repeat_pairs_all_methods_on_same_split():
    rows = _synthetic_ledger(n_acc=12, rows_per_acc=30)
    run = P.run_repeat(rows, seed=0, level=0.90)
    assert run["n_accessions"] == 12
    assert sum(run["split_sizes"].values()) == 12
    for m in ("pooled", "mondrian", "cqr", "isotonic"):
        assert run["methods"][m].ok
        assert run["methods"][m].n_test > 0
        assert run["methods"][m].bins  # per-bin n always reported


def test_run_evaluation_end_to_end_and_exclusions():
    rows = _synthetic_ledger()
    rows.append({"accession": "AAAA00", "plddt": 70.0, "lddt": None})   # null label
    rows.append({"accession": "AAAA00", "plddt": None, "lddt": 0.5})   # null feature
    rows.append({"accession": "AAAA00", "plddt": 150.0, "lddt": 0.5})  # out of grid
    art = P.run_evaluation(rows, seeds=range(3), levels=(0.90,))
    assert art["status"] == "EVALUATED"
    assert art["exclusions"] == {"null_label": 1, "nonfinite_label": 0,
                                 "null_feature": 1, "nonfinite_feature": 0,
                                 "feature_out_of_grid": 1,
                                 "null_descriptive_channel": 0}
    assert art["canonical_accessions"] == sorted(art["canonical_accessions"])
    agg = art["aggregates"]["0.9"]
    assert "r_eff" in agg and "per_bin" in agg
    preds = art["predictions"]
    assert set(preds) == {"p5", "p6", "p7a", "p7b"}
    outcomes = ("HOLD", "FAIL", "REPLICATE", "FALSIFY",
                "WIN", "MIRRORED_WIN", "NO_DECISION", P.DEFERRED)
    # P5 exists per level, independently (no combined P5 status)
    assert set(preds["p5"]) == {"0.9"}
    for lv_res in preds["p5"].values():
        assert lv_res["prediction"] == "P5-90"
        assert lv_res["outcome"] in outcomes
    for key in ("p6", "p7a", "p7b"):
        assert preds[key]["outcome"] in outcomes


def test_run_evaluation_deferred_when_pool_too_small():
    rows = [{"accession": "A", "plddt": 50.0, "lddt": 0.5},
            {"accession": "B", "plddt": 60.0, "lddt": 0.6}]
    art = P.run_evaluation(rows, seeds=range(2), levels=(0.90,))
    assert art["status"] == P.DEFERRED


def test_conformal_marginal_coverage_is_approximately_nominal():
    """Statistical sanity: pooled split conformal on exchangeable synthetic
    data must land near nominal marginal coverage (development-only data;
    no registered number is computed here)."""
    rows = _synthetic_ledger(n_acc=30, rows_per_acc=40, seed=7)
    art = P.run_evaluation(rows, seeds=range(10), levels=(0.90,))
    pooled = art["aggregates"]["0.9"]["marginal"]["pooled"]
    assert 0.80 <= pooled["mean"] <= 1.0
    assert pooled["r_eff"] == 10


# --- gate-(c) disposition tests (findings e422-c-01/-02/-03, C-01/-02/-03) --

def test_prep_rows_excludes_nonfinite_label():
    rows = [{"accession": "A", "plddt": 70.0, "lddt": float("nan"),
             "d_kabsch": 1.0, "b_factor_z": 0.0},
            {"accession": "A", "plddt": 70.0, "lddt": 0.5,
             "d_kabsch": 1.0, "b_factor_z": 0.0}]
    accessions, by_acc, excluded = P._prep_rows(rows)
    assert accessions == ["A"]
    assert excluded["nonfinite_label"] == 1
    assert len(by_acc["A"][0]) == 1


def test_prep_rows_excludes_missing_descriptive_channel():
    rows = [{"accession": "A", "plddt": 70.0, "lddt": 0.5,
             "d_kabsch": None, "b_factor_z": 0.0},
            {"accession": "A", "plddt": 70.0, "lddt": 0.5,
             "d_kabsch": 1.0, "b_factor_z": None}]
    _acc, _by, excluded = P._prep_rows(rows)
    assert excluded["null_descriptive_channel"] == 2
    # feature/grid problems take precedence over descriptive-channel checks
    rows2 = [{"accession": "A", "plddt": 150.0, "lddt": 0.5,
              "d_kabsch": None, "b_factor_z": None}]
    _a, _b, ex2 = P._prep_rows(rows2)
    assert ex2["feature_out_of_grid"] == 1
    assert ex2["null_descriptive_channel"] == 0


def test_k_above_ncal_returns_whole_evaluation_deferred():
    # the stat-audit demonstration case: 4 accessions x 2 rows ->
    # calibration 2 scores, k=3 > 2 at 0.90 -> typed DEFERRED, no crash
    rows = [{"accession": f"A{a}", "plddt": 70.0 + a, "lddt": 0.5 + 0.01 * a,
             "d_kabsch": 1.0, "b_factor_z": 0.0}
            for a in range(4) for _ in range(2)]
    art = P.run_evaluation(rows, seeds=range(3), levels=(0.90,))
    assert art["status"] == P.DEFERRED
    assert "calibration_support_below_registered_order_statistic" in art["detail"]


def test_registration_citation_is_frozen_gate_a_object():
    # gate (a) 2026-09-09: BOTH identifiers labeled (git-blob cat-file
    # resolvable; content-sha256 of the identical bytes) — finding 3B-4.
    assert "git-blob 6e1b9ce128ac3a1be4678463bdc777c25b50ca8d" \
        in P.REGISTRATION_CITATION
    assert "content-sha256 dc7a362773fbf3a03e759aa7d178125963a14d88d9bbd" \
        "61b66389a04ba8d7b3b" in P.REGISTRATION_CITATION
    assert "9dd2380b" in P.REGISTRATION_CITATION
    assert "FROZEN gate-(a)" in P.REGISTRATION_CITATION
    art = P.run_evaluation(
        _synthetic_ledger(n_acc=40, rows_per_acc=40, seed=3),
        seeds=range(3), levels=(0.90,))
    assert art["registration"] == P.REGISTRATION_CITATION


def test_bias_sign_is_registered_y_minus_center():
    # section 3D: mean signed residual y - method center (descriptive)
    stats = P._bin_stats(np.array([85.0, 85.0]), np.array([0.9, 0.9]),
                         np.array([0.0, 0.0]), np.array([0.5, 0.5]))
    assert stats[4]["bias"] == pytest.approx(0.9 - 0.25)


# --- section 3B joint grid (pLDDT bin x B-factor-z band) -------------------

def test_bz_band_boundaries_are_registered_inequality():
    assert P._bz_band(-1.0) == 0     # <= -1 is band 0 (inclusive)
    assert P._bz_band(-2.5) == 0
    assert P._bz_band(-0.999) == 1   # (-1, +1) open on both ends
    assert P._bz_band(0.0) == 1
    assert P._bz_band(0.999) == 1
    assert P._bz_band(1.0) == 2      # >= +1 is band 2 (inclusive)
    assert P._bz_band(3.0) == 2


def test_joint_cells_counts_coverage_and_proteins():
    # 4 residues in bin 3 (70-80), band 1: two from acc X, one from acc Y,
    # one uncovered -> n_rows=4, n_prot=2, cov=0.75
    p = np.array([75.0] * 4)
    y = np.array([0.5, 0.6, 0.7, 0.8])
    z = np.array([0.0, 0.2, -0.5, 0.9])
    acc = ["X", "X", "Y", "Z"]
    lo = np.full(4, 0.0)
    hi = np.full(4, 0.75)
    cells = P._joint_cells(p, y, z, acc, lo, hi)
    c = cells["3:1"]
    assert c["n_rows"] == 4
    assert c["n_prot"] == 3
    assert c["cov"] == pytest.approx(0.75)
    # z = -1 exactly is band 0, z = +1 exactly is band 2 (registered bands)
    z = np.array([-1.0, 1.0])
    cells = P._joint_cells(np.array([75.0, 75.0]), np.array([0.5, 0.5]), z,
                           ["X", "Y"], np.zeros(2), np.ones(2))
    assert cells["3:0"]["n_rows"] == 1
    assert cells["3:2"]["n_rows"] == 1
    assert cells["3:1"]["n_rows"] == 0


def test_joint_cells_inherit_bound_grid():
    rows = _synthetic_ledger(n_acc=8, rows_per_acc=30, seed=3)
    accessions, by_acc, _ = P._prep_rows(rows)
    p, y, z, acc = [], [], [], []
    for a in accessions:
        ps, ys, zs = by_acc[a]
        p.extend(ps); y.extend(ys); z.extend(zs); acc.extend([a] * len(ps))
    p, y, z = map(np.array, (p, y, z))
    restore = P.use_grid(((0.0, 80.0), (80.0, 90.0), (90.0, 100.01)))
    try:
        cells = P._joint_cells(p, y, z, acc, np.zeros(len(y)), np.ones(len(y)))
        assert len(cells) == 3 * 3  # merged-step-2 grid: 3 pLDDT bins x 3 bands
        assert {k.split(":")[0] for k in cells} == {"0", "1", "2"}
    finally:
        restore()


def test_aggregate_joint_registered_rule():
    def rep(n_rows, n_prot, cov):
        return {"joint_cells": {"2:1": {"n_rows": n_rows, "n_prot": n_prot,
                                        "cov": cov}}}
    # eligible in exactly 2 repeats -> REPORTED with MCSE over those 2
    repeats = [rep(150, 25, 0.9), rep(120, 21, 0.8), rep(10, 3, None)]
    out = P.aggregate_joint(repeats)
    entry = out["2:1"]
    assert entry["r_joint"] == 2
    assert entry["cov_mean"] == pytest.approx(0.85)
    assert entry["mcse"] is not None
    assert "status" not in entry
    assert entry["n_rows_mean_all_repeats"] == pytest.approx(280 / 3)
    # eligible in only 1 repeat -> UNDER_SUPPORTED_STRATUM, no coverage claim
    repeats = [rep(150, 25, 0.9), rep(50, 10, 0.5), rep(10, 3, None)]
    out = P.aggregate_joint(repeats)
    entry = out["2:1"]
    assert entry["status"] == "UNDER_SUPPORTED_STRATUM"
    assert "cov_mean" not in entry and "mcse" not in entry
    assert entry["r_joint"] == 1


def test_joint_cells_protein_floor_uses_distinct_accessions():
    # 150 rows but ALL from one accession -> n_prot=1 < 20 -> ineligible
    def rep(n_rows, n_prot):
        return {"joint_cells": {"0:0": {"n_rows": n_rows, "n_prot": n_prot,
                                        "cov": 0.9}}}
    out = P.aggregate_joint([rep(150, 1), rep(150, 1)])
    assert out["0:0"]["status"] == "UNDER_SUPPORTED_STRATUM"


def test_run_evaluation_emits_joint_grid_at_090_only():
    rows = _synthetic_ledger(n_acc=40, rows_per_acc=40, seed=1)
    art = P.run_evaluation(rows, seeds=range(2), levels=(0.90, 0.95))
    assert art["status"] == "EVALUATED"
    assert "joint_grid" in art
    assert len(art["joint_grid"]) == 6 * 3  # six-bin grid x 3 B-z bands
    # small synthetic ledger -> every cell UNDER_SUPPORTED_STRATUM (the
    # registered surface for low support; no coverage claim anywhere)
    assert all(e["status"] == "UNDER_SUPPORTED_STRATUM"
               for e in art["joint_grid"].values())
    assert all("cov_mean" not in e for e in art["joint_grid"].values())


def test_batch4_scope_strips_joint_grid():
    import copy
    from experiments import e422_evaldriver as drv
    ev = {"predictions": {"p5": {}, "p6": 1}, "joint_grid": {"0:0": {}},
          "aggregates": {}}
    scoped = drv._scope_batch4(copy.deepcopy(ev))
    assert "joint_grid" not in scoped
    assert "p6" not in scoped["predictions"]


def test_bz_band_derives_from_single_vector_source():
    # single source of truth: scalar derives from the vectorized helper
    grid = np.array([-3.0, -1.0, -0.999, 0.0, 0.999, 1.0, 5.0])
    vec = P._bz_bands_vec(grid)
    assert list(vec) == [0, 0, 1, 1, 1, 2, 2]
    assert all(P._bz_band(float(z)) == int(b) for z, b in zip(grid, vec))


def test_aggregate_joint_empty_repeats_is_typed_violation():
    with pytest.raises(P.E422SpecViolation):
        P.aggregate_joint([])


def test_joint_cells_computed_at_090_only():
    rows = _synthetic_ledger(n_acc=12, rows_per_acc=30, seed=5)
    run90 = P.run_repeat(rows, seed=0, level=0.90)
    run95 = P.run_repeat(rows, seed=0, level=0.95)
    assert len(run90["joint_cells"]) == 6 * 3
    assert run95["joint_cells"] == {}
