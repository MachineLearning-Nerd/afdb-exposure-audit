"""Revision-3 figures: F1 mechanism schematic, F2 defect scope, F3 exposure.

F2 implements REGISTRATION.md Addendum A4 "F-scope" (pLDDT-free): the pLDDT
field is dropped from every e422/e427 ledger row on load and never used.
F1 numbers are read from the 9qj6 audit JSON and from out/v7_trigger_prevalence.json
(which re-derives them from the raw SIFTS response and reference mmCIF bytes);
the script asserts they agree before drawing.

Usage: python -I make_figures_rev3.py <repo_root>
Outputs: paper/arxiv/figures/rev3/{f1_mechanism,f2_scope,f3_exposure}.{pdf,png}
         paper/arxiv/verification/out/figures_rev3_numbers.json
"""
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

INK, INK2, MUTED, GRID, SURF = "#1a1a19", "#52514e", "#8a8985", "#e4e3df", "#ffffff"
BLUE, ORANGE = "#2a78d6", "#eb6834"          # validated pair (CVD dE 24.7)
SEQ = ["#12366b", "#3d86de", "#9cc3f2", "#e4e3df"]  # ordinal: strongest -> none (lightness-monotone)
DELTA, ERR = 0.2, 0.60

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False,
    "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02})


def save(fig, fdir, name):
    fig.savefig(fdir / f"{name}.pdf")
    fig.savefig(fdir / f"{name}.png", dpi=300)
    plt.close(fig)


def load_rows(path):
    rows = []
    for r in json.loads(path.read_text())["rows"]:
        r = dict(r)
        r.pop("plddt", None)            # registered blinding (P6/P7): never read
        rows.append(r)
    return rows


def ledger(root, sub):
    rows = []
    for b in (1, 2, 3, 4):
        rows += load_rows(root / f"results/{sub}/batches/e422_batch0{b}_checkpoint.json")
    return rows


def finite(v):
    return v is not None and v == v


def seg_medians(rows):
    acc = defaultdict(list)
    for r in rows:
        if finite(r.get("lddt")):
            acc[(r["entry"], r["chain"], r["accession"])].append(float(r["lddt"]))
    return {k: statistics.median(v) for k, v in acc.items()}


# --------------------------------------------------------------------------- F1
def f1_numbers(root, vdir):
    audit = json.loads((root / "paper/arxiv/verification/audit/k7pq54_9qj6_audit.json").read_text())
    a = audit["chains"]["A"]
    s = a["sifts"]
    first = a["per_residue"][0]
    v7 = json.loads((vdir / "v7_trigger_prevalence.json").read_text())
    seg = next(t for t in v7["triggered_segments"] if t["entry"] == "9qj6" and t["chain_id"] == "A")
    n = {
        "entry": "9qj6", "chain": "A", "accession": "K7PQ54",
        "sifts_start_label_seq_id": s["start"]["residue_number"],
        "sifts_start_author_residue_number": s["start"]["author_residue_number"],
        "sifts_unp_start": s["unp_start"], "sifts_unp_end": s["unp_end"],
        "sifts_end_label_seq_id": s["end"]["residue_number"],
        "sifts_end_author_residue_number": s["end"]["author_residue_number"],
        "first_observed_label_seq_id": first["pdb_label_seq"],
        "first_observed_auth_seq_id": first["pdb_auth_seq"],
        "first_observed_correct_uniprot": first["uniprot_pos"],
        "ledger_first_uniprot_original": audit["e422_checkpoint_matches_summary"]["A"]["uniprot_range"][0],
    }
    n["author_offset"] = n["first_observed_auth_seq_id"] - n["first_observed_label_seq_id"]
    # faulty fallback: keys start at the label start, atoms looked up by author number
    n["faulty_uniprot_for_first_atom"] = (n["sifts_unp_start"] + n["first_observed_auth_seq_id"]
                                          - n["sifts_start_label_seq_id"])
    n["register_shift"] = n["faulty_uniprot_for_first_atom"] - n["first_observed_correct_uniprot"]
    # cross-checks against the raw-byte derivation in V7 and against the original ledger
    oc = seg["offset_check"]
    assert n["sifts_start_author_residue_number"] is None
    assert (oc["first_label_seq_id"], oc["first_auth_seq_id"], oc["offset"]) == (
        n["first_observed_label_seq_id"], n["first_observed_auth_seq_id"], n["author_offset"])
    assert seg["start_label"] == n["sifts_start_label_seq_id"] and seg["unp_start"] == n["sifts_unp_start"]
    assert n["faulty_uniprot_for_first_atom"] == n["ledger_first_uniprot_original"]
    assert n["first_observed_correct_uniprot"] == n["sifts_unp_start"] + (
        n["first_observed_label_seq_id"] - n["sifts_start_label_seq_id"])
    return n


def fig1(n, fdir):
    off_u = n["sifts_unp_start"] - n["sifts_start_label_seq_id"]   # UniProt = label + 34
    off_a = n["author_offset"]                                      # auth = label + 34
    labels = [1, 2, 3, 4, None, off_a + 1, off_a + 2, off_a + 3]   # label 35..37 = UniProt 69..71
    hit = n["first_observed_auth_seq_id"]                           # faulty key (label space)
    xs = np.arange(len(labels)) * 1.0
    y = {"uni": 2.0, "lab": 1.0, "auth": 0.0}
    fig, ax = plt.subplots(figsize=(7.0, 2.2))
    ax.set_xlim(-2.15, xs[-1] + 4.6)
    ax.set_ylim(-0.75, 2.75)
    ax.axis("off")
    for key, txt in (("uni", "UniProt position"), ("lab", "label_seq_id\n(SIFTS, mmCIF)"),
                     ("auth", "auth_seq_id\n(author)")):
        ax.text(-0.6, y[key], txt, ha="right", va="center", fontsize=8, color=INK)

    def box(x, yy, txt, observed=True, edge=INK2, lw=0.7, fc=SURF, weight="normal"):
        ax.add_patch(FancyBboxPatch((x - 0.36, yy - 0.2), 0.72, 0.4, boxstyle="round,pad=0,rounding_size=0.06",
                                    fc=fc, ec=edge, lw=lw, ls="-" if observed else (0, (2, 1.5))))
        ax.text(x, yy, txt, ha="center", va="center", fontsize=8, color=INK if observed else MUTED,
                fontweight=weight)

    for x, lab in zip(xs, labels):
        if lab is None:
            for key in y:
                ax.text(x, y[key], "…", ha="center", va="center", fontsize=9, color=INK2)
            continue
        observed = lab != n["sifts_start_label_seq_id"]
        uni, auth = lab + off_u, lab + off_a
        is_first = lab == n["first_observed_label_seq_id"]
        is_hit = lab == hit
        box(x, y["uni"], str(uni), edge=BLUE if is_first else (ORANGE if is_hit else INK2),
            lw=1.3 if (is_first or is_hit) else 0.7)
        box(x, y["lab"], str(lab), observed=observed, edge=ORANGE if is_hit else INK2,
            lw=1.3 if is_hit else 0.7)
        if observed:
            box(x, y["auth"], str(auth), edge=INK if is_first else INK2, lw=1.3 if is_first else 0.7,
                weight="bold" if is_first else "normal")
        else:
            box(x, y["auth"], "null", observed=False)
    x0 = xs[labels.index(n["first_observed_label_seq_id"])]
    xh = xs[labels.index(hit)]
    xu = xs[labels.index(n["sifts_start_label_seq_id"])]
    ax.text(xu, -0.43, "not modelled", ha="center", va="top", fontsize=7, color=MUTED)
    # SIFTS segment start bracket
    ax.annotate("", xy=(xu, y["lab"] + 0.22), xytext=(xu, y["uni"] - 0.22),
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8, ls=(0, (2, 1.5))))
    # correct mapping: auth 36 -> label 2 -> UniProt 36 (straight up)
    for ya, yb in ((y["auth"] + 0.22, y["lab"] - 0.22), (y["lab"] + 0.22, y["uni"] - 0.22)):
        ax.add_patch(FancyArrowPatch((x0, ya), (x0, yb), arrowstyle="-|>", mutation_scale=8, color=BLUE,
                                     lw=1.4, shrinkA=0, shrinkB=0, zorder=3))
    # faulty: auth 36 used as a label key -> key 36 -> UniProt 70
    ax.add_patch(FancyArrowPatch((x0 + 0.3, y["auth"] + 0.22), (xh - 0.15, y["lab"] - 0.22),
                                 arrowstyle="-|>", mutation_scale=8,
                                 color=ORANGE, lw=1.4, ls=(0, (4, 2)), shrinkA=0, shrinkB=0, zorder=3))
    ax.add_patch(FancyArrowPatch((xh, y["lab"] + 0.22), (xh, y["uni"] - 0.22), arrowstyle="-|>",
                                 mutation_scale=8, color=ORANGE, lw=1.4, ls=(0, (4, 2)), shrinkA=0, shrinkB=0,
                                 zorder=3))
    ax.annotate("", xy=(xh, y["uni"] + 0.32), xytext=(x0, y["uni"] + 0.32),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.7, shrinkA=0, shrinkB=0))
    ax.text((x0 + xh) / 2, y["uni"] + 0.4, f"{n['register_shift']}-residue register shift",
            ha="center", va="bottom", fontsize=7.5, color=INK)
    tx = xs[-1] + 0.75
    ax.text(tx, 2.45, f"PDB {n['entry']} chain {n['chain']} ({n['accession']})", fontsize=8, color=INK,
            fontweight="bold", va="top")
    ax.text(tx, 2.05, f"SIFTS segment: label {n['sifts_start_label_seq_id']}–{n['sifts_end_label_seq_id']}"
            f" ↔ UniProt {n['sifts_unp_start']}–{n['sifts_unp_end']};\nstart author number is null",
            fontsize=7.5, color=INK2, va="top", linespacing=1.25)
    ax.plot([tx, tx + 0.45], [1.15, 1.15], color=BLUE, lw=1.4)
    ax.text(tx + 0.6, 1.15, f"correct: atom auth {n['first_observed_auth_seq_id']} = label "
            f"{n['first_observed_label_seq_id']}\n→ UniProt {n['first_observed_correct_uniprot']}",
            fontsize=7.5, color=INK, va="center", linespacing=1.25)
    ax.plot([tx, tx + 0.45], [0.35, 0.35], color=ORANGE, lw=1.4, ls=(0, (4, 2)))
    ax.text(tx + 0.6, 0.35, f"faulty fallback: map keyed by label,\natom looked up by auth "
            f"{n['first_observed_auth_seq_id']} → key {hit}\n→ UniProt {n['faulty_uniprot_for_first_atom']}",
            fontsize=7.5, color=INK, va="center", linespacing=1.25)
    save(fig, fdir, "f1_mechanism")


# --------------------------------------------------------------------------- F2
def fig2(root, fdir):
    orig, corr = ledger(root, "e422"), ledger(root, "e427")
    mo, mc = seg_medians(orig), seg_medians(corr)
    both = sorted(set(mo) & set(mc))
    only_o, only_c = sorted(set(mo) - set(mc)), sorted(set(mc) - set(mo))
    x = np.array([mo[k] for k in both])
    yv = np.array([mc[k] for k in both])
    big = np.abs(yv - x) > DELTA
    changed = np.abs(yv - x) > 0
    ro = np.array([r["lddt"] for r in orig if finite(r.get("lddt"))], float)
    rc = np.array([r["lddt"] for r in corr if finite(r.get("lddt"))], float)
    nums = {
        "segment_key": "(entry, chain, accession); median residue lDDT per key; batch 1-4 cumulative ledgers",
        "segments_in_both": len(both), "segments_original_only": len(only_o),
        "segments_corrected_only": len(only_c),
        "segments_changed_any": int(changed.sum()),
        "segments_abs_delta_gt_0.2": int(big.sum()),
        "accessions_abs_delta_gt_0.2": len({both[i][2] for i in np.flatnonzero(big)}),
        "accessions_in_both": len({k[2] for k in both}),
        "residues_original_finite_lddt": int(ro.size), "residues_corrected_finite_lddt": int(rc.size),
        "rows_original_total": len(orig), "rows_corrected_total": len(corr),
        "share_lddt_lt_0.60_original": float((ro < ERR).mean()),
        "share_lddt_lt_0.60_corrected": float((rc < ERR).mean()),
    }
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.32})
    # (a) scatter
    a.plot([0, 1], [0, 1], color=MUTED, lw=0.7, ls=(0, (3, 2)), zorder=1)
    a.scatter(x[~big], yv[~big], s=9, color=BLUE, alpha=0.55, lw=0, zorder=2,
              label=f"|Δ| ≤ {DELTA} (n = {int((~big).sum())})")
    a.scatter(x[big], yv[big], s=9, color=ORANGE, alpha=0.8, lw=0, zorder=3,
              label=f"|Δ| > {DELTA} (n = {int(big.sum())})")
    if only_c:   # present only in the corrected ledger: drawn in a strip left of x = 0
        ys = [mc[k] for k in only_c]
        a.scatter([-0.1] * len(ys), ys, s=11, marker="D", facecolor="none", edgecolor=INK2, lw=0.7, zorder=3,
                  label=f"only in corrected (n = {len(only_c)})")
        a.axvline(-0.05, color=GRID, lw=0.8)
    if only_o:
        xs_ = [mo[k] for k in only_o]
        a.scatter(xs_, [-0.1] * len(xs_), s=11, marker="s", facecolor="none", edgecolor=INK2, lw=0.7, zorder=3,
                  label=f"only in original (n = {len(only_o)})")
        a.axhline(-0.05, color=GRID, lw=0.8)
    a.set_xlim(-0.15 if only_c else -0.02, 1.02)
    a.set_ylim(-0.15 if only_o else -0.02, 1.02)
    a.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    a.set_aspect("equal")
    a.set_xlabel("median lDDT, original mapping")
    a.set_ylabel("median lDDT, corrected mapping")
    a.legend(loc="lower right", frameon=False, handletextpad=0.2, borderaxespad=0.1, fontsize=7)
    a.text(0.0, 0.81, f"{int(big.sum())} of {len(both)} chain segments\nfrom {nums['accessions_abs_delta_gt_0.2']}"
           f" proteins\nmove by > {DELTA}", ha="left", va="top", fontsize=7.5, color=INK, linespacing=1.25)
    a.text(-0.3, 1.04, "a", transform=a.transAxes, fontsize=9, fontweight="bold", va="bottom")
    # (b) residue histograms (fraction of residues per bin)
    bins = np.linspace(0, 1, 41)
    b.hist(ro, bins=bins, weights=np.full(ro.size, 1 / ro.size), histtype="step", color=ORANGE, lw=1.3,
           ls=(0, (4, 2)), label=f"original ({ro.size:,} residues)")
    b.hist(rc, bins=bins, weights=np.full(rc.size, 1 / rc.size), histtype="step", color=BLUE, lw=1.3,
           label=f"corrected ({rc.size:,} residues)")
    b.axvline(ERR, color=INK2, lw=0.7, ls=(0, (3, 2)))
    b.set_yscale("log")
    b.set_ylim(1e-5, 2.0)
    b.text(ERR + 0.02, 0.97, f"lDDT < {ERR:.2f}\noriginal {100 * nums['share_lddt_lt_0.60_original']:.1f}%\n"
           f"corrected {100 * nums['share_lddt_lt_0.60_corrected']:.1f}%", ha="left", va="top", fontsize=7.5,
           color=INK, linespacing=1.3, transform=b.get_xaxis_transform())
    b.set_xlim(0, 1)
    b.set_xlabel("residue lDDT (Cα)")
    b.set_ylabel("fraction of residues (log scale)")
    b.yaxis.grid(True, color=GRID, lw=0.5)
    b.set_axisbelow(True)
    from matplotlib.lines import Line2D
    b.legend(handles=[Line2D([], [], color=ORANGE, lw=1.3, ls=(0, (4, 2))), Line2D([], [], color=BLUE, lw=1.3)],
             labels=[f"original\n({ro.size:,} residues)", f"corrected\n({rc.size:,} residues)"],
             loc="upper left", frameon=False, fontsize=7, handlelength=2.2, labelspacing=0.8)
    b.text(-0.18, 1.04, "b", transform=b.transAxes, fontsize=9, fontweight="bold", va="bottom")
    save(fig, fdir, "f2_scope")
    return nums


# --------------------------------------------------------------------------- F3
TRAIN = [("SAME_ACCESSION", "same accession"), ("SEQ95", "≥95% identity"), ("SEQ30", "≥30% identity"),
         ("NOVEL", "no relative ≥30%")]
TEMPL = [("SAME_ACC", "same-accession template"), ("OTHER", "other template"), ("NONE", "no template"),
         ("UNRESOLVED", "unresolved")]


def template_counts(vdir):
    """Template-class counts from V5 (out/v5_templates.json), or None if V5 has not been run."""
    p = vdir / "v5_templates.json"
    if not p.exists():
        return None
    return dict(json.loads(p.read_text())["results"]["template_class_counts"])


def fig3(vdir, fdir):
    v2 = json.loads((vdir / "v2_prior_exposure.json").read_text())
    recs = [r for r in v2["accessions"].values() if r["in_e422"]]
    train = {c: sum(r["training"]["class"] == c for r in recs) for c, _ in TRAIN}
    unres = len(recs) - sum(train.values())
    rows = [("training-era PDB\n(released ≤ 2018-04-30)", TRAIN, train, len(recs))]
    tc = template_counts(vdir)
    if tc is not None:
        rows.append(("AFDB model\ntemplates", TEMPL, {c: tc.get(c, 0) for c, _ in TEMPL}, sum(tc.values())))
    fig, axes = plt.subplots(len(rows), 1, figsize=(7.0, 0.62 + 0.55 * len(rows)), squeeze=False)
    for ax, (name, classes, counts, n) in zip(axes[:, 0], rows):
        left = 0.0
        shown = [(c, lab) for c, lab in classes if counts.get(c, 0) or c != "UNRESOLVED"]
        for j, (c, lab) in enumerate(shown):
            k = counts.get(c, 0)
            w = 100 * k / n
            ax.barh(0, w, left=left, color=SEQ[j], edgecolor=SURF, lw=1.5, height=0.62,
                    label=f"{lab} ({k})", hatch="////" if c == "UNRESOLVED" else None)
            if w >= 5:
                ax.text(left + w / 2, 0, str(k), ha="center", va="center", fontsize=7.5,
                        color="#ffffff" if j < 2 else INK)
            left += w
        ax.set_xlim(0, 100)
        ax.set_ylim(-0.45, 0.45)
        ax.set_yticks([0], [f"{name}\nn = {n}"])
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.legend(ncol=len(shown), frameon=False, loc="lower left", bbox_to_anchor=(0, 0.92),
                  handlelength=1.0, handletextpad=0.4, columnspacing=1.2, fontsize=7)
        ax.set_xticks([0, 25, 50, 75, 100])
        if ax is not axes[-1, 0]:
            ax.set_xticklabels([])
    axes[-1, 0].set_xlabel("share of proteins (%)")
    fig.subplots_adjust(hspace=1.1)
    save(fig, fdir, "f3_exposure")
    return {"n_proteins": len(recs), "training_class_counts": train, "training_unresolved": unres,
            "template_class_counts": tc,
            "template_panel": "drawn from out/v5_templates.json" if tc is not None
            else "STUB: out/v5_templates.json absent at figure time; template bar not drawn"}


def main(root):
    root = Path(root)
    vdir = root / "paper/arxiv/verification/out"
    fdir = root / "paper/arxiv/figures/rev3"
    fdir.mkdir(parents=True, exist_ok=True)
    n1 = f1_numbers(root, vdir)
    fig1(n1, fdir)
    n2 = fig2(root, fdir)
    n3 = fig3(vdir, fdir)
    out = {"schema": "figures-rev3-numbers-v1",
           "note": "Every number annotated on the rev3 figures. pLDDT is dropped on load from e422/e427 "
                   "ledgers and never used (registered blinding of P6/P7).",
           "f1_mechanism_9qj6": n1, "f2_scope": n2, "f3_exposure": n3}
    (vdir / "figures_rev3_numbers.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
