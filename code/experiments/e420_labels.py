"""e420 gate (g): label construction from the frozen channels (offline core).

Registration: docs/e420_registration.md sections 5-7, 10(g), 12.
Frozen rules: docs/g1_successor_label_rules.json (schema
g1-successor-label-rules-v1) -- thresholds, atom selection, alignment, and
missingness rules were frozen there BEFORE any threshold interpretation; this
module only implements them and never tunes them. Pool/floor helpers are
reused from experiments/e420_pool.py (evaluate_support_floors, pool_verdict,
base_accession, canonical_digest); SourceRole and the CIF tokenizer are
reused from experiments/e420_grammar_validate.py. Typed failures are frozen
dataclasses with codes in the registration section 7 style (GrammarFailure
precedent), raised via E420LabelError -- never silent, never repaired.

OFFLINE, PURE, DETERMINISTIC. No network, no source access, no prints, no
timestamps, no real coordinates (synthetic fixtures live in the tests).
numpy + stdlib only, plus the two sibling e420 modules. Same inputs ->
identical outputs.

Channels (registration 5):

E -- coordinate error (5.1): parse entry-files mmCIF CA rows (Atom_site
loop; label_atom_id == CA; label_comp_id in the standard amino-acid set;
pdbx_PDB_model_num == 1); registered global rigid alignment as a binary64
row-vector Kabsch fit with math.fsum centroids, the pinned backend call
numpy.linalg.svd(H, full_matrices=False, compute_uv=True, hermitian=False),
the deterministic reflection sign rule R = U diag(1,1,sign(det(UV.T))) V.T,
and NO per-residue refit. Distances d_i = ||mobile_i - ref_i|| AFTER
superposition. Complete four-threshold C-alpha lDDT: neighbor set
N_i = {j != i : reference squared CA distance strictly < 225.0 A^2};
thresholds t in {0.5, 1.0, 2.0, 4.0} A; a pair is preserved iff
|d_ij^mobile - d_ij^reference| < t (strict); lDDT_i = mean of the four
fractions (math.fsum). Residues with < 1 neighbor, nonfinite evidence, or
an incomplete/nonfinite neighbor are typed U. E=correct iff lddt_i >= 0.60
AND d_i <= 4.0; E=error iff lddt_i < 0.60 AND d_i > 4.0 (complementary
strict inequalities); mixed (threshold-boundary conflict), nonfinite, or
rank-degenerate evidence yields E=U with a retained typed reason. pLDDT,
B-factors, absorption proxies, sequence-only scores, and any AF-derived
feature are barred from E and are never read anywhere in this module.

F -- dynamics/flexibility (5.2): the mean of EXACTLY three finite,
nonnegative ATLAS replicate values, computed only as math.fsum(values)/3;
flexible iff mean >= 1.5 A; ordered iff mean <= 1.0 A; U for the open
interval (1.0, 1.5), a missing member/value, a nonfinite or negative
replicate, a nonzero flag, or an ambiguous correspondence. A two-replicate
or partial mean is forbidden.

Combined (5.3): (correct,ordered)=OC, (correct,flexible)=A,
(error,ordered)=B, (error,flexible)=O; any U endpoint yields U with a
retained reason. Callers cannot provide a state, class count, residue
count, or favorable digest -- states are recomputed here from evidence.

Gap policy (5.4): the ATLAS correspondence is rejected as a whole only for
malformed rows, a PDB number outside the chain range, ZERO joined rows, or
5 * unmapped_count > total_count (integer arithmetic; strictly > 20%).
Boundary equality (5 * unmapped == total) is ACCEPTED. SIFTS-absent rows
stay per-row U with full no-drop accounting.

Census: four_cell_census aggregates per-residue label rows into per-cell
unique (base_accession, uniprot_num) identity counts, per-protein and
per-class counts, per-fold support (folds from the frozen group map via
pool.evaluate_support_floors / pool.assign_folds), the registered floor
evaluation, and the verdict from pool.pool_verdict.
"""
from __future__ import annotations

import math
import numbers
from dataclasses import dataclass
from enum import Enum
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

import e420_grammar_validate as gv
import e420_pool as pool

__all__ = [
    "LABELS_SCHEMA", "CENSUS_SCHEMA", "EMITTED_BY",
    "STANDARD_AA3", "NEIGHBOR_CUTOFF_A2", "LDDT_THRESHOLDS_A",
    "LDDT_CORRECT_MIN", "DIST_CORRECT_MAX_A", "INCLUSION_RADIUS_A",
    "RMSF_FLEXIBLE_MIN_A", "RMSF_ORDERED_MAX_A", "RMSF_REPLICATE_COUNT",
    "GAP_REJECT_FACTOR", "CORRESP_FLAG_NAMES",
    "E_CORRECT", "E_ERROR", "F_FLEXIBLE", "F_ORDERED", "STATE_U",
    "CELL_OC", "CELL_A", "CELL_B", "CELL_O", "CELL_U", "ALL_CELLS",
    "CELL_OF_STATE", "LabelFailureCode", "LabelFailure", "E420LabelError",
    "parse_mmcif_ca", "ca_coords_array", "kabsch_fit", "lddt_per_residue",
    "e_residue_state", "f_residue_state", "combine_residue_labels",
    "gap_share_status", "four_cell_census",
]

# --------------------------------------------------------------------------
# Schema identity and frozen thresholds (docs/g1_successor_label_rules.json).
# --------------------------------------------------------------------------
LABELS_SCHEMA = "e420-labels-v1"
CENSUS_SCHEMA = "e420-four-cell-census-v1"
EMITTED_BY = "e420-labels-v1"

STANDARD_AA3 = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
})

NEIGHBOR_CUTOFF_A2 = 225.0        # strict < on reference squared CA distance
INCLUSION_RADIUS_A = 15.0         # sqrt of the frozen squared cutoff
LDDT_THRESHOLDS_A = (0.5, 1.0, 2.0, 4.0)
LDDT_CORRECT_MIN = 0.60
DIST_CORRECT_MAX_A = 4.0

RMSF_FLEXIBLE_MIN_A = 1.5
RMSF_ORDERED_MAX_A = 1.0
RMSF_REPLICATE_COUNT = 3
GAP_REJECT_FACTOR = 5             # rejected iff 5 * unmapped > total (strict)
CORRESP_FLAG_NAMES = ("gap", "hetatm", "corrected", "no_fullbb", "no_ca")

# Registered state / cell literals (registration 5.1-5.3).
E_CORRECT = "correct"
E_ERROR = "error"
F_FLEXIBLE = "flexible"
F_ORDERED = "ordered"
STATE_U = "U"
CELL_OC = "OC"
CELL_A = "A"
CELL_B = "B"
CELL_O = "O"
CELL_U = "U"
ALL_CELLS = (CELL_OC, CELL_A, CELL_B, CELL_O, CELL_U)
CELL_OF_STATE = {
    (E_CORRECT, F_ORDERED): CELL_OC,
    (E_CORRECT, F_FLEXIBLE): CELL_A,
    (E_ERROR, F_ORDERED): CELL_B,
    (E_ERROR, F_FLEXIBLE): CELL_O,
}

# Memory guard for the O(N^2) lDDT distance matrices (AGENTS section 5).
LDDT_MAX_RESIDUES = 20_000

# Required Atom_site tags for the registered CA selection (lowercased).
_MMCIF_REQUIRED_TAGS = (
    "label_atom_id",
    "label_comp_id",
    "auth_seq_id",
    "cartn_x",
    "cartn_y",
    "cartn_z",
    "pdbx_pdb_model_num",
)


# --------------------------------------------------------------------------
# Typed failures (registration section 7 style: coded, frozen, retained).
# --------------------------------------------------------------------------
class LabelFailureCode(str, Enum):
    """Every typed label-construction failure class this module can emit."""

    MMCIF_RAW_BYTES_NOT_BYTES = "mmcif_raw_bytes_not_bytes"
    MMCIF_RAW_BYTES_EMPTY = "mmcif_raw_bytes_empty"
    MMCIF_DECODE_FAILED = "mmcif_decode_failed"
    MMCIF_TOKEN_INVALID = "mmcif_token_invalid"
    MMCIF_ATOM_SITE_LOOP_MISSING = "mmcif_atom_site_loop_missing"
    MMCIF_LOOP_MALFORMED = "mmcif_loop_malformed"
    MMCIF_TAG_MISSING = "mmcif_tag_missing"
    MMCIF_COORD_INVALID = "mmcif_coord_invalid"
    MMCIF_RESIDUE_NUM_INVALID = "mmcif_residue_num_invalid"
    MMCIF_NO_CA = "mmcif_no_ca"
    KABSCH_INPUT_INVALID = "kabsch_input_invalid"
    LDDT_INPUT_INVALID = "lddt_input_invalid"
    STATE_INPUT_INVALID = "state_input_invalid"
    GAP_SHARE_INPUT_INVALID = "gap_share_input_invalid"
    CENSUS_INPUT_INVALID = "census_input_invalid"


@dataclass(frozen=True)
class LabelFailure:
    """Typed failure record (dataclass, never a silent fallback)."""

    code: LabelFailureCode
    role: Any = None                  # gv.SourceRole when channel-bound
    detail: str = ""


class E420LabelError(ValueError):
    """Typed label-construction error carrying its frozen LabelFailure."""

    def __init__(self, code: LabelFailureCode, detail: str = "", role: Any = None) -> None:
        self.failure = LabelFailure(code=code, role=role, detail=detail)
        self.code = code
        self.detail = detail
        super().__init__(f"{code.value}: {detail}")


def _fail(code: LabelFailureCode, detail: str = "", role: Any = None) -> Any:
    raise E420LabelError(code, detail, role=role)


# --------------------------------------------------------------------------
# Input helpers.
# --------------------------------------------------------------------------
def _require_points(name: str, value: Any, code: LabelFailureCode, *, finite: bool = True) -> np.ndarray:
    """Coerce to a nonempty float64 Nx3 array (typed failure otherwise)."""
    try:
        arr = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError):
        _fail(code, f"{name}: not convertible to a float64 array")
    if arr.ndim != 2 or arr.shape[1] != 3:
        _fail(code, f"{name}: shape {arr.shape} is not Nx3")
    if arr.shape[0] == 0:
        _fail(code, f"{name}: zero residues")
    if finite and not bool(np.isfinite(arr).all()):
        bad = int((~np.isfinite(arr)).sum())
        _fail(code, f"{name}: {bad} nonfinite entries")
    return arr


def ca_coords_array(rows: Sequence[Any]) -> np.ndarray:
    """[(resnum, x, y, z), ...] from parse_mmcif_ca -> Nx3 coordinate array."""
    try:
        arr = np.asarray([(float(r[1]), float(r[2]), float(r[3])) for r in rows],
                         dtype=np.float64)
    except (TypeError, ValueError, IndexError):
        _fail(LabelFailureCode.LDDT_INPUT_INVALID, "CA rows not (resnum, x, y, z) tuples")
    if arr.ndim != 2 or arr.shape[1] != 3 or arr.shape[0] == 0:
        _fail(LabelFailureCode.LDDT_INPUT_INVALID, f"CA rows decode to shape {arr.shape}")
    return arr


# --------------------------------------------------------------------------
# 1. Entry-files mmCIF CA parsing (E-channel reference / predictor input).
# --------------------------------------------------------------------------
def _extract_cif_loops(tokens: Sequence[str]) -> list[tuple[list[str], list[str]]]:
    """loops = [(tag_names, values), ...] in file order.

    Mirrors the registered tokenizer's grammar (gv._tokenize_cif): loop_
    headers are the run of '_' tokens after loop_; values run until the next
    loop_, stop_, data_ block, or a new '_' tag. A value token beginning with
    '_' truncates the loop, which surfaces as an arity mismatch (typed), the
    same ambiguity the grammar module documents for legacy files.
    """
    loops: list[tuple[list[str], list[str]]] = []
    i, n = 0, len(tokens)
    while i < n:
        tok = tokens[i]
        if tok == "loop_":
            i += 1
            tags: list[str] = []
            while i < n and tokens[i].startswith("_"):
                tags.append(tokens[i])
                i += 1
            values: list[str] = []
            while (i < n and tokens[i] != "loop_" and tokens[i] != "stop_"
                   and not tokens[i].startswith("_") and not tokens[i].startswith("data_")):
                values.append(tokens[i])
                i += 1
            loops.append((tags, values))
        elif tok.startswith("_"):
            i += 2  # scalar item name + its single value token
        else:       # data_ block header, stop_, or stray value token
            i += 1
    return loops


def parse_mmcif_ca(raw: Any) -> list[tuple[int, float, float, float]]:
    """CA rows of polymer (standard-amino-acid) residues from an entry-files
    mmCIF, in file order: [(auth residue number, x, y, z), ...].

    Selection (frozen): Atom_site loop rows with label_atom_id == 'CA',
    label_comp_id in the standard AA set, pdbx_PDB_model_num == '1'.
    residue_number is the PDB AUTHOR residue number (auth_seq_id), the
    registered mapping join key (registration 4.4). Chain scoping and the
    one-to-one mapping guards happen upstream on registered paths; this
    function never imputes, renumbers, or drops a selected row silently.

    Typed failures (LabelFailureCode): raw bytes class, decode, tokenizer,
    missing/multiple Atom_site loops, malformed loop (arity), missing
    required tags, unparseable/nonfinite coordinates or residue numbers on a
    selected row, and zero surviving CA rows.
    """
    role = gv.SourceRole.REFERENCE_COORDINATE
    if not isinstance(raw, (bytes, bytearray)):
        _fail(LabelFailureCode.MMCIF_RAW_BYTES_NOT_BYTES, f"type {type(raw).__name__}", role)
    data = bytes(raw)
    if not data:
        _fail(LabelFailureCode.MMCIF_RAW_BYTES_EMPTY, "empty input", role)
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        _fail(LabelFailureCode.MMCIF_DECODE_FAILED, "not utf-8", role)
    try:
        tokens = gv._tokenize_cif(text)
    except gv._GrammarError as exc:
        _fail(LabelFailureCode.MMCIF_TOKEN_INVALID, exc.detail, role)
    loops = _extract_cif_loops(tokens)
    atom_loops = [(tags, values) for tags, values in loops
                  if any(tag.lower().startswith("_atom_site.") for tag in tags)]
    if not atom_loops:
        _fail(LabelFailureCode.MMCIF_ATOM_SITE_LOOP_MISSING, "no _atom_site loop", role)
    if len(atom_loops) > 1:
        _fail(LabelFailureCode.MMCIF_LOOP_MALFORMED,
              f"multiple _atom_site loops ({len(atom_loops)})", role)
    tags, values = atom_loops[0]
    if len(tags) == 0:
        _fail(LabelFailureCode.MMCIF_LOOP_MALFORMED, "empty loop header", role)
    lower_tags = []
    for tag in tags:
        lowered = tag.lower()
        # category.attr -> attribute name (tags are case-insensitive per CIF).
        lower_tags.append(lowered.split(".", 1)[1] if "." in lowered else lowered)
    missing = sorted(set(_MMCIF_REQUIRED_TAGS) - set(lower_tags))
    if missing:
        _fail(LabelFailureCode.MMCIF_TAG_MISSING, f"missing {missing}", role)
    if len(values) % len(tags) != 0:
        _fail(LabelFailureCode.MMCIF_LOOP_MALFORMED,
              f"loop arity {len(values)} % {len(tags)} != 0", role)
    col = {name: pos for pos, name in enumerate(lower_tags)}
    out: list[tuple[int, float, float, float]] = []
    total_rows = len(values) // len(tags)
    for row_index in range(total_rows):
        row = values[row_index * len(tags):(row_index + 1) * len(tags)]
        if row[col["pdbx_pdb_model_num"]] != "1":
            continue
        if row[col["label_atom_id"]] != "CA":
            continue
        if row[col["label_comp_id"]].upper() not in STANDARD_AA3:
            continue
        seq_token = row[col["auth_seq_id"]]
        try:
            resnum = int(seq_token)
        except ValueError:
            _fail(LabelFailureCode.MMCIF_RESIDUE_NUM_INVALID,
                  f"row {row_index}: auth_seq_id={seq_token!r}", role)
        xyz: list[float] = []
        for tag in ("cartn_x", "cartn_y", "cartn_z"):
            token = row[col[tag]]
            try:
                val = float(token)
            except ValueError:
                _fail(LabelFailureCode.MMCIF_COORD_INVALID,
                      f"row {row_index}: {tag}={token!r}", role)
            if not math.isfinite(val):
                _fail(LabelFailureCode.MMCIF_COORD_INVALID,
                      f"row {row_index}: {tag}={token!r} nonfinite", role)
            xyz.append(val)
        out.append((resnum, xyz[0], xyz[1], xyz[2]))
    if not out:
        _fail(LabelFailureCode.MMCIF_NO_CA,
              f"0 of {total_rows} Atom_site rows are model-1 standard-AA CA", role)
    return out


# --------------------------------------------------------------------------
# 2. Registered Kabsch fit (E-channel alignment; registration 5.1).
# --------------------------------------------------------------------------
def kabsch_fit(mobile: Any, ref: Any) -> dict[str, Any]:
    """Registered global rigid alignment (one fit for the whole chain; no
    per-residue refit).

    Row-vector convention: H = (mobile - mu_mobile)^T (ref - mu_ref);
    centroids via math.fsum per column; SVD via the pinned backend call
    numpy.linalg.svd(H, full_matrices=False, compute_uv=True,
    hermitian=False); deterministic reflection sign rule
    R = U diag(1, 1, sign(det(UV.T))) V.T (R is always a proper rotation,
    det + 1). superposed_i = (mobile_i - mu_mobile) R + mu_ref and
    d_i = ||mobile_i - ref_i|| AFTER superposition (i.e. against the
    reference). Rank info is returned (rank_ok / singular_values) so the
    caller can map tie-degenerate fits to E=U per the frozen rule.
    """
    code = LabelFailureCode.KABSCH_INPUT_INVALID
    P = _require_points("mobile", mobile, code)
    Q = _require_points("ref", ref, code)
    if P.shape != Q.shape:
        _fail(code, f"shape mismatch mobile={P.shape} ref={Q.shape}")
    n = P.shape[0]
    mu_mobile = np.array([math.fsum(float(v) for v in P[:, c]) / n for c in range(3)],
                         dtype=np.float64)
    mu_ref = np.array([math.fsum(float(v) for v in Q[:, c]) / n for c in range(3)],
                      dtype=np.float64)
    P_centered = P - mu_mobile
    Q_centered = Q - mu_ref
    H = P_centered.T @ Q_centered
    U, S, Vt = np.linalg.svd(H, full_matrices=False, compute_uv=True, hermitian=False)
    uv_t = U @ Vt
    det_sign = 1.0 if float(np.linalg.det(uv_t)) >= 0.0 else -1.0
    R = U @ np.diag([1.0, 1.0, det_sign]) @ Vt
    superposed = P_centered @ R + mu_ref
    distances = np.linalg.norm(superposed - Q, axis=1)
    singular = np.asarray(S, dtype=np.float64)
    rank_ok = bool(singular.size == 3 and singular[-1] > 0.0)
    return {
        "rotation": R,
        "superposed": superposed,
        "distances": distances,
        "centroid_mobile": mu_mobile,
        "centroid_ref": mu_ref,
        "singular_values": singular,
        "determinant_sign": det_sign,
        "rank_ok": rank_ok,
        "n_residues": int(n),
        "residual_max": float(np.max(distances)) if n else 0.0,
    }


# --------------------------------------------------------------------------
# 3. Complete four-threshold C-alpha lDDT (E-channel residue metric).
# --------------------------------------------------------------------------
def lddt_per_residue(
    superposed_mobile: Any,
    reference: Any,
    *,
    inclusion_radius_a: float = INCLUSION_RADIUS_A,
    thresholds_a: Sequence[float] = LDDT_THRESHOLDS_A,
) -> dict[str, Any]:
    """Per-index lDDT-C-alpha of a superposed mobile chain against the
    reference (both Nx3, same order; no re-matching here).

    Frozen rule: N_i = {j != i : ||ref_i - ref_j||^2 strictly < cutoff^2}
    (default 15.0 A -> 225.0 A^2); for each threshold t in thresholds_a
    (default (0.5, 1.0, 2.0, 4.0) A) the preserved fraction of N_i is
    #{j in N_i : |d_ij^mobile - d_ij^reference| < t} / |N_i| (STRICT <);
    lDDT_i = math.fsum of the four fractions / 4. Residues with < 1
    neighbor are typed U ("u_insufficient_neighbors"); nonfinite own rows
    ("u_nonfinite") and finite-reference neighbors whose mobile row is
    nonfinite ("u_nonfinite_neighbor") are typed U -- incomplete-neighbor
    evidence is never silently included or dropped. lDDT is np.nan exactly
    at U positions. Pure vectorized numpy; deterministic; O(N^2) memory
    guarded by LDDT_MAX_RESIDUES.
    """
    code = LabelFailureCode.LDDT_INPUT_INVALID
    M = _require_points("superposed_mobile", superposed_mobile, code, finite=False)
    R = _require_points("reference", reference, code, finite=False)
    if M.shape != R.shape:
        _fail(code, f"shape mismatch superposed_mobile={M.shape} reference={R.shape}")
    n = M.shape[0]
    if n > LDDT_MAX_RESIDUES:
        _fail(code, f"n={n} exceeds the {LDDT_MAX_RESIDUES}-residue memory guard")
    cutoff = float(inclusion_radius_a)
    if not math.isfinite(cutoff) or cutoff <= 0.0:
        _fail(code, f"inclusion_radius_a={inclusion_radius_a!r}")
    cutoff2 = cutoff * cutoff
    thr = tuple(float(t) for t in thresholds_a)
    if not thr or any(not math.isfinite(t) or t <= 0.0 for t in thr):
        _fail(code, f"thresholds_a={thresholds_a!r}")

    finite_m = np.isfinite(M).all(axis=1)
    finite_r = np.isfinite(R).all(axis=1)
    with np.errstate(invalid="ignore"):
        # Neighbor rule on SQUARED reference distances (strict < 225.0 A^2);
        # preservation on distances, strictly < threshold (frozen rule).
        D2_ref = ((R[:, None, :] - R[None, :, :]) ** 2).sum(axis=-1)
        D2_mob = ((M[:, None, :] - M[None, :, :]) ** 2).sum(axis=-1)
        D_ref = np.sqrt(D2_ref)
        D_mob = np.sqrt(D2_mob)
        diffs = np.abs(D_mob - D_ref)
        preserved = {t: diffs < t for t in thr}  # NaN diffs compare False
    geo_neighbors = ((D2_ref < cutoff2) & ~np.eye(n, dtype=bool)
                     & finite_r[:, None] & finite_r[None, :])

    lddt = np.full(n, np.nan, dtype=np.float64)
    statuses: list[str] = []
    neighbor_counts: list[int] = []
    for i in range(n):
        if not (bool(finite_m[i]) and bool(finite_r[i])):
            statuses.append("u_nonfinite")
            neighbor_counts.append(0)
            continue
        nb = geo_neighbors[i]
        count = int(nb.sum())
        if count == 0:
            statuses.append("u_insufficient_neighbors")
            neighbor_counts.append(0)
            continue
        if bool((~finite_m[nb]).any()):
            statuses.append("u_nonfinite_neighbor")
            neighbor_counts.append(count)
            continue
        fractions = [float(preserved[t][i][nb].mean()) for t in thr]
        lddt[i] = math.fsum(fractions) / len(fractions)
        statuses.append("ok")
        neighbor_counts.append(count)
    return {
        "lddt": lddt,
        "statuses": statuses,
        "neighbor_counts": neighbor_counts,
        "n_residues": int(n),
        "thresholds_a": thr,
        "inclusion_radius_a": cutoff,
        "neighbor_cutoff_a2": cutoff2,
    }


# --------------------------------------------------------------------------
# 4. Per-residue E / F states and combined four-cell labels (5.1-5.3).
# --------------------------------------------------------------------------
def e_residue_state(lddt: Any, distance: Any) -> tuple[str, str]:
    """Frozen E decision per residue (registration 5.1).

    correct iff lddt >= 0.60 AND d <= 4.0; error iff lddt < 0.60 AND
    d > 4.0 (complementary strict inequalities); anything else -- including
    the mixed threshold-boundary conflict and nonfinite evidence -- is
    (U, typed_reason). Rank/tie-degenerate fits are mapped to U by the
    caller via kabsch_fit()["rank_ok"] before this function is used for a
    claim-bearing row.
    """
    for name, value in (("lddt", lddt), ("distance", distance)):
        if isinstance(value, bool) or not isinstance(value, numbers.Real):
            _fail(LabelFailureCode.STATE_INPUT_INVALID, f"{name}={value!r} is not a real number")
    lddt_f = float(lddt)
    distance_f = float(distance)
    if not math.isfinite(lddt_f):
        return STATE_U, "u_nonfinite_lddt"
    if not math.isfinite(distance_f):
        return STATE_U, "u_nonfinite_distance"
    if lddt_f >= LDDT_CORRECT_MIN and distance_f <= DIST_CORRECT_MAX_A:
        return E_CORRECT, ""
    if lddt_f < LDDT_CORRECT_MIN and distance_f > DIST_CORRECT_MAX_A:
        return E_ERROR, ""
    return STATE_U, "u_e_conflict"


def f_residue_state(values: Any, flags: Any = None) -> tuple[str, str]:
    """Frozen F decision per residue (registration 5.2).

    values: EXACTLY RMSF_REPLICATE_COUNT (=3) replicate values; the mean is
    computed ONLY as math.fsum(values) / 3. flexible iff mean >= 1.5 A;
    ordered iff mean <= 1.0 A; U for the open interval (1.0, 1.5), a
    non-numeric/nonfinite/negative replicate, or a wrong replicate count (a
    two-replicate or partial mean is forbidden). flags: optional mapping of
    the five corresp flags; any value that is not the literal 0/1, or any
    nonzero flag, yields U with a typed reason.
    """
    if not isinstance(values, (list, tuple)) or len(values) != RMSF_REPLICATE_COUNT:
        return STATE_U, "u_replicates_not_three"
    parsed: list[float] = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, numbers.Real):
            return STATE_U, "u_replicate_type_invalid"
        number = float(value)
        if not math.isfinite(number):
            return STATE_U, "u_replicate_nonfinite"
        if number < 0.0:
            return STATE_U, "u_replicate_negative"
        parsed.append(number)
    if flags is not None:
        if not isinstance(flags, Mapping):
            return STATE_U, "u_flags_not_a_mapping"
        for name in CORRESP_FLAG_NAMES:
            if name not in flags:
                continue
            raw = flags[name]
            if isinstance(raw, bool) or not isinstance(raw, (int, str)):
                return STATE_U, f"u_flag_literal_invalid:{name}"
            literal = raw.strip() if isinstance(raw, str) else raw
            if literal not in ("0", "1", 0, 1):
                return STATE_U, f"u_flag_literal_invalid:{name}"
            if literal in ("1", 1):
                return STATE_U, f"u_nonzero_flag:{name}"
    mean = math.fsum(parsed) / RMSF_REPLICATE_COUNT
    if mean >= RMSF_FLEXIBLE_MIN_A:
        return F_FLEXIBLE, ""
    if mean <= RMSF_ORDERED_MAX_A:
        return F_ORDERED, ""
    return STATE_U, "u_rmsf_ambiguous_interval"


def _normalize_endpoint(item: Any, allowed: Sequence[str], code: LabelFailureCode,
                        side: str) -> tuple[str, str]:
    state, reason = item, ""
    if isinstance(item, (tuple, list)):
        if len(item) != 2 or not isinstance(item[1], str):
            _fail(code, f"{side} endpoint {item!r} is not (state, reason)")
        state, reason = item
    if not isinstance(state, str) or state not in allowed:
        _fail(code, f"{side} state {state!r} outside {sorted(allowed)}")
    return state, reason


def combine_residue_labels(e_states: Sequence[Any], f_states: Sequence[Any]) -> list[dict[str, str]]:
    """Registered combination (registration 5.3): per residue exactly one of
    OC / A / B / O / U, recomputed here from the two endpoint states.

    e_states / f_states: per-residue state strings ("correct"/"error"/"U"
    and "flexible"/"ordered"/"U") or (state, typed_reason) pairs; use the
    outputs of e_residue_state / f_residue_state. Any U endpoint yields
    ("U", combined reason retaining both endpoint reasons). Malformed
    inputs are typed failures; callers cannot inject a cell directly.
    Returns a list of {"cell", "e_state", "f_state", "e_reason", "f_reason",
    "reason"} dicts in input order.
    """
    code = LabelFailureCode.STATE_INPUT_INVALID
    if not isinstance(e_states, (list, tuple)) or not isinstance(f_states, (list, tuple)):
        _fail(code, "e_states / f_states must be sequences")
    if len(e_states) != len(f_states):
        _fail(code, f"length mismatch e={len(e_states)} f={len(f_states)}")
    out: list[dict[str, str]] = []
    for e_item, f_item in zip(e_states, f_states):
        e_state, e_reason = _normalize_endpoint(e_item, (E_CORRECT, E_ERROR, STATE_U), code, "e")
        f_state, f_reason = _normalize_endpoint(f_item, (F_FLEXIBLE, F_ORDERED, STATE_U), code, "f")
        if e_state == STATE_U or f_state == STATE_U:
            parts: list[str] = []
            if e_state == STATE_U:
                parts.append(f"u_e_endpoint:{e_reason or 'unspecified'}")
            if f_state == STATE_U:
                parts.append(f"u_f_endpoint:{f_reason or 'unspecified'}")
            cell, reason = CELL_U, ";".join(parts)
        else:
            cell, reason = CELL_OF_STATE[(e_state, f_state)], ""
        out.append({
            "cell": cell,
            "e_state": e_state,
            "f_state": f_state,
            "e_reason": e_reason,
            "f_reason": f_reason,
            "reason": reason,
        })
    return out


# --------------------------------------------------------------------------
# 5. Gap policy (registration 5.4): 5 * unmapped > total, strict.
# --------------------------------------------------------------------------
def gap_share_status(unmapped_count: Any, total_count: Any) -> dict[str, Any]:
    """Integer gap-share rule, boundary-equality accepting.

    The ATLAS correspondence is rejected as a whole iff ZERO joined rows or
    5 * unmapped_count > total_count (integer arithmetic; strictly > 20%).
    5 * unmapped == total is ACCEPTED. Inputs must be nonnegative integers
    (bool excluded); anything else is a typed failure.
    """
    code = LabelFailureCode.GAP_SHARE_INPUT_INVALID
    for name, value in (("unmapped_count", unmapped_count), ("total_count", total_count)):
        if isinstance(value, bool) or not isinstance(value, numbers.Integral):
            _fail(code, f"{name}={value!r} is not an integer")
        if int(value) < 0:
            _fail(code, f"{name}={value!r} is negative")
    unmapped = int(unmapped_count)
    total = int(total_count)
    if total == 0:
        return {"accepted": False, "detail": "zero_joined_rows",
                "unmapped_count": unmapped, "total_count": total}
    if GAP_REJECT_FACTOR * unmapped > total:
        return {"accepted": False,
                "detail": f"{GAP_REJECT_FACTOR}*unmapped={GAP_REJECT_FACTOR * unmapped} > total={total}",
                "unmapped_count": unmapped, "total_count": total}
    return {"accepted": True, "detail": "",
            "unmapped_count": unmapped, "total_count": total}


# --------------------------------------------------------------------------
# 6. Four-cell census (registration 2.1 / 6 floors; verdict via pool module).
# --------------------------------------------------------------------------
def four_cell_census(
    label_rows: Any,
    group_map: Any,
    *,
    cap_bound: bool = False,
    fold_of_group: Any = None,
) -> dict[str, Any]:
    """Aggregate per-residue label rows into the registered four-cell census.

    label_rows: mappings with exactly the required keys base_accession,
    uniprot_num (int), cell (one of OC/A/B/O/U); extra ledger keys are
    allowed. U rows are retained for accounting but never enter a cell.
    group_map: frozen pool map accession -> group_id. Folds come from the
    group map (pool.evaluate_support_floors applies the registered
    assign_folds rule unless fold_of_group is given). The floor evaluation
    and the verdict are the pool module's own functions
    (pool.evaluate_support_floors + pool.pool_verdict with cap_bound);
    per-cell counts here key unique (base_accession, uniprot_num)
    identities exactly as registered. Deterministic; carries a
    self-digest over the payload without that field (registration 12 rule).
    """
    code = LabelFailureCode.CENSUS_INPUT_INVALID
    if not isinstance(label_rows, (list, tuple)):
        _fail(code, "label_rows must be a sequence of row mappings")
    if not isinstance(group_map, Mapping):
        _fail(code, "group_map must be a mapping accession -> group_id")
    identity_sets: dict[str, set[tuple[str, int]]] = {c: set() for c in ALL_CELLS}
    row_counts: dict[str, int] = {c: 0 for c in ALL_CELLS}
    for pos, row in enumerate(label_rows):
        if not isinstance(row, Mapping):
            _fail(code, f"row {pos} is not a mapping")
        for key in ("base_accession", "uniprot_num", "cell"):
            if key not in row:
                _fail(code, f"row {pos} missing key {key!r}")
        if not isinstance(row["base_accession"], str) or not row["base_accession"].strip():
            _fail(code, f"row {pos}: base_accession={row['base_accession']!r}")
        num = row["uniprot_num"]
        if isinstance(num, bool) or not isinstance(num, int):
            _fail(code, f"row {pos}: uniprot_num={num!r} is not an int")
        cell = row["cell"]
        if not isinstance(cell, str) or cell not in ALL_CELLS:
            _fail(code, f"row {pos}: cell={cell!r} outside {sorted(ALL_CELLS)}")
        acc = pool.base_accession(row["base_accession"])
        identity_sets[cell].add((acc, int(num)))
        row_counts[cell] += 1
    floor_eval = pool.evaluate_support_floors(list(label_rows), dict(group_map),
                                              fold_of_group=fold_of_group)
    verdict = pool.pool_verdict(floor_eval, cap_bound=cap_bound)
    cell_residue_min = int(floor_eval["parameters"]["cell_residue_min"])
    cells_out: dict[str, dict[str, Any]] = {}
    for c in ALL_CELLS:
        # pool floors cover the four registered cells; U has no floor.
        floor_cell = floor_eval["cells"].get(
            c, {"proteins": 0, "groups": 0, "unique_residues": 0})
        cells_out[c] = {
            "row_count": row_counts[c],
            "unique_residue_identities": len(identity_sets[c]),
            "unique_residue_identities_meet_floor": len(identity_sets[c]) >= cell_residue_min,
            "floor_unique_residues": floor_cell["unique_residues"],
            "floor_proteins": floor_cell["proteins"],
            "floor_groups": floor_cell["groups"],
        }
    payload: dict[str, Any] = {
        "schema": CENSUS_SCHEMA,
        "emitted_by": EMITTED_BY,
        "registration": "docs/e420_registration.md",
        "label_rules": "docs/g1_successor_label_rules.json",
        "row_count": len(label_rows),
        "group_count": len({gid for gid in dict(group_map).values()}),
        "cells": cells_out,
        "proteins": floor_eval["proteins"],
        "groups": floor_eval["groups"],
        "fold_of_group": floor_eval["fold_of_group"],
        "folds": floor_eval["folds"],
        "floor_parameters": floor_eval["parameters"],
        "floor_checks": floor_eval["checks"],
        "floors_satisfied": floor_eval["satisfied"],
        "anomalies": floor_eval["anomalies"],
        "cap_bound": bool(cap_bound),
        "verdict": verdict,
    }
    payload["census_digest"] = pool.canonical_digest(payload)
    return payload
