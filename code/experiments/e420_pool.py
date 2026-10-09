"""e420 gate (f): predictor-blind pool construction over the availability manifest.

Registration: docs/e420_registration.md sections 4.2-4.5, 6, 8, 10(f), 12.
Selection: docs/g1_successor_analysis_plan.md section 3 (accession-level
deterministic order, cap 600). Fold rules: docs/g1_successor_fold_rules.json
(schema g1-successor-fold-rules-v1). Manifest row schema:
experiments/e420_census.py (e420-census-availability-manifest-v1). Fold
preimage helpers: experiments/e420_grammar_validate.py (fold_hex,
build_fold_preimage_row, aggregate_fold_preimage_hash, validate_fold_preimage).

OFFLINE, PURE, DETERMINISTIC. No network, no source access, no timestamps in
any produced artifact. Every function is a pure function of its arguments:
same inputs -> byte-identical outputs.

FROZEN-PARAMETER FLAG (registered prefilter). The Smith-Waterman edge rule is
registered exactly (match=2, mismatch=-1, gap=-2, leftmost tie break; edge iff
identity > 0.25 AND aligned_columns >= 100). To keep all-pairs grouping inside
a CPU-minutes envelope on ~600-2000 distinct sequences, Smith-Waterman runs
ONLY on candidate pairs chosen by a cheap amino-acid k-mer prefilter. The
prefilter threshold is A FROZEN ANALYSIS-PLAN PARAMETER owned by main
(registration section 13 O5 class): default PREFILTER_JACCARD_MIN = 0.10 with
k = PREFILTER_K = 6. Sensitivity note bound into this flag: an aa 6-mer
Jaccard of 0.10 requires roughly >= 0.69 full-sequence identity
(0.69**6 ~= 0.107), while the registered edge rule can accept local
alignments at lower identity (and shared short domains inside otherwise
distant sequences), so the DEFAULT threshold is a declared OPEN parameter that
main MUST freeze (possibly lowering it or switching prefilter_mode to
"shared_kmer", or auditing with prefilter_mode="none") BEFORE gate (f) seals
POOL_FREEZE. The prefilter is recorded in every freeze artifact together with
its per-pair audit rows so the freeze is falsifiable.

SEQUENCE UNITS. Smith-Waterman clustering runs on DISTINCT sequences only:
exact-duplicate canonical sequence strings are deduplicated first (same
sequence string => one unit carrying all of its member base accessions), so
the DP never runs twice over an identical pair. Accession -> sequence pairs
come from AFDB-metadata-successful census rows; because the census manifest
retains identity fields only (registration 4.1 parse whitelist), the canonical
sequence lookup is injected by the gate-(f) runner from the census's already
captured AFDB raw bytes (no new fetch; a canonical sequence is candidate
identity, never a predictor value, so this respects fence 11.6).

GROUPS / FOLDS / FLOORS (registration section 6, carried rules):
  * edges -> union-find connected components; group_id = lexical component
    minimum over member base accessions; groups ordered lexically;
  * fold order = SHA256("g1-successor-fold-v1:" || group_id) lexical hex
    (helpers reused from e420_grammar_validate);
  * fold assignment (amended 2026-09-04, review m9) is WITHIN EACH PRIMARY
    CLASS (A, then B): sort that class's groups by fold hex, assign sorted
    position i modulo 5 (guarantees 4+4 per-fold floors when the cohort floor
    holds); OC and O groups: global hex order with one shared counter, i mod 5;
  * POOL_FREEZE seals the label-independent preimage (selection, group map,
    per-group fold hex, rule text, typed exclusions); per-class fold NUMBERS
    are computed post-label by the pure function assign_folds, which has no
    freedom once classes are known;
  * support-floor evaluation carries the e418 section 5.1 member rules: a
    protein gets a cell only with >= N_RES_PRIMARY unique non-U residue
    identities all in one cell (else protein_ambiguous / unresolved); a group
    counts only if every member protein is mandatory, unambiguous and shares
    the group cell (else group_unresolved; the whole group counts in no cell
    and a member can never be selectively dropped).

Tie-break interpretation (registered "leftmost", documented here because the
fold-rules JSON names it without elaboration): forward-score argmax picks the
leftmost column, then the topmost row; traceback prefers diagonal, then left,
then up among predecessors achieving the cell value. Deterministic and
recorded verbatim in the freeze artifact; flagged as an interpretation OPEN
item for main's freeze confirmation.

Typed verdicts (registration section 8): pool_verdict maps floor evaluation +
cap binding to E420_H1_SUPPORT_PASS / E420_UNDERPOWERED_POOL (deterministic
cap-induced shortage; no post-hoc enlargement permitted) /
E420_INFEASIBLE_CURRENT_POOL. Typed exclusion records (never silent):
sequence_unavailable, sequence_invalid, alignment_too_large, group_unresolved,
prior-pool overlap assertion VIOLATED.
"""
from __future__ import annotations

import hashlib
import json
import math
from collections import Counter, defaultdict
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

import e420_grammar_validate as gv

# --------------------------------------------------------------------------
# Registered constants (registration 4.2/6; analysis plan 3; fold rules v1).
# --------------------------------------------------------------------------
POOL_SCHEMA = "g1-successor-pool-v1"
FREEZE_SCHEMA = "e420-pool-freeze-v1"
PRIOR_POOL_SCHEMA = "e420-prior-pool-comparison-v1"
FLOORS_SCHEMA = "e420-support-floor-evaluation-v1"
EMITTED_BY = "e420-pool-v1"

SELECTION_CAP = 600
SELECTION_SEED_PREFIX = "g1-successor-pool-v1:seed:"

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_SCORE = -2
EDGE_IDENTITY_GT = 0.25
EDGE_ALIGNED_COLUMNS_GTE = 100

PREFILTER_K = 6
PREFILTER_JACCARD_MIN = 0.10  # OPEN frozen analysis-plan parameter (see docstring)
PREFILTER_JACCARD_MODE = "jaccard"
PREFILTER_SHARED_KMER_MODE = "shared_kmer"
PREFILTER_NONE_MODE = "none"
PREFILTER_SHARED_KMER_MIN = 8  # alternative mode default; also OPEN pending freeze
# Hard memory guard (AGENTS section 5: a runaway dies alone). int16 matrix.
SW_MAX_CELLS = 64_000_000

N_RES_PRIMARY = 10
COHORT_FLOOR = {"A": {"proteins": 20, "groups": 20}, "B": {"proteins": 20, "groups": 20}}
PER_FOLD_MIN_GROUPS = {"A": 4, "B": 4}
FOLD_COUNT = 5
CELL_RESIDUE_FLOOR = 200
PRIMARY_CLASSES = ("A", "B")
ALL_CELLS = ("OC", "A", "B", "O")

# Minimal census manifest row contract (experiments/e420_census.py ROW_FIELDS;
# rows are availability evidence: identity fields only, never values).
ROW_REQUIRED_KEYS = frozenset({"row_key", "status", "availability", "accessions"})

VERDICT_H1_SUPPORT_PASS = "E420_H1_SUPPORT_PASS"
VERDICT_UNDERPOWERED = "E420_UNDERPOWERED_POOL"
VERDICT_INFEASIBLE = "E420_INFEASIBLE_CURRENT_POOL"

EXCL_SEQUENCE_UNAVAILABLE = "sequence_unavailable"
EXCL_SEQUENCE_INVALID = "sequence_invalid"
EXCL_ALIGNMENT_TOO_LARGE = "alignment_too_large"
EXCL_GROUP_UNRESOLVED = "group_unresolved"

TIE_BREAK_RULE = ("leftmost: forward argmax = leftmost column then topmost row; "
                  "traceback precedence diagonal > left > up")

FOLD_ASSIGNMENT_RULE = ("WITHIN EACH PRIMARY CLASS (A, then B): sort that class's groups by "
                        "fold_hex=SHA256('g1-successor-fold-v1:'||group_id) lexical order, assign "
                        "sorted position i modulo 5; OC and O groups: global hex order with one "
                        "shared counter, i mod 5 (presence floors only)")


class PoolError(ValueError):
    """Typed pool-construction input error (never a silent fallback)."""

    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}")


# --------------------------------------------------------------------------
# Canonical helpers (same rule as registration section 12 / grammar module).
# --------------------------------------------------------------------------
def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_digest(value: Any) -> str:
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                           allow_nan=False)
    return sha256_hex(canonical.encode("utf-8"))


def base_accession(accession: Any) -> str:
    """Canonical UniProt base accession: uppercase, one trailing -N isoform
    suffix removed (e418 section 5.1 normalization precedent)."""
    if not isinstance(accession, str) or not accession.strip():
        raise PoolError("E420_POOL_INPUT_INVALID", f"accession={accession!r}")
    acc = accession.strip().upper()
    tail = acc.rsplit("-", 1)
    if len(tail) == 2 and tail[1].isdigit() and tail[1] != "":
        acc = tail[0]
    return acc


def derive_selection_seed(checkpoint_sha256: str) -> str:
    """Registered: SHA256('g1-successor-pool-v1:seed:' || checkpoint_sha256)."""
    if not isinstance(checkpoint_sha256, str) or len(checkpoint_sha256) != 64:
        raise PoolError("E420_POOL_INPUT_INVALID", "checkpoint_sha256 must be 64-hex")
    return sha256_hex(f"{SELECTION_SEED_PREFIX}{checkpoint_sha256}".encode("utf-8"))


# --------------------------------------------------------------------------
# 1. Sequence set construction (exact-duplicate dedup, distinct sequences only).
# --------------------------------------------------------------------------
def afdb_success_accessions(rows: Sequence[Mapping[str, Any]]) -> dict[str, list[str]]:
    """accession -> sorted row_keys of AFDB-metadata-successful manifest rows.

    AFDB-metadata-successful = row status "ok" AND availability.afdb_prediction
    true (the census row schema of experiments/e420_census.py).
    """
    if not isinstance(rows, (list, tuple)):
        raise PoolError("E420_POOL_INPUT_INVALID", "manifest rows must be a sequence")
    out: dict[str, list[str]] = {}
    for pos, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise PoolError("E420_POOL_INPUT_INVALID", f"manifest row {pos} not a mapping")
        missing_keys = sorted(ROW_REQUIRED_KEYS - set(row))
        if missing_keys:
            raise PoolError("E420_POOL_INPUT_INVALID",
                            f"manifest row {pos} missing {missing_keys} (census ROW_FIELDS)")
        availability = row.get("availability") or {}
        if row.get("status") != "ok" or not availability.get("afdb_prediction"):
            continue
        for raw_acc in row.get("accessions") or []:
            acc = base_accession(raw_acc)
            out.setdefault(acc, [])
            row_key = row.get("row_key")
            if isinstance(row_key, str) and row_key and row_key not in out[acc]:
                out[acc].append(row_key)
    return {acc: sorted(keys) for acc, keys in sorted(out.items())}


def build_sequence_units(
    rows: Sequence[Mapping[str, Any]],
    canonical_sequences: Mapping[str, str],
) -> tuple[list[dict[str, Any]], list[dict[str, str]]]:
    """Unique (base_accession -> canonical sequence) pairs from AFDB-metadata-
    successful rows; EXACT-duplicate sequences are deduplicated first (same
    sequence string => one unit carrying all member accessions). Smith-Waterman
    clustering runs on these DISTINCT sequences only.

    canonical_sequences: runner-injected lookup from the census's retained
    AFDB raw bytes (identity evidence, not a predictor value).

    Returns (units, typed_exclusions). Units are sorted by sequence sha256 for
    determinism; each unit = {unit_id, sequence, sequence_sha256, accessions}.
    Typed exclusions (registration 4.5; retained, never silent):
    sequence_unavailable (AFDB-successful accession with no injected
    sequence), sequence_invalid (empty/non-A-Z after uppercase normalization).
    """
    success = afdb_success_accessions(rows)
    exclusions: list[dict[str, str]] = []
    by_sequence: dict[str, list[str]] = defaultdict(list)
    for acc in sorted(success):
        seq = canonical_sequences.get(acc)
        if seq is None:
            exclusions.append({"reason": EXCL_SEQUENCE_UNAVAILABLE, "accession": acc,
                               "detail": "no canonical sequence in injected lookup"})
            continue
        if not isinstance(seq, str):
            exclusions.append({"reason": EXCL_SEQUENCE_INVALID, "accession": acc,
                               "detail": f"sequence type {type(seq).__name__}"})
            continue
        norm = "".join(seq.split()).upper()
        if not norm or not all("A" <= ch <= "Z" for ch in norm):
            exclusions.append({"reason": EXCL_SEQUENCE_INVALID, "accession": acc,
                               "detail": "empty or non-A-Z after normalization"})
            continue
        by_sequence[norm].append(acc)
    units: list[dict[str, Any]] = []
    for seq in sorted(by_sequence):
        accs = sorted(by_sequence[seq])
        units.append({"unit_id": sha256_hex(seq.encode("ascii")),
                      "sequence": seq,
                      "sequence_sha256": sha256_hex(seq.encode("ascii")),
                      "accessions": accs})
    return units, exclusions


# --------------------------------------------------------------------------
# 2a. Conservative k-mer prefilter (FROZEN-PARAMETER; see module docstring).
# --------------------------------------------------------------------------
def kmer_set(sequence: str, k: int) -> frozenset[str]:
    if k <= 0:
        raise PoolError("E420_POOL_INPUT_INVALID", "k must be positive")
    if len(sequence) < k:
        return frozenset()
    return frozenset(sequence[i:i + k] for i in range(len(sequence) - k + 1))


def prefilter_candidates(
    units: Sequence[Mapping[str, Any]],
    *,
    k: int = PREFILTER_K,
    threshold: float = PREFILTER_JACCARD_MIN,
    mode: str = PREFILTER_JACCARD_MODE,
) -> tuple[list[tuple[int, int, float, int]], dict[str, Any]]:
    """Cheap conservative candidate generation BEFORE Smith-Waterman.

    Returns (candidates, stats). candidates = sorted (i, j, score, shared)
    with i < j unit indices; score is the mode-specific k-mer similarity.
    Modes: "jaccard" (|A∩B| / |A∪B| >= threshold), "shared_kmer"
    (|A∩B| >= threshold), "none" (every pair; score 0.0; registered
    all-pairs semantics for audit reproduction).

    Pairs whose shorter sequence is below EDGE_ALIGNED_COLUMNS_GTE can never
    satisfy the registered aligned_columns >= 100 edge rule, so they are
    dropped before Smith-Waterman in every mode (exact, not heuristic).
    """
    if mode not in (PREFILTER_JACCARD_MODE, PREFILTER_SHARED_KMER_MODE, PREFILTER_NONE_MODE):
        raise PoolError("E420_POOL_INPUT_INVALID", f"prefilter mode {mode!r}")
    n = len(units)
    lengths = [len(u["sequence"]) for u in units]
    kmers = [kmer_set(u["sequence"], k) for u in units] if mode != PREFILTER_NONE_MODE else []
    candidates: list[tuple[int, int, float, int]] = []
    if mode == PREFILTER_NONE_MODE:
        for i in range(n):
            for j in range(i + 1, n):
                if min(lengths[i], lengths[j]) < EDGE_ALIGNED_COLUMNS_GTE:
                    continue
                candidates.append((i, j, 0.0, 0))
        stats = {"mode": mode, "k": k, "threshold": threshold,
                 "candidate_pairs": len(candidates), "distinct_sequences": n,
                 "pruned_by_min_length": None}
        return candidates, stats
    postings: dict[str, list[int]] = defaultdict(list)
    for idx, ks in enumerate(kmers):
        for km in ks:
            postings[km].append(idx)
    shared: Counter = Counter()
    for plist in postings.values():
        if len(plist) < 2:
            continue
        for a_pos in range(len(plist)):
            i = plist[a_pos]
            for j in plist[a_pos + 1:]:
                shared[(i, j)] += 1
    for (i, j), sh in sorted(shared.items()):
        if min(lengths[i], lengths[j]) < EDGE_ALIGNED_COLUMNS_GTE:
            continue
        union = len(kmers[i]) + len(kmers[j]) - sh
        if mode == PREFILTER_JACCARD_MODE:
            score = (sh / union) if union else 0.0
            if score >= threshold:
                candidates.append((i, j, score, sh))
        else:
            score = float(sh)
            if sh >= threshold:
                candidates.append((i, j, score, sh))
    stats = {"mode": mode, "k": k, "threshold": threshold,
             "candidate_pairs": len(candidates),
             "distinct_sequences": n,
             "all_pairs": n * (n - 1) // 2,
             "min_length_gate": EDGE_ALIGNED_COLUMNS_GTE}
    return candidates, stats


# --------------------------------------------------------------------------
# 2b. Registered Smith-Waterman (match 2, mismatch -1, gap -2, leftmost ties).
# --------------------------------------------------------------------------
def smith_waterman(seq_a: str, seq_b: str) -> dict[str, Any]:
    """Registered local alignment. Returns score, matches, aligned_columns
    (columns where BOTH sequences contribute a residue), identity
    (matches / aligned_columns), and end cells. Deterministic tie-break per
    TIE_BREAK_RULE. Raises PoolError(E420_ALIGNMENT_TOO_LARGE) above
    SW_MAX_CELLS (typed, never a silent skip)."""
    n, m = len(seq_a), len(seq_b)
    if n == 0 or m == 0:
        return {"score": 0, "matches": 0, "aligned_columns": 0, "identity": 0.0,
                "i_end": 0, "j_end": 0}
    if n * m > SW_MAX_CELLS:
        raise PoolError("E420_ALIGNMENT_TOO_LARGE", f"cells={n * m}")
    dtype = np.int16 if 2 * min(n, m) + abs(GAP_SCORE) * max(n, m) < 32767 else np.int32
    c1 = np.frombuffer(seq_a.encode("ascii"), dtype=np.uint8)
    c2 = np.frombuffer(seq_b.encode("ascii"), dtype=np.uint8)
    h = np.zeros((n + 1, m + 1), dtype=dtype)
    neg_inf = -np.iinfo(np.int64).max // 4
    j_col = np.arange(1, m + 1, dtype=np.int64)
    coef = j_col * (-GAP_SCORE)  # A0[k] - gap*k (gap negative -> plus 2k)
    for i in range(1, n + 1):
        prev = h[i - 1]
        eq = c2 == c1[i - 1]
        sub = prev[:-1].astype(np.int64) + np.where(eq, MATCH_SCORE, MISMATCH_SCORE)
        up = prev[1:].astype(np.int64) + GAP_SCORE
        a0 = np.maximum(np.maximum(sub, up), 0)
        run = np.maximum.accumulate(a0 + coef)
        prev_run = np.empty(m, dtype=np.int64)
        prev_run[0] = neg_inf
        prev_run[1:] = run[:-1]
        left = prev_run + GAP_SCORE * j_col  # R[j0-1] + gap*j, j = j0+1
        h[i, 1:] = np.maximum(a0, left)
    col_max = h.max(axis=0)
    j_end = int(col_max.argmax())          # leftmost column of the max
    i_end = int(h[:, j_end].argmax())      # topmost row within that column
    score = int(h[i_end, j_end])
    matches = aligned = 0
    i, j = i_end, j_end
    while i > 0 and j > 0 and h[i, j] > 0:
        cell = int(h[i, j])
        diag = int(h[i - 1, j - 1]) + (MATCH_SCORE if c1[i - 1] == c2[j - 1] else MISMATCH_SCORE)
        left_v = int(h[i, j - 1]) + GAP_SCORE
        up_v = int(h[i - 1, j]) + GAP_SCORE
        if diag == cell:
            aligned += 1
            if c1[i - 1] == c2[j - 1]:
                matches += 1
            i, j = i - 1, j - 1
        elif left_v == cell:
            j -= 1
        elif up_v == cell:
            i -= 1
        else:  # pragma: no cover - DP invariant
            raise PoolError("E420_ALIGNMENT_INTERNAL", f"traceback stuck at {(i, j)}")
    identity = (matches / aligned) if aligned else 0.0
    return {"score": score, "matches": matches, "aligned_columns": aligned,
            "identity": identity, "i_end": i_end, "j_end": j_end}


def edge_decision(alignment: Mapping[str, Any]) -> bool:
    """Registered edge rule: identity > 0.25 AND aligned_columns >= 100."""
    return (bool(alignment["aligned_columns"]) and alignment["aligned_columns"] >= EDGE_ALIGNED_COLUMNS_GTE
            and alignment["identity"] > EDGE_IDENTITY_GT)


# --------------------------------------------------------------------------
# 2c. Homology grouping: SW over candidates + union-find components.
# --------------------------------------------------------------------------
class _UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            if rb < ra:
                ra, rb = rb, ra
            self.parent[rb] = ra


def homology_groups(
    units: Sequence[Mapping[str, Any]],
    *,
    prefilter_mode: str = PREFILTER_JACCARD_MODE,
    prefilter_k: int = PREFILTER_K,
    prefilter_threshold: float = PREFILTER_JACCARD_MIN,
) -> dict[str, Any]:
    """Union-find connected components over DISTINCT sequences under the
    registered SW edge rule. group_id = lexical component minimum over member
    base accessions; groups ordered lexically. Returns group map, per-group
    members, SW edge audit rows (candidate, evaluated, edge), prefilter stats,
    and typed exclusions for pairs that exceed the memory cell guard
    (conservative: a skipped pair can only fail to MERGE, never merges)."""
    candidates, pre_stats = prefilter_candidates(
        units, k=prefilter_k, threshold=prefilter_threshold, mode=prefilter_mode)
    uf = _UnionFind(len(units))
    edges: list[dict[str, Any]] = []
    exclusions: list[dict[str, str]] = []
    for i, j, pre_score, shared in candidates:
        acc_i, acc_j = units[i]["accessions"][0], units[j]["accessions"][0]
        pair_key = sorted((acc_i, acc_j))
        try:
            aln = smith_waterman(units[i]["sequence"], units[j]["sequence"])
        except PoolError as exc:
            if exc.code == "E420_ALIGNMENT_TOO_LARGE":
                exclusions.append({"reason": EXCL_ALIGNMENT_TOO_LARGE,
                                   "accessions": ",".join(pair_key),
                                   "detail": exc.detail})
                continue
            raise
        is_edge = edge_decision(aln)
        if is_edge:
            uf.union(i, j)
        edges.append({"accession_a": pair_key[0], "accession_b": pair_key[1],
                      "score": aln["score"], "matches": aln["matches"],
                      "aligned_columns": aln["aligned_columns"],
                      "identity": aln["identity"], "edge": is_edge,
                      "prefilter_score": pre_score, "shared_kmers": shared})
    comp_members: dict[int, list[str]] = defaultdict(list)
    for idx, unit in enumerate(units):
        comp_members[uf.find(idx)].extend(unit["accessions"])
    groups: dict[str, list[str]] = {}
    group_of_accession: dict[str, str] = {}
    for members in comp_members.values():
        gid = min(members)  # lexical component minimum
        groups[gid] = sorted(set(members))
    for gid in sorted(groups):
        for acc in groups[gid]:
            group_of_accession[acc] = gid
    groups = {gid: groups[gid] for gid in sorted(groups)}
    edges.sort(key=lambda e: (e["accession_a"], e["accession_b"]))
    evaluated = sum(1 for e in edges)
    return {
        "group_of_accession": group_of_accession,
        "groups": groups,
        "edges": edges,
        "prefilter_stats": pre_stats,
        "typed_exclusions": exclusions,
        "sw_pair_stats": {"candidate_pairs": len(candidates), "evaluated_pairs": evaluated,
                          "edge_pairs": sum(1 for e in edges if e["edge"]),
                          "distinct_sequences": len(units)},
    }


# --------------------------------------------------------------------------
# 3. Fold order and within-class fold assignment (registration section 6).
# --------------------------------------------------------------------------
def fold_sorted_group_ids(group_ids: Iterable[str]) -> list[str]:
    """Groups sorted by fold_hex=SHA256('g1-successor-fold-v1:'||group_id)."""
    ids = list(group_ids)
    if any(not isinstance(g, str) or not g for g in ids):
        raise PoolError("E420_POOL_INPUT_INVALID", "group_id must be nonempty str")
    return sorted(set(ids), key=gv.fold_hex)


def assign_folds(group_ids: Iterable[str],
                 class_of_group: Mapping[str, str]) -> dict[str, int]:
    """Pure post-label fold assignment (registration section 6 amendment,
    review m9). Primary classes A and B: within-class fold-hex sort, position
    i modulo 5. OC and O: ONE shared global counter over the combined
    fold-hex sort, i modulo 5. Deterministic; no freedom once classes are
    known. Classes outside {A,B,OC,O} are a typed input error."""
    ids = sorted(set(group_ids))
    missing = [g for g in ids if g not in class_of_group]
    if missing:
        raise PoolError("E420_POOL_INPUT_INVALID", f"class missing for groups {missing[:5]}")
    bad = sorted({class_of_group[g] for g in ids} - set(ALL_CELLS))
    if bad:
        raise PoolError("E420_POOL_INPUT_INVALID", f"class outside {{OC,A,B,O}}: {bad[:5]}")
    folds: dict[str, int] = {}
    for cls in PRIMARY_CLASSES:
        members = fold_sorted_group_ids(g for g in ids if class_of_group[g] == cls)
        for i, gid in enumerate(members):
            folds[gid] = i % FOLD_COUNT
    presence = fold_sorted_group_ids(g for g in ids
                                     if class_of_group[g] in ("OC", "O"))
    for i, gid in enumerate(presence):
        folds[gid] = i % FOLD_COUNT
    return folds


# --------------------------------------------------------------------------
# 4. Deterministic selection (analysis plan section 3; registration 4.2).
# --------------------------------------------------------------------------
def selection_order_key(selection_seed: str, accession: str) -> str:
    """SHA256(selection_seed || '|' || base_accession)."""
    if not isinstance(selection_seed, str) or not selection_seed:
        raise PoolError("E420_POOL_INPUT_INVALID", "selection_seed empty")
    return sha256_hex(f"{selection_seed}|{base_accession(accession)}".encode("utf-8"))


def select_pool(accessions: Iterable[str], selection_seed: str,
                cap: int = SELECTION_CAP) -> dict[str, Any]:
    """Lexical order_key over eligible accessions; admit in order up to the
    registered 600-accession cap. ALL census rows of an admitted accession
    travel with it (group-atomicity is enforced by folds, not selection).
    cap_bound=True means the cap truncated the eligible set."""
    if cap < 1:
        raise PoolError("E420_POOL_INPUT_INVALID", "cap must be >= 1")
    ordered = sorted(((selection_order_key(selection_seed, acc), acc)
                      for acc in set(accessions)), key=lambda t: (t[0], t[1]))
    admitted = [acc for _, acc in ordered[:cap]]
    order_keys = {acc: key for key, acc in ordered}
    return {"admitted": admitted, "order_keys": order_keys,
            "eligible_count": len(ordered), "cap": cap,
            "cap_bound": len(ordered) > cap,
            "overflow_count": max(0, len(ordered) - cap)}


# --------------------------------------------------------------------------
# 5. Prior-pool exclusion interface (registration 4.3).
# --------------------------------------------------------------------------
def prior_pool_intersection(pool_identifiers: Iterable[str],
                            prior_pools: Mapping[str, Iterable[str]]) -> dict[str, Any]:
    """Intersection-check receipt against the preregistered prior-pool
    identifier sets (e413/CASP targets, e394 attempt identifiers, T035
    fixtures). Normalization precedent: docs/e418_prior_pool_sentinel.json
    (uppercase; base-accession folding via base_accession). The registered
    assertion is an EMPTY intersection; a violation is a typed no-go signal
    recorded in the receipt (prior_pool_overlap_violated), never repaired."""
    pool_ids = sorted({i.strip().upper() for i in pool_identifiers if isinstance(i, str) and i.strip()})
    pool_base = {base_accession(i) for i in pool_ids}
    hits: list[dict[str, str]] = []
    counts: dict[str, int] = {}
    for source in sorted(prior_pools):
        raw_ids = sorted({i.strip().upper() for i in prior_pools[source]
                          if isinstance(i, str) and i.strip()})
        counts[source] = len(raw_ids)
        raw_set = set(raw_ids)
        for pid in pool_ids:
            if pid in raw_set or base_accession(pid) in {base_accession(r) for r in raw_set}:
                hits.append({"source": source, "identifier": pid})
    hits.sort(key=lambda h: (h["source"], h["identifier"]))
    receipt = {
        "schema": PRIOR_POOL_SCHEMA,
        "pool_identifier_count": len(pool_ids),
        "prior_pool_identifier_counts": counts,
        "intersection": hits,
        "intersection_count": len(hits),
        "empty": len(hits) == 0,
        "empty_intersection_assertion": "HELD" if not hits else "VIOLATED",
        "normalization": "uppercase + base-accession folding (e418 sentinel precedent)",
        "projection_sha256": canonical_digest({"pool": pool_ids, "priors": counts}),
    }
    return receipt


# --------------------------------------------------------------------------
# 5b. Group-atomic member rule: group_unresolved on unexplained members.
# --------------------------------------------------------------------------
def unresolved_group_members(rows: Sequence[Mapping[str, Any]],
                             admitted_accessions: Iterable[str],
                             group_of_accession: Mapping[str, str]) -> list[dict[str, str]]:
    """Registration section 6 (carried e418 5.1): every selected member is
    mandatory. A manifest row of an admitted accession that cannot be
    explained into the group structure (accession without a group, or a row
    resolving to an admitted accession that has no group) is a typed
    group_unresolved exclusion. Records are deduplicated and sorted."""
    admitted = set(admitted_accessions)
    records: set[tuple[str, str, str]] = set()
    has_rows: dict[str, bool] = {acc: False for acc in admitted}
    for row in rows:
        accs = [base_accession(a) for a in (row.get("accessions") or [])]
        hit = False
        for acc in accs:
            if acc in admitted:
                has_rows[acc] = True
                hit = True
        if not hit:
            continue
        row_key = str(row.get("row_key") or "")
        for acc in sorted(set(accs)):
            if acc in admitted and acc not in group_of_accession:
                records.add((EXCL_GROUP_UNRESOLVED, acc, f"member_row_unexplained:{row_key}"))
    for acc in sorted(admitted):
        if not has_rows[acc]:
            records.add((EXCL_GROUP_UNRESOLVED, acc, "no_manifest_rows"))
        elif acc not in group_of_accession:
            records.add((EXCL_GROUP_UNRESOLVED, acc, "no_sequence_group"))
    return [{"reason": r, "accession": a, "detail": d} for r, a, d in sorted(records)]


# --------------------------------------------------------------------------
# 6. Support-floor evaluator (registration section 6; carried e418 5.1 rules).
# --------------------------------------------------------------------------
def evaluate_support_floors(
    label_rows: Iterable[Mapping[str, Any]],
    group_map: Mapping[str, str],
    *,
    class_of_group: Mapping[str, str] | None = None,
    fold_of_group: Mapping[str, int] | None = None,
    n_res_primary: int = N_RES_PRIMARY,
    cohort_min: int = None,
    per_fold_min: Mapping[str, int] | None = None,
    cell_residue_min: int = CELL_RESIDUE_FLOOR,
) -> dict[str, Any]:
    """Evaluate the registered support floors over per-residue label rows.

    label_rows: {base_accession, uniprot_num, cell} with cell in
    {OC,A,B,O,U}; U rows are ignored for floors (counted as U).
    group_map: accession -> group_id (the frozen pool group map).
    Folds: give fold_of_group directly, or class_of_group and the evaluator
    applies the pure assign_folds rule. Carried e418 5.1 member rules:
    a protein resolves to a cell only with >= n_res_primary unique non-U
    identities all in ONE cell (mixed -> protein_ambiguous, counts nowhere);
    a group counts only if EVERY member resolves to that same cell and the
    group holds >= n_res_primary unique non-U identities (else
    group_unresolved; no member may be selectively dropped)."""
    if cohort_min is None:
        cohort_min = COHORT_FLOOR["A"]["proteins"]
    if per_fold_min is None:
        per_fold_min = dict(PER_FOLD_MIN_GROUPS)
    identities: dict[tuple[str, str], set[int]] = defaultdict(set)
    identity_cells: dict[tuple[str, int], set[str]] = defaultdict(set)
    anomalies: set[tuple[str, str]] = set()
    for row in label_rows:
        acc = base_accession(row["base_accession"])
        num = int(row["uniprot_num"])
        cell = str(row["cell"]).upper()
        if acc not in group_map:
            anomalies.add(("label_row_unmapped", f"{acc}:{num}"))
            continue
        if cell == "U":
            continue
        if cell not in ALL_CELLS:
            raise PoolError("E420_POOL_INPUT_INVALID", f"cell={cell!r}")
        identities[(acc, cell)].add(num)
        identity_cells[(acc, num)].add(cell)
    proteins: dict[str, Any] = {}
    for acc in sorted(group_map):
        cells_present = [c for c in ALL_CELLS if identities.get((acc, c))]
        unique_non_u = {c: len(identities.get((acc, c)) or ()) for c in ALL_CELLS}
        if len(cells_present) > 1:
            status, cell = "protein_ambiguous", None
        elif (len(cells_present) == 1
                and len(identities[(acc, cells_present[0])]) >= n_res_primary):
            status, cell = "resolved", cells_present[0]
        else:
            status, cell = "unresolved", None
        proteins[acc] = {"status": status, "cell": cell, "unique_non_u": unique_non_u}
    members_by_group: dict[str, list[str]] = defaultdict(list)
    for acc, gid in group_map.items():
        members_by_group[gid].append(acc)
    groups_out: dict[str, Any] = {}
    cell_of_group: dict[str, str] = {}
    for gid in sorted(members_by_group):
        members = sorted(members_by_group[gid])
        member_cells = {proteins[a]["cell"] for a in members}
        member_set = set(members)
        group_unique = len({(a, n) for (a, n), cs in identity_cells.items()
                            if a in member_set and cs})
        if (len(member_cells) == 1 and None not in member_cells
                and group_unique >= n_res_primary):
            cell = next(iter(member_cells))
            cell_of_group[gid] = cell
            groups_out[gid] = {"cell": cell, "members": members, "status": "resolved",
                               "unique_non_u_identities": group_unique}
        else:
            cell_of_group[gid] = None
            groups_out[gid] = {"cell": None, "members": members, "status": EXCL_GROUP_UNRESOLVED,
                               "unique_non_u_identities": group_unique}
    if fold_of_group is None:
        if class_of_group is None:
            resolved_classes = {gid: cell_of_group[gid] for gid in cell_of_group
                                if cell_of_group[gid] in ALL_CELLS}
        else:
            resolved_classes = {gid: class_of_group[gid] for gid in cell_of_group
                                if gid in class_of_group and class_of_group[gid] in ALL_CELLS}
        fold_of_group = assign_folds(sorted(resolved_classes), resolved_classes)
    cells_summary: dict[str, dict[str, int]] = {}
    for c in ALL_CELLS:
        cells_summary[c] = {
            "proteins": sum(1 for p in proteins.values() if p["cell"] == c),
            "groups": sum(1 for g in groups_out.values() if g["cell"] == c),
            "unique_residues": len({(a, n) for (a, n), cs in identity_cells.items()
                                    if cs == {c} and proteins.get(a, {}).get("cell") == c}),
        }
    folds_summary: dict[str, dict[str, Any]] = {}
    per_fold_ok = True
    for f in range(FOLD_COUNT):
        a_groups = sum(1 for gid, c in cell_of_group.items()
                       if c == "A" and fold_of_group.get(gid) == f)
        b_groups = sum(1 for gid, c in cell_of_group.items()
                       if c == "B" and fold_of_group.get(gid) == f)
        ok = a_groups >= per_fold_min.get("A", 4) and b_groups >= per_fold_min.get("B", 4)
        per_fold_ok = per_fold_ok and ok
        folds_summary[str(f)] = {"A_groups": a_groups, "B_groups": b_groups, "ok": ok}
    checks = {
        "cohort_proteins_A": cells_summary["A"]["proteins"] >= cohort_min,
        "cohort_proteins_B": cells_summary["B"]["proteins"] >= cohort_min,
        "cohort_groups_A": cells_summary["A"]["groups"] >= cohort_min,
        "cohort_groups_B": cells_summary["B"]["groups"] >= cohort_min,
        "per_fold_ab_support": per_fold_ok,
        "oc_presence": (cells_summary["OC"]["proteins"] >= 1 and cells_summary["OC"]["groups"] >= 1),
        "o_presence": (cells_summary["O"]["proteins"] >= 1 and cells_summary["O"]["groups"] >= 1),
        "cell_residue_floor_A": cells_summary["A"]["unique_residues"] >= cell_residue_min,
        "cell_residue_floor_B": cells_summary["B"]["unique_residues"] >= cell_residue_min,
        "cell_residue_floor_OC": cells_summary["OC"]["unique_residues"] >= cell_residue_min,
        "cell_residue_floor_O": cells_summary["O"]["unique_residues"] >= cell_residue_min,
    }
    return {
        "schema": FLOORS_SCHEMA,
        "parameters": {"n_res_primary": n_res_primary, "cohort_min": cohort_min,
                       "per_fold_min": dict(per_fold_min), "cell_residue_min": cell_residue_min},
        "proteins": proteins,
        "groups": groups_out,
        "fold_of_group": {g: int(f) for g, f in sorted(fold_of_group.items())},
        "cells": cells_summary,
        "folds": folds_summary,
        "anomalies": sorted(anomalies),
        "checks": checks,
        "satisfied": all(checks.values()),
    }


# --------------------------------------------------------------------------
# 8. Typed verdict mapping (registration section 8).
# --------------------------------------------------------------------------
def pool_verdict(floor_evaluation: Mapping[str, Any], *, cap_bound: bool) -> str:
    """E420_H1_SUPPORT_PASS when every registered floor holds;
    E420_UNDERPOWERED_POOL on a deterministic cap-induced shortage (the cap
    truncated the eligible set AND floors fail; no post-hoc enlargement is
    permitted); E420_INFEASIBLE_CURRENT_POOL when floors fail without the cap
    binding."""
    if floor_evaluation.get("satisfied"):
        return VERDICT_H1_SUPPORT_PASS
    return VERDICT_UNDERPOWERED if cap_bound else VERDICT_INFEASIBLE


# --------------------------------------------------------------------------
# 7. POOL_FREEZE artifact builder (registration 4.2/10(f)/12; deterministic).
# --------------------------------------------------------------------------
def build_pool_freeze(
    manifest_rows: Sequence[Mapping[str, Any]],
    canonical_sequences: Mapping[str, str],
    selection_seed: str,
    manifest_sha256: str,
    prior_pools: Mapping[str, Iterable[str]],
    *,
    cap: int = SELECTION_CAP,
    prefilter_mode: str = PREFILTER_JACCARD_MODE,
    prefilter_k: int = PREFILTER_K,
    prefilter_threshold: float = PREFILTER_JACCARD_MIN,
) -> dict[str, Any]:
    """Deterministic POOL_FREEZE artifact over the availability manifest.

    Pipeline (registration 4.2-4.5): AFDB-successful accessions -> prior-pool
    exclusion receipt (empty-intersection assertion) -> sequence units with
    exact-duplicate dedup -> deterministic seeded selection under the cap ->
    homology groups over DISTINCT sequences (conservative prefilter + SW) ->
    group-atomic member check -> fold preimage via the registered grammar
    helpers -> canonical compact sorted-key JSON digest. No timestamps; no
    predictor values; the fold preimage aggregate uses
    e420_grammar_validate.aggregate_fold_preimage_hash. The artifact carries
    its own pool_freeze_digest over the payload WITHOUT that field."""
    if not isinstance(manifest_sha256, str) or len(manifest_sha256) != 64:
        raise PoolError("E420_POOL_INPUT_INVALID", "manifest_sha256 must be 64-hex")
    success = afdb_success_accessions(manifest_rows)
    receipt = prior_pool_intersection(success.keys(), prior_pools)
    excluded_accs = {h["identifier"] for h in receipt["intersection"]}
    eligible_set = {acc for acc in success if acc not in excluded_accs}
    units_all, seq_exclusions_all = build_sequence_units(manifest_rows, canonical_sequences)
    seq_exclusions = [e for e in seq_exclusions_all if e.get("accession") in eligible_set]
    # Restrict unit membership to eligible / admitted accessions; exact-duplicate
    # dedup is by sequence string, so restriction never merges or drops a unit.
    def _restrict(units_in: list[dict[str, Any]], allowed: set[str]) -> list[dict[str, Any]]:
        out_units: list[dict[str, Any]] = []
        for unit in units_in:
            members = [a for a in unit["accessions"] if a in allowed]
            if members:
                out_units.append(dict(unit, accessions=members))
        return out_units

    units = _restrict(units_all, eligible_set)
    unit_accs = [acc for unit in units for acc in unit["accessions"]]
    selection = select_pool(unit_accs, selection_seed, cap=cap)
    admitted = selection["admitted"]
    admitted_set = set(admitted)
    admitted_units = _restrict(units, admitted_set)
    grouping = homology_groups(admitted_units, prefilter_mode=prefilter_mode,
                               prefilter_k=prefilter_k, prefilter_threshold=prefilter_threshold)
    group_of_accession = grouping["group_of_accession"]
    unresolved = unresolved_group_members(manifest_rows, admitted, group_of_accession)
    fold_rows = [gv.build_fold_preimage_row(group_of_accession[acc], acc)
                 for acc in sorted(admitted) if acc in group_of_accession]
    fold_rows.sort(key=lambda r: (r["group_id"], r["accession"]))
    aggregate = gv.aggregate_fold_preimage_hash(fold_rows)
    exclusions = sorted(seq_exclusions + grouping["typed_exclusions"] + unresolved,
                        key=lambda e: (e.get("reason", ""), e.get("accession", ""),
                                       e.get("accessions", ""), e.get("detail", "")))
    payload = {
        "schema": FREEZE_SCHEMA,
        "pool_schema": POOL_SCHEMA,
        "emitted_by": EMITTED_BY,
        "registration": "docs/e420_registration.md",
        "analysis_plan": "docs/g1_successor_analysis_plan.md",
        "fold_rules": "docs/g1_successor_fold_rules.json",
        "manifest_sha256": manifest_sha256,
        "selection_seed": selection_seed,
        "selection": {"cap": selection["cap"], "eligible_count": selection["eligible_count"],
                      "admitted": admitted, "order_keys": selection["order_keys"],
                      "cap_bound": selection["cap_bound"],
                      "overflow_count": selection["overflow_count"]},
        "prior_pool_receipt": receipt,
        "sequence_units": {"unit_count": len(units),
                           "distinct_sequence_count": len(admitted_units),
                           "afdb_success_accessions": len(success),
                           "exact_duplicate_dedup_note":
                               "SW clustering runs on distinct sequences only; same sequence "
                               "string => same unit carrying all member accessions"},
        "alignment": {"algorithm": "Smith-Waterman", "match": MATCH_SCORE,
                      "mismatch": MISMATCH_SCORE, "gap": GAP_SCORE,
                      "tie_break": TIE_BREAK_RULE,
                      "edge_rule": {"identity_gt": EDGE_IDENTITY_GT,
                                    "aligned_columns_gte": EDGE_ALIGNED_COLUMNS_GTE}},
        "prefilter": {"mode": prefilter_mode, "k": prefilter_k, "threshold": prefilter_threshold,
                      "status": "OPEN_FROZEN_PARAMETER_PENDING_MAIN_FREEZE",
                      "stats": grouping["prefilter_stats"]},
        "group_map": group_of_accession,
        "groups": grouping["groups"],
        "fold_assignment_rule": FOLD_ASSIGNMENT_RULE,
        "fold_preimage_rows": fold_rows,
        "fold_preimage_aggregate_sha256": aggregate,
        "sw_pair_stats": grouping["sw_pair_stats"],
        "edge_audit": grouping["edges"],
        "typed_exclusions": exclusions,
        "floors_config": {"n_res_primary": N_RES_PRIMARY,
                          "cohort": COHORT_FLOOR, "per_fold": PER_FOLD_MIN_GROUPS,
                          "folds": FOLD_COUNT, "cell_residue_floor": CELL_RESIDUE_FLOOR},
    }
    payload["pool_freeze_digest"] = canonical_digest(payload)
    return payload


def verify_pool_freeze(artifact: Mapping[str, Any]) -> bool:
    """Recompute the digest and the fold-preimage aggregate of a freeze
    artifact (mutation-check helper for the audit)."""
    payload = {k: v for k, v in artifact.items() if k != "pool_freeze_digest"}
    if canonical_digest(payload) != artifact.get("pool_freeze_digest"):
        return False
    rows = artifact.get("fold_preimage_rows") or []
    return gv.aggregate_fold_preimage_hash(rows) == artifact.get("fold_preimage_aggregate_sha256")
