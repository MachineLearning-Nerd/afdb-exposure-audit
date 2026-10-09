"""e424 binding-state integrity gate for e422 scheduled evaluations.

Allocated 2026-09-10 by w1:pG under main's standing delegation
(DECISIONS 2026-09-10T03:55+05:30).  Trigger: the binding state machine
audit found the ONE unverified link in the registered evaluation chain —
the FROZEN e422_evaldriver reads the binding sidecar with plain
json.loads and never checks its "sha256" (every other link — the
cumulative checkpoint ledger, the freeze receipt, the evaluation and
decision artifacts — is digest-enforced), and tolerates a MISSING
sidecar at batches 8/12 as an empty prior.  A tampered, stale, or
wrong-path sidecar would silently feed forged prior bindings into the
day-56/day-84 artifacts — the silent-acceptance class e423 repaired for
the census path, here on the binding path.

Per the gate-(a) freeze discipline the frozen driver is NOT touched
(any change to a bound object is a new E-number): this gate VERIFIES
before delegating and RE-VERIFIES after.

  - the registered T0 identity is ENFORCED (detached review
    2026-09-10, wf_d8bff144: the t0 argument is otherwise
    procedural-only — any t0 yields self-consistently "closed"
    windows); the registered entry points pin the frozen bind.
  - anchored evaluation outputs are NEVER overwritten: an existing
    evaluation artifact, floor-decision receipt, or execution receipt
    for the batch refuses typed (the silent-replacement class
    mechanically killed on the dispatch path by E423_CHECKPOINT_EXISTS;
    re-run = dated DECISIONS entry + manual removal).
  - a supplied sidecar path is REQUIRED at every scheduled batch: an
    omitted path is a typed refusal (at batches >= 8 the frozen driver
    would silently start a fresh study; at batch 4 it would silently
    never write the sidecar — the chain must be established, not
    implied).  Batch 4's registered fresh study is an absent FILE with
    a supplied PATH.
  - batch < 8 (day-28): the sidecar is created by this evaluation; a
    pre-existing sidecar must digest-verify AND be from a strictly
    earlier chain position — the driver's binding_map is
    first-binding-wins, so a same-or-later-position sidecar would
    silently pass its bindings into a "day-28" artifact (chain-position
    hole found by the detached review and reproduced live).
  - batch >= 8 (day-56/day-84): the sidecar file is REQUIRED — a
    missing file is a typed refusal (a chain break must never look
    like a fresh study); present -> digest-verified from a strictly
    earlier chain position, with only registered binding keys.
  - after the frozen driver returns, the (re)written sidecar is
    re-read and re-verified at its own chain position (write-path
    check).

Accepted TOCTOU window between verify and the driver's own read: the
gate catches corruption/staleness in single-operator local machinery;
it is not an adversarial control.
"""
from __future__ import annotations

import enum
import json
import os

import e422_evaldriver as driver
import e422_manifest as manifest_mod

BINDING_STATE_SCHEMA = "e422-binding-state-v1"
REQUIRED_FROM_BATCH = 8
REGISTERED_T0 = "2026-09-11T00:00:00Z"  # frozen registration bind (gate a)
# the five registered prediction keys (driver _binding_view); anything
# else in a sidecar "binding" object is not a registered chain state
REGISTERED_BINDING_KEYS = frozenset(
    {"p5_90", "p5_95", "p6", "p7a", "p7b"})


class GateError(RuntimeError):
    """Typed gate failure (``CODE: detail``), mirroring DriverError."""

    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}" if detail else code)


class GateFailureCode(str, enum.Enum):
    BINDING_STATE_MISSING = "E424_BINDING_STATE_MISSING"
    BINDING_STATE_DIGEST_MISMATCH = "E424_BINDING_STATE_DIGEST_MISMATCH"
    BINDING_STATE_SCHEMA = "E424_BINDING_STATE_SCHEMA"
    BINDING_STATE_WRITE_UNVERIFIED = "E424_BINDING_STATE_WRITE_UNVERIFIED"
    BINDING_STATE_CHAIN_POSITION = "E424_BINDING_STATE_CHAIN_POSITION"
    EVALUATION_EXISTS = "E424_EVALUATION_EXISTS"
    T0_MISMATCH = "E424_T0_MISMATCH"


def _batch_artifact_paths(result_dir: str, batch_index: int) -> list[str]:
    """Every anchored output a scheduled evaluation for this batch may
    write (evaluation artifact, batch-8 floor-decision receipt, e425
    execution receipt).  Any of them existing means the evaluation
    already ran."""
    names = [f"e422_batch{batch_index:02d}_evaluation.json"]
    if batch_index in driver.SCHEDULE and batch_index >= REQUIRED_FROM_BATCH:
        names.append(f"e422_batch{batch_index:02d}_floor_decision.json")
    names.append(f"e422_batch{batch_index:02d}_execution_receipt.json")
    return [os.path.join(result_dir, name) for name in names]


def load_verified_binding_state(
    path: str, *, required: bool, max_batch_index: int,
) -> dict | None:
    """Parse + digest-verify a binding-state sidecar (or refuse typed).

    The document must carry an integer ``batch_index`` STRICTLY BELOW
    ``max_batch_index`` (pre-driver verification passes the evaluated
    batch; post-write re-verification passes evaluated_batch + 1, since
    the driver legitimately writes the sidecar AT the evaluated batch).
    Binding keys must be a subset of the five registered prediction
    keys.  Returns the verified document, or None when not required and
    absent.
    """
    if not os.path.exists(path):
        if required:
            raise GateError(
                GateFailureCode.BINDING_STATE_MISSING.value,
                f"binding sidecar missing at {path} — batches >= "
                f"{REQUIRED_FROM_BATCH} chain a REQUIRED binding state; a "
                "missing file is a chain break, never a fresh study "
                "(e424; runbook)")
        return None
    try:
        with open(path, "rb") as handle:
            doc = json.loads(handle.read().decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise GateError(
            GateFailureCode.BINDING_STATE_SCHEMA.value,
            f"{path}: undecodable or malformed JSON "
            f"({type(exc).__name__}: {exc})") from exc
    if not isinstance(doc, dict) or doc.get("schema") != BINDING_STATE_SCHEMA:
        raise GateError(
            GateFailureCode.BINDING_STATE_SCHEMA.value,
            f"{path}: not a {BINDING_STATE_SCHEMA} document")
    digest = doc.get("sha256")
    body = {k: v for k, v in doc.items() if k != "sha256"}
    if (not isinstance(digest, str)
            or manifest_mod.canonical_sha256(body) != digest):
        raise GateError(
            GateFailureCode.BINDING_STATE_DIGEST_MISMATCH.value,
            f"{path}: sha256 does not verify over canonical content — "
            "tampered or stale sidecar, refusing (e424)")
    if not isinstance(doc.get("binding"), dict):
        raise GateError(
            GateFailureCode.BINDING_STATE_SCHEMA.value,
            f"{path}: verified document carries no binding object")
    position = doc.get("batch_index")
    if (not isinstance(position, int) or isinstance(position, bool)
            or position < 1 or position >= max_batch_index):
        raise GateError(
            GateFailureCode.BINDING_STATE_CHAIN_POSITION.value,
            f"{path}: sidecar chain position {position!r} is not a "
            f"registered prior state (position must be a batch index "
            f"strictly below {max_batch_index}; the driver's "
            "first-binding-wins map would silently import foreign "
            "bindings — detached review 2026-09-10)")
    junk = sorted(set(doc["binding"]) - REGISTERED_BINDING_KEYS)
    if junk:
        raise GateError(
            GateFailureCode.BINDING_STATE_CHAIN_POSITION.value,
            f"{path}: unregistered binding keys {junk} — only "
            f"{sorted(REGISTERED_BINDING_KEYS)} are chain state")
    return doc


def run_gated_evaluation(
    t0_utc: str,
    batch_index: int,
    out_dir: str,
    result_dir: str,
    *,
    seeds=None,
    prior_binding_path: str | None = None,
    now_iso: str | None = None,
) -> dict:
    """Digest-gated wrapper around the FROZEN run_scheduled_evaluation.

    Same signature and return; adds the verification layers above.
    Refusals precede every side effect (the driver is not called, no
    file is written).
    """
    if t0_utc != REGISTERED_T0:
        raise GateError(
            GateFailureCode.T0_MISMATCH.value,
            f"t0 {t0_utc!r} is not the frozen registration bind "
            f"{REGISTERED_T0!r} — the registered entry points pin the "
            "T0 identity (frozen text dc7a3627; detached review "
            "2026-09-10); dev shakedowns call the raw driver directly")
    existing = [p for p in _batch_artifact_paths(result_dir, batch_index)
                if os.path.exists(p)]
    if existing:
        raise GateError(
            GateFailureCode.EVALUATION_EXISTS.value,
            f"batch {batch_index}: anchored evaluation output already "
            f"exists: {[os.path.basename(p) for p in existing]} — an "
            "evaluation is never silently re-run (its bytes carry the "
            "evaluated/decided instants); a documented re-run is a "
            "dated DECISIONS entry plus MANUAL removal of the prior "
            "outputs (E423_CHECKPOINT_EXISTS doctrine)")
    if prior_binding_path is None:
        # an unsupplied path must refuse exactly like a missing file:
        # at batches >= 8 the frozen driver would silently start a fresh
        # study, and at batch 4 it would silently never WRITE the
        # sidecar (the driver only writes when the path is supplied) —
        # either way the chain is broken with nothing recording it
        raise GateError(
            GateFailureCode.BINDING_STATE_MISSING.value,
            f"batch {batch_index}: no prior_binding_path supplied — the "
            "binding sidecar path is REQUIRED at every scheduled batch "
            "(batch 4's registered fresh study is an absent FILE with a "
            "supplied PATH — the driver then creates the sidecar; an "
            "omitted argument is a chain break, never a fresh study "
            "— detached review 2026-09-10)")
    if prior_binding_path is not None:
        load_verified_binding_state(
            prior_binding_path, required=batch_index >= REQUIRED_FROM_BATCH,
            max_batch_index=batch_index)
    artifact = driver.run_scheduled_evaluation(
        t0_utc, batch_index, out_dir, result_dir, seeds=seeds,
        prior_binding_path=prior_binding_path, now_iso=now_iso)
    if prior_binding_path is not None:
        # the driver (re)wrote the sidecar AT this batch's chain
        # position — the write path verifies too
        try:
            load_verified_binding_state(
                prior_binding_path, required=True,
                max_batch_index=batch_index + 1)
        except GateError as exc:
            raise GateError(
                GateFailureCode.BINDING_STATE_WRITE_UNVERIFIED.value,
                f"{prior_binding_path}: post-write re-verification failed "
                f"({exc.code})") from exc
    return artifact
