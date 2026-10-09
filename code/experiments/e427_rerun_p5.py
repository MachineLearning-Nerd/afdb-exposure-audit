"""e427 §4.3: corrected-ledger batch-4 evaluation with the REGISTERED,
unchanged evaluator (e425 run_logged_evaluation -> e424 gate -> e422 driver).

Called directly (not via the e425 CLI) because the program pause blocks the
scheduled CLI; this is the one evaluation docs/e427_registration.md authorizes
during the pause.  Inputs: results/e427/batches (corrected ledger, committed
and OTS-stamped before this run).  Outputs: results/e427/evals/.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import e425_execreceipt as e425  # noqa: E402

if __name__ == "__main__":
    art = e425.run_logged_evaluation(
        e425.REGISTERED_T0, 4, "results/e427/batches", "results/e427/evals",
        prior_binding_path="results/e427/batches/e422_binding_state.json")
    for key in ("p5_90", "p5_95"):
        print(key, art.get("binding", {}).get(key))
