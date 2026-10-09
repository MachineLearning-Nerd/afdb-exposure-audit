"""Assemble a scratch "work root" that reproduces the original project layout.

The registered scripts were frozen with repository-relative paths such as
results/e427/batches, paper/arxiv/verification/out and docs/e422_registration.md,
and several compute their root from their own location. Rather than edit
registered code, this tool copies the released files into the layout they
expect:

  <work>/experiments/            <- code/experiments/*.py
  <work>/tests/                  <- code/tests/*.py
  <work>/docs/                   <- registrations/*  (registration texts)
  <work>/results/ ...            <- data/derived/results/ ...
  <work>/data/paper_verification <- data/derived/data/paper_verification (RCSB cache)
  <work>/paper/arxiv/verification/{v*.py,make_figures_rev3.py,REGISTRATION.md}
  <work>/paper/arxiv/verification/{out,cache,audit}/  <- data/derived/paper/arxiv/verification/...
  <work>/paper/arxiv/verification/audit/k7pq54_mapping_audit.py
  <work>/results/e427/audit/independent_label_sample.py
  <work>/results/e427/audit/phaseB_subagent/{recompute_p5,sensitivity_p5}.py

Everything is copied (no symlinks), so running scripts in the work root can
never modify the released files. Raw third-party inputs (results/e422/raw,
data/e420/*, data/e427/e421_raw, the 9qj6/K7PQ54 mmCIF files) are not
released; see data/fetch/README.md for re-fetching them into <work>.

Usage: python -I code/tools/make_workroot.py [<work_dir>]   (default: _work)
"""
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def copy(src, dst):
    src, dst = REPO / src, Path(dst)
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def main(work):
    work = Path(work).resolve()
    if work.exists() and any(work.iterdir()):
        sys.exit(f"refusing to write into non-empty {work}")
    work.mkdir(parents=True, exist_ok=True)
    copy("data/derived", work)
    copy("code/experiments", work / "experiments")
    copy("code/tests", work / "tests")
    copy("registrations", work / "docs")
    v = work / "paper/arxiv/verification"
    for p in (REPO / "code/verification").glob("*.py"):
        copy(p.relative_to(REPO), v / p.name)
    copy("registrations/verification/REGISTRATION.md", v / "REGISTRATION.md")
    copy("code/audit/k7pq54_mapping_audit.py", v / "audit/k7pq54_mapping_audit.py")
    copy("code/audit/independent_label_sample.py", work / "results/e427/audit/independent_label_sample.py")
    for name in ("recompute_p5.py", "sensitivity_p5.py"):
        copy(f"code/audit/e427_phaseB/{name}", work / f"results/e427/audit/phaseB_subagent/{name}")
    copy("paper/figures", work / "paper/arxiv/figures")
    print(f"work root ready: {work}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "_work")
