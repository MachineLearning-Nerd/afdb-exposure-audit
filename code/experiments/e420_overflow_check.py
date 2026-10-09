
"""e420 amendment 15.3: bounded temporal-check-only capture over the 269 overflow
accessions (288 census entries). For each entry: exactly one reference cif GET
(entry-files lane); extract initial-release date; compare vs the accession's
AFDB modelCreatedDate (from ALREADY-CAPTURED census raw bytes — no AFDB fetch).
Output: measured temporal pass/fail per accession -> mechanically resolves
INFEASIBLE vs UNDERPOWERED. Schedule/receipts per analysis plan."""
import hashlib, json, os, sys, time, random, datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import requests
import e420_grammar_validate as gv

ov = json.load(open(os.path.join(ROOT, "results/e420/overflow_set.json")))
m = json.load(open(os.path.join(ROOT, "results/e420/source_availability_manifest.json")))
# modelCreatedDate per accession from census raw bytes
mcd = {}
for r in m["rows"]:
    refs = r.get("raw_refs") or {}
    ref = refs.get("E420_AFDB_METADATA_V1") or {}
    dig = ref.get("sha256")
    if not dig: continue
    path = os.path.join(ROOT, "data/e420/census_raw", dig)
    if not os.path.exists(path): continue
    try:
        arr = json.load(open(path))
        obj = sorted(arr, key=lambda o: (o.get("modelCreatedDate") or "", o.get("modelEntityId") or ""))[0]
        for a in (r.get("accessions") or []):
            mcd.setdefault(a, obj.get("modelCreatedDate"))
    except Exception:
        pass
OUT = os.path.join(ROOT, "results/e420/overflow_temporal_check.json")
results = {}
count = {"n": 0}
def late(): time.sleep(max(0.0, 0.25 + random.uniform(-0.1, 0.1)))
def get_cif(entry):
    if count["n"] >= 400: return None
    count["n"] += 1; late()
    url = f"https://www.ebi.ac.uk/pdbe/entry-files/{entry}.cif"
    resp = requests.get(url, headers={"User-Agent": "moluq-e420-overflow/1.0"}, timeout=60)
    h = hashlib.sha256(resp.content).hexdigest()
    raw_path = os.path.join(ROOT, "data/e420/labels_raw", f"overflow_{entry}_{h[:12]}.cif")
    with open(raw_path, "wb") as f: f.write(resp.content)
    return resp, raw_path, url

entries = ov["entries"]
print("overflow entries to check:", len(entries), "| accessions:", len(ov["overflow_accessions"]))
temporal_pass, temporal_fail, no_date, fetch_fail = [], [], [], []
for i, entry in enumerate(sorted(entries)):
    try:
        rr, raw_path, url = get_cif(entry)
    except Exception as e:
        fetch_fail.append({"entry": entry, "error": str(e)[:120]}); continue
    if rr is None: break
    if rr.status_code != 200:
        fetch_fail.append({"entry": entry, "status": rr.status_code}); continue
    try:
        outcome = gv.validate_reference_cif(rr.content)
        if not outcome.ok:
            no_date.append({"entry": entry, "reason": outcome.failure.code.value}); continue
        date = outcome.payload.get("initial_release_date")
    except Exception as e:
        no_date.append({"entry": entry, "error": str(e)[:120]}); continue
    accs = [a for a, es in ov["acc_entries"].items() if entry in es]
    for acc in accs:
        mcd_val = mcd.get(acc)
        if not mcd_val:
            no_date.append({"entry": entry, "accession": acc, "reason": "no modelCreatedDate in census bytes"}); continue
        created = (mcd_val or "")[:10]
        passed = bool(date and created < date)  # strict: created < release
        (temporal_pass if passed else temporal_fail).append(
            {"accession": acc, "entry": entry, "model_created": created, "release": date})
    if (i + 1) % 25 == 0:
        print(f"  {i+1}/{len(entries)} entries checked | pass: {len(temporal_pass)} | fail: {len(temporal_fail)}")

summary = {
  "schema": "e420-overflow-temporal-check-v1",
  "overflow_accessions": len(ov["overflow_accessions"]),
  "entries_checked": len(entries), "requests": count["n"],
  "temporal_pass_accessions": len({r['accession'] for r in temporal_pass}),
  "temporal_fail_accessions": len({r['accession'] for r in temporal_fail}),
  "no_date_rows": len(no_date), "fetch_failures": len(fetch_fail),
  "temporal_pass": temporal_pass, "temporal_fail": temporal_fail, "no_date": no_date, "fetch_failures": fetch_fail,
}
with open(OUT, "w") as f: json.dump(summary, f, indent=1, sort_keys=True)
print(json.dumps({k: v for k, v in summary.items() if not isinstance(v, list)}, indent=1))
