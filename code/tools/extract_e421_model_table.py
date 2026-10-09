"""Write the derived AFDB model table of the post-2022 (e421) census relabel.

The relabel (code/experiments/e427_e421_relabel.py) downloaded the AFDB model
files of the post-2022 census again in October 2026; they are not
redistributed. This tool records, for every AFDB model file that its receipts
show as downloaded (HTTP 200), the URL, SHA-256 (verified against the
receipt), _ma_model_list.model_group_name, the _software name/version rows,
the census accession(s) whose qualifying entries point to the file and the
model's _ma_target_ref_db_details.db_accession, so that code/paper_numbers.py
can count model types (AlphaFold Monomer v2.0 vs ColabFold) offline. No
coordinate or pLDDT value is read.

Usage (needs the files, e.g. after
`python -I data/fetch/fetch_raw.py --set e421relabel --work _work`):
  python -I code/tools/extract_e421_model_table.py <derived_root> <raw_dir> <out_json>
  e.g. data/derived _work/data/e427/e421_raw
       data/derived/results/e427/e421_relabel/afdb_model_table.json
"""
import hashlib
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_v5_template_table import categories  # noqa: E402  (mmCIF reader written for these tools)

AFDB = re.compile(r"^https://alphafold\.ebi\.ac\.uk/files/AF-[^/]+\.cif$")  # UniProt-keyed and ColabFold (numeric) ids


def main(derived, raw_dir, out_path):
    receipts = json.loads((Path(derived) / "results/e427/e421_relabel/receipts.json").read_text())
    qual = json.loads((Path(derived) / "results/e421/qualifying_entries.json").read_text())
    acc_of = {}
    for q in qual:
        acc_of.setdefault(q["cif_url"], set()).add(q["afdb_accession"])
    files = {}
    for r in receipts:
        m = AFDB.match(r["url"])
        if not m or r.get("status") != 200:
            continue
        raw = (Path(raw_dir) / r["sha256"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != r["sha256"]:
            raise SystemExit(f"SHA-256 mismatch for {r['url']}")
        if r["url"] in files:
            if files[r["url"]]["sha256"] != r["sha256"]:
                raise SystemExit(f"two different files for {r['url']}")
            continue
        c = categories(raw.decode("utf-8"), want=("_ma_model_list", "_software", "_ma_target_ref_db_details"))
        files[r["url"]] = {"accessions_in_census": sorted(acc_of.get(r["url"], [])), "sha256": r["sha256"],
                           "target_db_accession": (c["_ma_target_ref_db_details"] or [{}])[0].get("db_accession"),
                           "model_group_name": [x.get("model_group_name") for x in c["_ma_model_list"]],
                           "software": [{"name": s.get("name"), "version": s.get("version")} for s in c["_software"]]}
    res = {"schema": "e421-afdb-model-table-v1",
           "source": "AFDB model files downloaded by the post-2022 census relabel (October 2026), selected via "
                     "results/e427/e421_relabel/receipts.json; AlphaFold DB data, CC BY 4.0",
           "n_files": len(files), "files": dict(sorted(files.items()))}
    Path(out_path).write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(f"{len(files)} model files -> {out_path}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
