"""Write the derived AFDB template table used by code/paper_numbers.py (Table 2, templates).

V5 (code/verification/v5_templates.py) reads the 138 AFDB model files of the
preregistered census roster from results/e422/raw. Those files (38.6 MB) are
not redistributed; this tool extracts from them only what the template
classification needs, so that paper_numbers.py can recompute the template
classes offline from this table plus the shipped PDBe/RCSB lookup cache
(data/derived/paper/arxiv/verification/cache/v5/http):

  per accession: the AFDB file's receipt URL, raw name and SHA-256 (verified
  against the census receipts), _ma_target_ref_db_details.db_accession,
  _ma_model_list.model_group_name, _software name/version, and every
  _ma_template_ref_db_details row with db_name = PDB together with the
  template_auth_asym_id from _ma_template_details.

No pLDDT (_ma_qa_metric_*) or coordinate is read. The mmCIF reader below is
written for this tool (it does not import V5's parser), so the table also
serves as an independent check of V5's per-accession template lists.

Usage (needs the AFDB model files, e.g. after
`python -I data/fetch/fetch_raw.py --set e422 --work _work`):
  python -I code/tools/extract_v5_template_table.py <work_root> <out_json>
The released table is
data/derived/paper/arxiv/verification/v5_afdb_template_table.json.
"""
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

AFDB = re.compile(r"^https://alphafold\.ebi\.ac\.uk/files/AF-(.+?)-F(\d+)-model_v(\d+)\.cif$")
WANT = ("_ma_target_ref_db_details", "_ma_model_list", "_software",
        "_ma_template_details", "_ma_template_ref_db_details")
# CIF 1.1 value tokens: a quoted string ends at a matching quote followed by whitespace/end
TOKEN = re.compile(r"""'(?:[^']|'(?=\S))*'(?=\s|$)|"(?:[^"]|"(?=\S))*"(?=\s|$)|\S+""")


def tokens(text):
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith(";"):  # multi-line text field
            buf = [ln[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            yield ("value", "\n".join(buf).strip())
            i += 1
            continue
        for t in TOKEN.findall(ln):
            if t.startswith("#"):
                break  # comment to end of line
            if len(t) >= 2 and t[0] in "'\"" and t[-1] == t[0]:
                yield ("value", t[1:-1])
            else:
                yield ("word", t)
        i += 1


def categories(text, want=WANT):
    """{category: [row dict]} for the requested categories (loops and key-value pairs)."""
    out = defaultdict(list)
    toks = list(tokens(text))
    j = 0
    while j < len(toks):
        kind, t = toks[j]
        if kind == "word" and t == "loop_":
            j += 1
            names = []
            while j < len(toks) and toks[j][0] == "word" and toks[j][1].startswith("_"):
                names.append(toks[j][1])
                j += 1
            vals = []
            while j < len(toks) and not (toks[j][0] == "word" and (toks[j][1].startswith("_")
                                                                 or toks[j][1] in ("loop_",)
                                                                 or toks[j][1].startswith("data_"))):
                vals.append(toks[j][1])
                j += 1
            cat = names[0].split(".")[0]
            if cat in want:
                if len(vals) % len(names):
                    raise ValueError(f"{cat}: {len(vals)} values for {len(names)} items")
                for k in range(0, len(vals), len(names)):
                    out[cat].append({n.split(".", 1)[1]: vals[k + m] for m, n in enumerate(names)})
            continue
        if kind == "word" and t.startswith("_") and "." in t:
            cat, item = t.split(".", 1)
            if cat in want:
                if not out[cat]:
                    out[cat].append({})
                out[cat][0][item] = toks[j + 1][1]
            j += 2
            continue
        j += 1
    return out


def main(work, out_path):
    work = Path(work)
    roster, recs = set(), defaultdict(dict)
    for b in (1, 2, 3, 4):
        d = json.loads((work / f"results/e422/batches/e422_batch0{b}_checkpoint.json").read_text())
        roster |= {r["accession"] for r in d["rows"]}
        for r in d["receipts"]:
            m = AFDB.match(r["url"])
            if m:
                recs[m.group(1)][r["raw_name"]] = r
    table = {}
    for acc in sorted(roster):
        names = sorted(recs[acc])
        if len(names) != 1:
            raise SystemExit(f"{acc}: expected one AFDB model file, found {len(names)}")
        rec = recs[acc][names[0]]
        raw = (work / "results/e422/raw" / names[0]).read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != rec["sha256"] or digest != names[0]:
            raise SystemExit(f"{acc}: SHA-256 mismatch for {rec['url']}")
        c = categories(raw.decode("utf-8"))
        chain_of = {t.get("template_id"): t.get("template_auth_asym_id") for t in c["_ma_template_details"]}
        table[acc] = {
            "url": rec["url"], "raw_name": names[0], "sha256": digest,
            "target_db_accession": (c["_ma_target_ref_db_details"] or [{}])[0].get("db_accession"),
            "model_group_name": [m.get("model_group_name") for m in c["_ma_model_list"]],
            "software": [{"name": s.get("name"), "version": s.get("version")} for s in c["_software"]],
            "templates": [{"template_id": t.get("template_id"), "pdb_id": t["db_accession_code"].upper(),
                           "template_auth_asym_id": chain_of.get(t.get("template_id"))}
                          for t in c["_ma_template_ref_db_details"] if t.get("db_name") == "PDB"],
            "non_pdb_template_refs": sum(t.get("db_name") != "PDB" for t in c["_ma_template_ref_db_details"]),
        }
    res = {"schema": "v5-afdb-template-table-v1",
           "source": "AFDB model files of the e422 census roster (results/e422/raw), selected via the batch 1-4 "
                     "checkpoint receipts; AlphaFold DB data, CC BY 4.0",
           "n_accessions": len(table), "accessions": table}
    Path(out_path).write_text(json.dumps(res, indent=1, sort_keys=True) + "\n")
    print(f"{len(table)} accessions, {sum(len(v['templates']) for v in table.values())} PDB templates -> {out_path}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
