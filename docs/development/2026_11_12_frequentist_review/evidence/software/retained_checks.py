"""Read-only recorded identity/checksum verification; no sampler imports.

python -B retained_checks.py SOURCE CACHE
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("repo", type=Path)
ap.add_argument("cache", type=Path)
args = ap.parse_args()
base = args.repo / "docs/development/2026_11_10_early_toxicity"
ids = {}
for p in (base / "step5_review").glob("*/*/manifest.json"):
    ids[p.parent.relative_to(base).as_posix()] = json.loads(p.read_text())["identity"][:16]
cache_manifest = json.loads((args.cache / "cache_manifest.json").read_text())
cache_matches = []
for p in args.cache.glob("*.npz"):
    with np.load(p, allow_pickle=False) as z:
        identity = str(z["identity"])
    ids["step6_review/cache/" + p.name] = identity[:16]
    entry = cache_manifest[p.stem]
    cache_matches.append(entry["identity"] == identity and entry["sha256"] == hashlib.sha256(p.read_bytes()).hexdigest())
records = list(csv.DictReader((base / "chain_recheck_2026_11_11.csv").open()))
same = all(ids[r["fit"]] == r["identity"] for r in records)
print("recheck identities match", len(records), same)
print("cache manifest entries", len(cache_matches), "all file identities and bytes match", all(cache_matches))
report_root = base / "step6_review"
report = json.loads((report_root / "report_manifest.json").read_text())["files"]
report_ok = all(hashlib.sha256((report_root / p).read_bytes()).hexdigest() == h for p, h in report.items())
print("published report hashes match", len(report), report_ok)
assert same and all(cache_matches) and report_ok
