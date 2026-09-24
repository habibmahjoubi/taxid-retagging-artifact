#!/usr/bin/env python3
"""Post-publication check (Corrigendum, 2026-09-24): were the sequences that have no merge
predecessor really absent before the transition?

The merge-aware audit (m1_quantify_merge_fix.py / snapshot_aware_validation.py) only sees
identifier changes recorded in merged.dmp. It cannot see a second mechanism: a sequence being
reassigned to a child node created *later* beneath a taxon that already existed (widespread in
2024, when NCBI separated virus names from binomial species names). This script checks, for
each accession/current-taxid pair given:

  1. the current node's rank, parent and creation date (NCBI E-utilities efetch);
  2. whether the parent node already existed before a cutoff date (default: end of 2023, the
     MSL39 snapshot);
  3. whether the genome was already in a real, contemporary reference -- the official Kraken2
     viral index of 9 October 2023 (library_report.tsv) -- and which node its minimizers were
     assigned to there (inspect.txt).

Usage: check_child_node_insertions.py [CANDIDATES_CSV] [CUTOFF_YYYY/MM/DD]
  CANDIDATES_CSV  columns: accession,taxid   (default: data/no_merge_predecessor_17.csv)
  CUTOFF          default: 2023/12/31

Standard library only; needs internet access (NCBI E-utilities, genome-idx.s3.amazonaws.com).
"""
import csv
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&retmode=xml&id="
K2_2023 = "https://genome-idx.s3.amazonaws.com/kraken/viral_20231009/"


def fetch(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode("utf-8", errors="replace")


def efetch_taxa(taxids):
    """Return {requested_taxid: dict(taxid, name, rank, created, parent)}; NCBI resolves a
    merged taxid to its current successor, so the returned taxid may differ."""
    out = {}
    ids = sorted(set(taxids))
    for i in range(0, len(ids), 100):
        root = ET.fromstring(fetch(EUTILS + ",".join(ids[i:i + 100])))
        for t in root.findall("Taxon"):
            rec = {"taxid": t.findtext("TaxId"), "name": t.findtext("ScientificName"),
                   "rank": t.findtext("Rank"), "created": (t.findtext("CreateDate") or "")[:10],
                   "parent": t.findtext("ParentTaxId"),
                   "names": [t.findtext("ScientificName") or ""] + [e.text or "" for e in (t.find("OtherNames") if t.find("OtherNames") is not None else [])
                                                                  if e.tag in ("EquivalentName", "Synonym", "GenbankCommonName", "CommonName")]}
            aka = [a.text for a in t.findall("AkaTaxIds/TaxId")]
            for key in [rec["taxid"]] + aka:
                out[key] = rec
        time.sleep(0.4)
    return out


def main():
    here = Path(__file__).resolve().parent.parent
    cand_path = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "data" / "no_merge_predecessor_17.csv"
    cutoff = sys.argv[2] if len(sys.argv) > 2 else "2023/12/31"
    rows = list(csv.DictReader(open(cand_path, encoding="utf-8")))

    library = fetch(K2_2023 + "library_report.tsv")
    inspect = {}
    for line in fetch(K2_2023 + "inspect.txt").splitlines():
        p = line.split("\t")
        if len(p) >= 6:
            inspect[p[4]] = (int(p[1]), int(p[2]), p[5].strip())  # clade, direct, name
    by_name = {}
    for tid, (_, _, name) in inspect.items():
        by_name.setdefault(name.lower(), []).append(tid)

    current = efetch_taxa([r["taxid"] for r in rows])
    parents = efetch_taxa([current[r["taxid"]]["parent"] for r in rows if r["taxid"] in current])

    # A node of the 2023 index bearing one of the current node's names, which NCBI now resolves to
    # the current taxid, is a predecessor merged after the study's merge-history snapshot.
    old_by_row = {}
    for r in rows:
        cur = current.get(r["taxid"])
        if cur:
            old = {t for n in cur["names"] for t in by_name.get(n.lower(), []) if t != cur["taxid"]}
            old_by_row[r["accession"]] = old
    resolved = efetch_taxa([t for old in old_by_row.values() for t in old]) if any(old_by_row.values()) else {}

    counts = {}
    print("accession\tcurrent_taxid\trank\tcreated\tparent\tparent_created\tin_2023_index\tcategory")
    for r in rows:
        cur = current.get(r["taxid"])
        par = parents.get(cur["parent"]) if cur else None
        in_lib = r["accession"].split(".")[0] in library
        new_node = cur is not None and cur["created"] > cutoff
        parent_old = par is not None and par["created"] <= cutoff
        if not in_lib:
            cat = "genome absent from 2023 reference (genuine composition change)"
        elif cur is not None and any(resolved.get(t, {}).get("taxid") == cur["taxid"]
                                     for t in old_by_row.get(r["accession"], ())):
            old = ",".join(t for t in old_by_row[r["accession"]] if resolved.get(t, {}).get("taxid") == cur["taxid"])
            cat = f"predecessor {old} (in 2023 index) merged into current taxid after the merge-history snapshot"
        elif new_node and parent_old:
            cat = "child node created after cutoff beneath a pre-existing node"
        else:
            cat = "present in 2023 reference (inspect manually)"
        counts[cat] = counts.get(cat, 0) + 1
        print("\t".join([r["accession"], r["taxid"], cur["rank"] if cur else "?",
                         cur["created"] if cur else "?",
                         f"{par['name']} ({par['taxid']})" if par else "?",
                         par["created"] if par else "?", "yes" if in_lib else "no", cat]))
    print("\nSummary:")
    for cat, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {n:3d}  {cat}")
    print("\nNote: a sequence whose current taxid is a node created after the cutoff, beneath a"
          " parent that already existed, was represented before the transition under that parent"
          " (confirm with the 2023 inspect.txt node assignment); it is apparent, not genuine, novelty.")


if __name__ == "__main__":
    main()
