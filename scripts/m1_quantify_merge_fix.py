#!/usr/bin/env python3
"""Library-wide merge-aware audit: does treating a taxid as 'existing at snapshot T' only if
it (or a predecessor id that later merged into it) appears in nodes.dmp(T) change the
per-version library-size trajectory and the 13-way fixed-library intersection? This is the
central quantification behind the article's Table 1 (218-269 recovered identifiers per
snapshot; 13-way intersection growing from 5,270 to 5,518 once merge-corrected).

Uses the LATEST available merged.dmp (from the most recent snapshot, MSL41 here) as the
canonical, cumulative old->new taxid mapping (NCBI keeps this collapsed to the current valid
id, not chained hop-by-hop -- verified below by checking no old_id also appears as a new_id
in the same file).

Expects a directory layout of the form:
  $KDCR_BASE/kraken2_dbs/<msl>/taxonomy/nodes.dmp    (one per snapshot, msl29..msl41)
  $KDCR_BASE/kraken2_dbs/msl41/taxonomy/merged.dmp   (latest cumulative merge table)
  $KDCR_BASE/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna  (tagged bulk FASTA)
  $KDCR_BASE/kraken2_dbs_fixed/fixed_taxids.txt      (naive 13-way intersection, for comparison)
"""
import os
import re

BASE = os.environ.get("KDCR_BASE", ".")
MSLS = [f"msl{n}" for n in range(29, 42)]
TAGGED_FASTA = f"{BASE}/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna"
LATEST_MERGED = f"{BASE}/kraken2_dbs/msl41/taxonomy/merged.dmp"

HEADER_RE = re.compile(r"^>(\S+)\|kraken:taxid\|(\d+)")


def load_valid_taxids(nodes_dmp_path):
    valid = set()
    with open(nodes_dmp_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            taxid = line.split("\t|\t", 1)[0].strip()
            if taxid.isdigit():
                valid.add(int(taxid))
    return valid


def load_merged(merged_dmp_path):
    m = {}
    with open(merged_dmp_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split("\t|\t")
            if len(parts) >= 2:
                old = parts[0].strip()
                new = parts[1].split("\t")[0].strip().rstrip("|").strip()
                if old.isdigit() and new.isdigit():
                    m[int(old)] = int(new)
    return m


def main():
    tagged_taxids = set()
    with open(TAGGED_FASTA, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                mm = HEADER_RE.match(line.strip())
                if mm:
                    tagged_taxids.add(int(mm.group(2)))
    print(f"Tagged taxids in bulk FASTA: {len(tagged_taxids)}")

    merged = load_merged(LATEST_MERGED)
    print(f"merged.dmp (latest snapshot) entries: {len(merged)}")
    # Sanity check: is this already collapsed (no old_id also a new_id)?
    news = set(merged.values())
    olds = set(merged.keys())
    chained = olds & news
    print(f"IDs appearing as both old and new (chain check): {len(chained)} "
          f"({'collapsed/cumulative' if len(chained) == 0 else 'NOT fully collapsed -- needs transitive closure'})")

    # reverse map: current taxid -> set of predecessor ids that ever merged into it
    reverse = {}
    for old, new in merged.items():
        reverse.setdefault(new, set()).add(old)

    tagged_with_history = sum(1 for t in tagged_taxids if t in reverse)
    print(f"Tagged taxids that have >=1 predecessor id in merged.dmp: {tagged_with_history}")

    naive_sets = {}
    corrected_sets = {}
    recovered_examples = []
    for msl in MSLS:
        nodes_path = f"{BASE}/kraken2_dbs/{msl}/taxonomy/nodes.dmp"
        valid = load_valid_taxids(nodes_path)
        naive = {t for t in tagged_taxids if t in valid}
        corrected = set(naive)
        for t in tagged_taxids - naive:
            preds = reverse.get(t)
            if preds and (preds & valid):
                corrected.add(t)
                if len(recovered_examples) < 15:
                    recovered_examples.append((msl, t, sorted(preds & valid)[:3]))
        naive_sets[msl] = naive
        corrected_sets[msl] = corrected
        print(f"{msl}: naive={len(naive)}  corrected={len(corrected)}  recovered={len(corrected) - len(naive)}")

    print("\nExamples of recovered taxids (would have been wrongly excluded):")
    for msl, t, preds in recovered_examples:
        print(f"  {msl}: taxid {t} recovered via predecessor(s) {preds}")

    fixed_naive = set(tagged_taxids)
    fixed_corrected = set(tagged_taxids)
    for msl in MSLS:
        fixed_naive &= naive_sets[msl]
        fixed_corrected &= corrected_sets[msl]
    print(f"\n13-way fixed-library intersection: naive={len(fixed_naive)}  corrected={len(fixed_corrected)}")
    print(f"Newly-recovered stable taxids (in corrected fixed lib but not naive): {len(fixed_corrected - fixed_naive)}")

    # Does the corrected fixed library still match the naive intersection already reported
    # elsewhere for this project? Just report identity for transparency.
    fixed_taxids_file = f"{BASE}/kraken2_dbs_fixed/fixed_taxids.txt"
    if os.path.exists(fixed_taxids_file):
        with open(fixed_taxids_file) as f:
            original_fixed_file = {int(x) for x in f.read().split()}
        print(f"\noriginal_fixed_file ({fixed_taxids_file}) size: {len(original_fixed_file)}")
        print(f"matches naive recompute: {original_fixed_file == fixed_naive}")


if __name__ == "__main__":
    main()
