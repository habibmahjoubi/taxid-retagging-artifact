#!/usr/bin/env python3
"""Snapshot-specific merge-aware audit -- generates the manuscript's Table 2 (recommendation 2).

An earlier version of this analysis (m1_quantify_merge_fix.py) resolved each identifier's
predecessor via the single most recent merged.dmp (MSL41), cumulative but blind to any merge
that was later reversed (e.g. taxid 35320, merged into 3052452 by MSL39 then reinstated by
MSL41 -- invisible in a cumulative file compiled after the reversal). This script instead
checks each snapshot S against its OWN historical record, using two directions:

  (1) forward / self-merge: is the identifier itself a merge source in merged.dmp(S), with a
      valid target there? If so it was retired by S into some then-current target -- this is
      what recovers reversal cases like taxid 35320's, and single-snapshot cases like taxid
      2651918's (Xanthomonas phage XaF13, naively absent even at the most recent snapshot),
      neither of which direction (2) alone can see.
  (2) reverse / predecessor: does merged.dmp(S_j) for some snapshot S_j after S record a
      predecessor of the identifier that was itself valid at S? This is the standard
      renumbering case (e.g. NC_034381.1, taxid 1195365 -> 2870378).

Run on an independently obtained population of 14,607 current NCBI RefSeq viral genome taxids
(not the same 19,582 sequences tagged for the 56-candidate panel and Kraken2 experiment, which
this repository does not redistribute) against all thirteen snapshots, this is the manuscript's
Table 2: naive valid / corrected valid / recovered per snapshot, and the resulting 13-way
fixed-library intersection (5,231 naive -> 5,481 corrected, +4.8%). Comparability between this
population and the original panel's, under the cumulative-only method used for both, is also
reported (2,446 vs 2,455 recovered; 24.9% vs 26.8% of checks naively absent; 5.18% vs 4.59%
recovery rate among them).

Expects, under $KDCR_BASE:
  extracted/<msl>/nodes.dmp    per-snapshot taxonomy nodes, msl29..msl41
  extracted/<msl>/merged.dmp   per-snapshot merge table, msl29..msl41
  all_taxids.txt               one taxid per line, the identifiers to audit

Reproduce the population used in the manuscript with:
  curl -sL -o assembly_summary.txt \\
    https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
  awk -F'\\t' 'NR>2 {print $6}' assembly_summary.txt | sort -n | uniq > all_taxids.txt

Reproduce the thirteen per-snapshot taxonomy files (not redistributed here given their
combined size, and because NCBI already archives them permanently at a stable URL) with:
  declare -A MSL_DATE=(
    [msl29]=2014-12-01 [msl30]=2015-12-01 [msl31]=2016-12-01 [msl32]=2017-12-01
    [msl33]=2018-06-01 [msl34]=2018-12-01 [msl35]=2019-12-01 [msl36]=2020-12-01
    [msl37]=2021-12-01 [msl38]=2022-12-01 [msl39]=2023-12-01 [msl40]=2024-12-01
    [msl41]=2025-12-01
  )
  for msl in "${!MSL_DATE[@]}"; do
    date="${MSL_DATE[$msl]}"
    curl -sL -o "taxdmp_${date}.zip" \\
      "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_${date}.zip"
    mkdir -p "extracted/$msl"
    unzip -o -q "taxdmp_${date}.zip" nodes.dmp merged.dmp -d "extracted/$msl"
  done
"""
import os

BASE = os.environ.get("KDCR_BASE", ".")
MSLS = [f"msl{n}" for n in range(29, 42)]
MSL_LABEL = {
    "msl29": "MSL29 (2014)", "msl30": "MSL30 (2015)", "msl31": "MSL31 (2016)",
    "msl32": "MSL32 (2017)", "msl33": "MSL33 (2018a)", "msl34": "MSL34 (2018b)",
    "msl35": "MSL35 (2019)", "msl36": "MSL36 (2020)", "msl37": "MSL37 (2021)",
    "msl38": "MSL38 (2022)", "msl39": "MSL39 (2023)", "msl40": "MSL40 (2024)",
    "msl41": "MSL41 (2025)",
}

# Cumulative-only totals reported for the original 15,359-identifier tagged library
# (the panel's own library, via m1_quantify_merge_fix.py), used only for the
# comparability check below -- not recomputed here since that exact library is not
# redistributed in this repository.
ORIGINAL_LIBRARY_SIZE = 15359
ORIGINAL_CUMULATIVE_RECOVERED = 2455
ORIGINAL_NAIVE_ABSENT_CHECKS = 53504


def load_valid(msl):
    valid = set()
    with open(f"{BASE}/extracted/{msl}/nodes.dmp", encoding="utf-8", errors="replace") as f:
        for line in f:
            t = line.split("\t|\t", 1)[0].strip()
            if t.isdigit():
                valid.add(int(t))
    return valid


def load_merged(msl):
    m = {}
    with open(f"{BASE}/extracted/{msl}/merged.dmp", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split("\t|\t")
            if len(parts) >= 2:
                old = parts[0].strip()
                new = parts[1].split("\t")[0].strip().rstrip("|").strip()
                if old.isdigit() and new.isdigit():
                    m[int(old)] = int(new)
    return m


def main():
    valid_sets = {msl: load_valid(msl) for msl in MSLS}
    merged_maps = {msl: load_merged(msl) for msl in MSLS}
    order = MSLS

    # Cumulative-only reverse map, matching m1_quantify_merge_fix.py: from the single
    # latest merged.dmp (MSL41) only. Used for the comparability check, not for Table 2.
    reverse_cumulative = {}
    for old, new in merged_maps["msl41"].items():
        reverse_cumulative.setdefault(new, set()).add(old)

    # Snapshot-specific reverse map, direction (2): from every snapshot strictly AFTER S.
    reverse_snapshot = {}
    for i, S in enumerate(order):
        rev = {}
        for j in range(i + 1, len(order)):
            for old, new in merged_maps[order[j]].items():
                rev.setdefault(new, set()).add(old)
        reverse_snapshot[S] = rev

    with open(f"{BASE}/all_taxids.txt") as f:
        taxids = [int(x) for x in f.read().split()]
    print(f"Population: {len(taxids)} distinct current viral RefSeq taxids\n")

    def cumulative_recovers(t, S):
        return bool(reverse_cumulative.get(t, set()) & valid_sets[S])

    def snapshot_recovers(t, S):
        target = merged_maps[S].get(t)  # direction (1): self as merge source at S
        if target is not None and target in valid_sets[S]:
            return True
        return bool(reverse_snapshot[S].get(t, set()) & valid_sets[S])  # direction (2)

    naive_sets, corrected_sets = {}, {}
    tot_naive_absent = tot_cumulative_recovered = 0
    print(f"{'Snapshot':16s} {'Naive valid':12s} {'Corrected valid':16s} {'Recovered':10s}")
    for S in order:
        valid_S = valid_sets[S]
        naive = {t for t in taxids if t in valid_S}
        naive_absent = set(taxids) - naive
        corrected = set(naive)
        cumulative_recovered_here = 0
        for t in naive_absent:
            if snapshot_recovers(t, S):
                corrected.add(t)
            if cumulative_recovers(t, S):
                cumulative_recovered_here += 1
        naive_sets[S], corrected_sets[S] = naive, corrected
        print(f"{MSL_LABEL[S]:16s} {len(naive):12d} {len(corrected):16d} "
              f"{len(corrected) - len(naive):10d}")
        tot_naive_absent += len(naive_absent)
        tot_cumulative_recovered += cumulative_recovered_here

    naive_13way = set(taxids)
    corrected_13way = set(taxids)
    for S in order:
        naive_13way &= naive_sets[S]
        corrected_13way &= corrected_sets[S]
    print(f"\n13-way naive intersection: {len(naive_13way)}")
    print(f"13-way corrected intersection: {len(corrected_13way)} "
          f"(+{len(corrected_13way) - len(naive_13way)}, "
          f"{100 * (len(corrected_13way) - len(naive_13way)) / len(naive_13way):.2f}%)")

    print("\n--- Comparability with the original 15,359-identifier tagged library "
          "(cumulative-only method, both populations) ---")
    share_here = 100 * tot_naive_absent / (len(taxids) * len(order))
    rate_here = 100 * tot_cumulative_recovered / tot_naive_absent
    share_orig = 100 * ORIGINAL_NAIVE_ABSENT_CHECKS / (ORIGINAL_LIBRARY_SIZE * len(order))
    rate_orig = 100 * ORIGINAL_CUMULATIVE_RECOVERED / ORIGINAL_NAIVE_ABSENT_CHECKS
    print(f"this population:  {tot_cumulative_recovered} recovered, "
          f"{share_here:.1f}% of checks naively absent, {rate_here:.2f}% recovery rate")
    print(f"original library: {ORIGINAL_CUMULATIVE_RECOVERED} recovered, "
          f"{share_orig:.1f}% of checks naively absent, {rate_orig:.2f}% recovery rate")


if __name__ == "__main__":
    main()
