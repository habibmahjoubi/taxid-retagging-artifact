#!/usr/bin/env python3
"""Snapshot-specific validation of the merge-aware audit (Limitations, recommendation 2).

The library-wide audit behind Table 2 (m1_quantify_merge_fix.py) resolves each identifier's
predecessor via the single most recent merged.dmp (MSL41), cumulative but blind to any merge
that was later reversed (e.g. taxid 35320, merged into 3052452 by MSL39 then reinstated by
MSL41 -- invisible in a cumulative file compiled after the reversal). This script instead
checks each snapshot S against its OWN historical record, using two directions:

  (1) forward / self-merge: is the identifier itself a merge source in merged.dmp(S)? If so
      it was retired by S into some current-at-S target -- this is what recovers reversal
      cases like taxid 35320's, which direction (2) alone cannot see.
  (2) reverse / predecessor: does merged.dmp(S_j) for some snapshot S_j after S record a
      predecessor of the identifier that was itself valid at S? This is the standard
      renumbering case (e.g. NC_034381.1, taxid 1195365 -> 2870378).

Run on a freshly-drawn sample of current NCBI RefSeq viral genome taxids (not the exact
15,359-identifier tagged library behind Table 2, which this repository does not redistribute)
against all thirteen snapshots, this recovers 2,488 of 47,207 naively-absent checks versus
2,446 under the cumulative-only method -- a 1.72% relative increase, 0.089% discordance rate
(95% CI 0.066-0.120%, n=47,207). It independently reproduces the taxid 35320 case and resolves
Table 2's own noted anomaly for taxid 2651918 (Xanthomonas phage XaF13).

Expects, under $KDCR_BASE:
  extracted/<msl>/nodes.dmp    per-snapshot taxonomy nodes, msl29..msl41
  extracted/<msl>/merged.dmp   per-snapshot merge table, msl29..msl41
  all_taxids.txt               one taxid per line, the identifiers to audit

Reproduce the sampling frame used in the manuscript with:
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
import math
import os

BASE = os.environ.get("KDCR_BASE", ".")
MSLS = [f"msl{n}" for n in range(29, 42)]


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

    # OLD (cumulative-only) method, matching m1_quantify_merge_fix.py: reverse map from the
    # single latest merged.dmp (MSL41) only.
    reverse_old = {}
    for old, new in merged_maps["msl41"].items():
        reverse_old.setdefault(new, set()).add(old)

    # NEW method direction (2): reverse map from every snapshot strictly AFTER S.
    reverse_new = {}
    for i, S in enumerate(order):
        rev = {}
        for j in range(i + 1, len(order)):
            for old, new in merged_maps[order[j]].items():
                rev.setdefault(new, set()).add(old)
        reverse_new[S] = rev

    with open(f"{BASE}/all_taxids.txt") as f:
        taxids = [int(x) for x in f.read().split()]

    def old_recovers(t, S):
        return bool(reverse_old.get(t, set()) & valid_sets[S])

    def new_recovers(t, S):
        target = merged_maps[S].get(t)  # direction (1): self as merge source at S
        if target is not None and target in valid_sets[S]:
            return True
        return bool(reverse_new[S].get(t, set()) & valid_sets[S])  # direction (2)

    tot_na = tot_old = tot_new = tot_disc = 0
    discordant = []
    for S in order:
        valid_S = valid_sets[S]
        naive_absent = [t for t in taxids if t not in valid_S]
        old_r = {t for t in naive_absent if old_recovers(t, S)}
        new_r = {t for t in naive_absent if new_recovers(t, S)}
        disc = new_r - old_r
        for t in disc:
            discordant.append((S, t, merged_maps[S].get(t)))
        print(f"{S}: naive_absent={len(naive_absent)} old={len(old_r)} new={len(new_r)} "
              f"discordant={len(disc)}")
        tot_na += len(naive_absent)
        tot_old += len(old_r)
        tot_new += len(new_r)
        tot_disc += len(disc)

    print(f"\nTOTAL naive_absent={tot_na} old_recovered={tot_old} new_recovered={tot_new} "
          f"discordant={tot_disc}")
    if tot_old:
        print(f"relative increase: {100 * tot_disc / tot_old:.2f}%")
    p = tot_disc / tot_na
    z = 1.96
    denom = 1 + z * z / tot_na
    center = (p + z * z / (2 * tot_na)) / denom
    halfwidth = z * math.sqrt(p * (1 - p) / tot_na + z * z / (4 * tot_na * tot_na)) / denom
    print(f"discordance rate: {100 * p:.3f}% (95% CI {100 * (center - halfwidth):.3f}-"
          f"{100 * (center + halfwidth):.3f}%, n={tot_na})")
    print(f"\ndistinct discordant taxids: {sorted(set(t for _, t, _ in discordant))}")


if __name__ == "__main__":
    main()
