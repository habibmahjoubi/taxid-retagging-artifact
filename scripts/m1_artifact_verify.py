#!/usr/bin/env python3
"""Directly test the retagging hypothesis on the 56-taxon candidate panel: for the subset of
accessions whose CURRENT taxid did not exist in nodes.dmp(MSL39) but whose taxonomic concept
did, under a predecessor id later merged into the current one (per merged.dmp), does a
period-accurate library (tagging by the taxid actually valid at MSL39) classify these reads
correctly at MSL39 -- confirming the 'before' 0% in the naive panel is a tagging artifact, not
real absence of the taxonomic concept from the knowledge base? This produces the article's
central experimental result: 99.3% correct classification under period-accurate retagging
versus 0.0% under naive current-taxid tagging (Figure 1).

Builds a small supplementary Kraken2 DB: MSL39's existing filtered library + the 39 "artifact"
sequences, retagged with their MSL39-valid predecessor taxid instead of the current one.

Requires: kraken2 (tested with v2.17.1) on PATH; simulated paired-end reads per accession
(e.g. via wgsim) already generated under FASTA_DIR as <accession>_1.fastq / <accession>_2.fastq.

Expects, under $KDCR_BASE:
  results/positive_control_56taxa.csv   (accession, taxid, description columns; the 56-taxon panel)
  positive_control_56/<accession>.fna, <accession>_1.fastq, <accession>_2.fastq
  kraken2_dbs/msl39/viral.filtered.fna
  kraken2_dbs/msl39/taxonomy/nodes.dmp
  kraken2_dbs/msl41/taxonomy/merged.dmp
"""
import csv
import os
import re
import shutil
import subprocess

BASE = os.environ.get("KDCR_BASE", ".")
POS56 = f"{BASE}/results/positive_control_56taxa.csv"
FASTA_DIR = f"{BASE}/positive_control_56"
MSL39_FILTERED = f"{BASE}/kraken2_dbs/msl39/viral.filtered.fna"
MSL39_TAXONOMY = f"{BASE}/kraken2_dbs/msl39/taxonomy"
NEWDB = f"{BASE}/kraken2_dbs_msl39_artifact_test"


def load_valid_taxids(p):
    valid = set()
    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            t = line.split("\t|\t", 1)[0].strip()
            if t.isdigit():
                valid.add(int(t))
    return valid


def load_merged(p):
    m = {}
    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.split("\t|\t")
            if len(parts) >= 2:
                old = parts[0].strip()
                new = parts[1].split("\t")[0].strip().rstrip("|").strip()
                if old.isdigit() and new.isdigit():
                    m[int(old)] = int(new)
    return m


def main():
    merged = load_merged(f"{BASE}/kraken2_dbs/msl41/taxonomy/merged.dmp")
    reverse = {}
    for old, new in merged.items():
        reverse.setdefault(new, set()).add(old)
    nodes39 = load_valid_taxids(f"{MSL39_TAXONOMY}/nodes.dmp")

    rows = list(csv.DictReader(open(POS56)))
    artifact_rows = []
    for r in rows:
        t = int(r["taxid"])
        preds = reverse.get(t, set()) & nodes39
        if preds:
            pred = sorted(preds)[0]
            artifact_rows.append((r["accession"], t, pred, r["description"]))

    print(f"Artifact rows: {len(artifact_rows)}")

    os.makedirs(NEWDB, exist_ok=True)
    # 1. copy msl39's taxonomy (predecessor ids already valid there, no taxonomy change needed)
    if not os.path.exists(f"{NEWDB}/taxonomy"):
        shutil.copytree(MSL39_TAXONOMY, f"{NEWDB}/taxonomy")

    # 2. build combined filtered fasta: msl39's existing library + 39 artifact seqs retagged
    combined = f"{NEWDB}/combined.fna"
    with open(combined, "w") as fout:
        with open(MSL39_FILTERED, encoding="utf-8", errors="replace") as fin:
            fout.write(fin.read())
        for acc, cur_taxid, pred_taxid, desc in artifact_rows:
            fna = f"{FASTA_DIR}/{acc}.fna"
            with open(fna, encoding="utf-8", errors="replace") as fin:
                lines = fin.readlines()
            fout.write(f">{acc}|kraken:taxid|{pred_taxid} {desc}\n")
            for line in lines[1:]:
                fout.write(line)
    print(f"Combined library written: {combined}")

    # 3. build kraken2 db
    subprocess.run(["kraken2-build", "--add-to-library", combined, "--db", NEWDB], check=True)
    subprocess.run(["kraken2-build", "--build", "--db", NEWDB, "--threads", "8"], check=True)
    print("DB build complete.")

    # 4. classify each artifact accession's existing simulated reads, report correct-classification %
    out_csv = f"{BASE}/results/m1_artifact_verification.csv"
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as fout:
        w = csv.writer(fout)
        w.writerow(["accession", "current_taxid", "predecessor_taxid", "n_reads",
                    "n_correct_to_predecessor", "pct_correct_predecessor",
                    "pct_any_classified", "description"])
        for acc, cur_taxid, pred_taxid, desc in artifact_rows:
            r1 = f"{FASTA_DIR}/{acc}_1.fastq"
            r2 = f"{FASTA_DIR}/{acc}_2.fastq"
            out_report = f"{NEWDB}/{acc}_report.txt"
            out_output = f"{NEWDB}/{acc}_output.txt"
            subprocess.run(
                ["kraken2", "--db", NEWDB, "--paired", "--threads", "4",
                 "--report", out_report, "--output", out_output, r1, r2],
                check=True, capture_output=True,
            )
            n_total = 0
            n_correct = 0
            n_any = 0
            with open(out_output) as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    n_total += 1
                    status, taxid_field = parts[0], parts[2]
                    m = re.search(r"\(taxid (\d+)\)", taxid_field) or re.match(r"^(\d+)$", taxid_field)
                    assigned = int(m.group(1)) if m else None
                    if status != "U":
                        n_any += 1
                    if assigned == pred_taxid:
                        n_correct += 1
            pct_correct = 100.0 * n_correct / n_total if n_total else None
            pct_any = 100.0 * n_any / n_total if n_total else None
            print(f"{acc} (cur={cur_taxid} pred={pred_taxid}): n={n_total} "
                  f"correct_to_pred={n_correct} ({pct_correct:.2f}%) any_classified={pct_any:.2f}%")
            w.writerow([acc, cur_taxid, pred_taxid, n_total, n_correct,
                        f"{pct_correct:.2f}", f"{pct_any:.2f}", desc])
    print(f"\nSaved: {out_csv}")


if __name__ == "__main__":
    main()
