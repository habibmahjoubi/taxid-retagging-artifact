#!/usr/bin/env python3
"""Map every accession in the bulk RefSeq viral release FASTA to its NCBI taxid via batched
NCBI eutils esummary calls (fast: ~100 requests total instead of one FTP fetch per genome),
then rewrite the FASTA headers with the |kraken:taxid|TAXID tag Kraken2 expects for a custom
--add-to-library file, so kraken2-build --build can proceed without needing the huge generic
nucl_gb/nucl_wgs accession2taxid maps.
"""
import json
import os
import time
import urllib.request
import urllib.parse

BASE = os.environ.get("KDCR_BASE", ".")
IN_FASTA = f"{BASE}/kraken2_dbs/current/bulk/viral.1.1.genomic.fna"
OUT_FASTA = f"{BASE}/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna"
MAP_TSV = f"{BASE}/kraken2_dbs/current/bulk/accession_taxid_map.tsv"
BATCH = 200
SLEEP = 0.4

def load_accessions(path):
    accs = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                accs.append(line[1:].split()[0])
    return accs

def esummary_batch(accs):
    ids = ",".join(accs)
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "nuccore", "id": ids, "retmode": "json"}
    )
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.load(resp)
    out = {}
    result = data.get("result", {})
    for uid in result.get("uids", []):
        entry = result.get(uid, {})
        acc = entry.get("accessionversion") or entry.get("caption")
        taxid = entry.get("taxid")
        if acc and taxid:
            out[acc] = taxid
    return out

def main():
    accs = load_accessions(IN_FASTA)
    unique = sorted(set(accs))
    print(f"Total sequences: {len(accs)}, unique accessions: {len(unique)}")

    acc2taxid = {}
    if os.path.exists(MAP_TSV):
        with open(MAP_TSV, "r", encoding="utf-8") as f:
            for line in f:
                a, t = line.rstrip("\n").split("\t")
                acc2taxid[a] = int(t)
        print(f"Loaded {len(acc2taxid)} cached accession->taxid pairs from {MAP_TSV}")

    remaining = [a for a in unique if a not in acc2taxid]
    unique_to_query = remaining
    for i in range(0, len(unique_to_query), BATCH):
        batch = unique_to_query[i:i + BATCH]
        for attempt in range(3):
            try:
                acc2taxid.update(esummary_batch(batch))
                break
            except Exception as e:
                print(f"  batch {i}: retry after error {e}")
                time.sleep(2)
        print(f"  mapped {len(acc2taxid)}/{len(unique)} so far (batch starting at {i})")
        time.sleep(SLEEP)

    missing = [a for a in unique if a not in acc2taxid]
    print(f"Resolved {len(acc2taxid)}/{len(unique)} accessions; {len(missing)} missing")

    with open(MAP_TSV, "w", encoding="utf-8") as f:
        for acc, taxid in sorted(acc2taxid.items()):
            f.write(f"{acc}\t{taxid}\n")

    n_written = 0
    n_skipped = 0
    with open(IN_FASTA, "r", encoding="utf-8", errors="replace") as fin, \
         open(OUT_FASTA, "w", encoding="utf-8") as fout:
        write_seq = False
        for line in fin:
            if line.startswith(">"):
                rest = line[1:].rstrip("\n")
                parts = rest.split(None, 1)
                acc = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                taxid = acc2taxid.get(acc)
                if taxid is None:
                    write_seq = False
                    n_skipped += 1
                    continue
                write_seq = True
                n_written += 1
                # kraken:taxid tag must directly follow the sequence ID token, before the
                # free-text description, or kraken2-build's header parser will not detect it
                # and silently falls back to ACCNUM-based resolution (which then requires the
                # accession2taxid maps we deliberately skipped).
                fout.write(f">{acc}|kraken:taxid|{taxid} {description}\n")
            elif write_seq:
                fout.write(line)

    print(f"Wrote {n_written} tagged sequences, skipped {n_skipped} (no taxid resolved)")


if __name__ == "__main__":
    main()
