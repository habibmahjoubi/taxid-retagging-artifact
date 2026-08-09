#!/usr/bin/env python3
"""Filter a kraken:taxid-tagged FASTA to only sequences whose taxid exists in a given
historical NCBI nodes.dmp snapshot. This operationalizes "knowledge drift" for the taxonomy
case study: a viral genome/species not yet present in an older taxonomy tree is excluded
from that year's Kraken2 library, exactly as it would have been absent from any reference
database built at that time. This is the naive tagging mechanism this study's article
identifies as producing false novelty when the sequence's taxid was later renumbered.

Usage: filter_fasta_by_taxid_existence.py IN_FASTA NODES_DMP OUT_FASTA
"""
import re
import sys

def load_valid_taxids(nodes_dmp_path):
    valid = set()
    with open(nodes_dmp_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            taxid = line.split("\t|\t", 1)[0].strip()
            if taxid.isdigit():
                valid.add(int(taxid))
    return valid

HEADER_TAXID_RE = re.compile(r"\|kraken:taxid\|(\d+)")

def main():
    in_fasta, nodes_dmp, out_fasta = sys.argv[1], sys.argv[2], sys.argv[3]
    valid_taxids = load_valid_taxids(nodes_dmp)
    print(f"Loaded {len(valid_taxids)} valid taxids from {nodes_dmp}")

    n_kept = 0
    n_dropped = 0
    dropped_taxids = set()
    with open(in_fasta, "r", encoding="utf-8", errors="replace") as fin, \
         open(out_fasta, "w", encoding="utf-8") as fout:
        write_seq = False
        for line in fin:
            if line.startswith(">"):
                m = HEADER_TAXID_RE.search(line)
                taxid = int(m.group(1)) if m else None
                if taxid is not None and taxid in valid_taxids:
                    write_seq = True
                    n_kept += 1
                    fout.write(line)
                else:
                    write_seq = False
                    n_dropped += 1
                    if taxid is not None:
                        dropped_taxids.add(taxid)
            elif write_seq:
                fout.write(line)

    print(f"Kept {n_kept} sequences, dropped {n_dropped} (taxid not yet present in this "
          f"snapshot) -- {len(dropped_taxids)} distinct dropped taxids")

if __name__ == "__main__":
    main()
