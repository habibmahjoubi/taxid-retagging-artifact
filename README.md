# taxid-retagging-artifact

Code and analysis scripts for:

> Mahjoubi H. **When Current Taxids Make Historical Taxa Look New.** Critical
> Comments, Journal of Bioinformatics and Computational Biology. Preprint:
> Zenodo, https://doi.org/10.5281/zenodo.21861802.

## Summary

Reconstructing a historical reference database by resolving each sequence's identifier via a
current-day lookup can make an already-existing taxon look newly created whenever the
identifier has been renumbered -- apparent novelty caused by identifier evolution, not
biological novelty. Among 56 viral sequences flagged as newly created between two NCBI
taxonomy snapshots (ICTV MSL39 -> MSL40), 39 existed at the earlier snapshot under a
predecessor taxid later merged into today's identifier. Retagging these 39 sequences with
their period-valid identifier recovers 99.3% correct Kraken2 classification, versus 0.0% under
naive current-taxid tagging. Extended library-wide, the same mechanism recovers 222-272
identifiers per snapshot for ten of thirteen dated snapshots studied (2014-2022).

This repository contains the scripts that produced that quantification (Table 2) and the
experimental verification (Figure 1). It does **not** contain the full pipeline for the
larger revalidation-framework study this note is drawn from (in preparation, cited above as
the companion manuscript) -- only the self-contained taxid-retagging analysis.

## Repository contents

```
scripts/
  tag_viral_fasta_with_taxid.py          Tag a bulk RefSeq viral FASTA with |kraken:taxid|N
                                          headers via NCBI E-utilities (esummary).
  filter_fasta_by_taxid_existence.py     Filter a tagged FASTA to taxids valid in a given
                                          historical nodes.dmp snapshot (the naive-tagging
                                          mechanism this study audits).
  find_newly_added_taxon.py              Identify candidate sequences whose taxid is present
                                          in a newer snapshot but absent from an older one
                                          (source of the 56-sequence candidate panel).
  m1_quantify_merge_fix.py               Cumulative-only merge-aware audit (predecessor via
                                          the single latest merged.dmp only). Superseded by
                                          snapshot_aware_validation.py for Table 2; kept for
                                          the comparability check reported there and in the
                                          manuscript's Analysis Protocol.
  m1_artifact_verify.py                  Builds a supplementary Kraken2 database with the 39
                                          artifact sequences retagged to their period-valid
                                          predecessor identifier and reclassifies their
                                          simulated reads -> Figure 1's 99.3%/0.0% result.
  paperA_figure.py                       Renders Figure 1 from the reported percentages.
  snapshot_aware_validation.py           Snapshot-specific merge-aware audit -> Table 2:
                                          checks each snapshot's own nodes.dmp/merged.dmp
                                          directly, in both directions, instead of only the
                                          latest cumulative merged.dmp -- recovers reversal
                                          cases (e.g. taxid 35320) and single-snapshot cases
                                          (e.g. taxid 2651918 at MSL41) that the
                                          cumulative-only method misses.
```

## Dependencies

- Python 3.8+ (standard library only for all scripts except `paperA_figure.py`)
- `matplotlib`, `numpy` (see `requirements.txt`) -- only for `paperA_figure.py`
- [Kraken2](https://github.com/DerrickWood/kraken2) (tested with v2.17.1) on `PATH`, for
  `m1_artifact_verify.py`
- [wgsim](https://github.com/lh3/wgsim) or equivalent, to simulate the paired-end reads
  consumed by `m1_artifact_verify.py` (not included here; any read simulator producing
  `<accession>_1.fastq` / `<accession>_2.fastq` per sequence works)
- Internet access to NCBI E-utilities, for `tag_viral_fasta_with_taxid.py`

Install the Python dependencies with:

```
pip install -r requirements.txt
```

## Data

This study uses exclusively public data, none of which is redistributed in this repository:

- NCBI taxonomy archive snapshots (`nodes.dmp`, `names.dmp`, `merged.dmp`, `delnodes.dmp` per
  release): https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/
- ICTV Master Species Lists MSL29-MSL41 (used to date each taxonomy snapshot to an ICTV
  release): https://ictv.global/msl
- NCBI RefSeq viral genome release (the bulk FASTA tagged and filtered by the scripts above):
  https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/

## Reproducing the results

All scripts read their inputs from a base directory given by the `KDCR_BASE` environment
variable (default: current directory), expected to contain one subdirectory per ICTV MSL
snapshot with that snapshot's `taxonomy/nodes.dmp` (and, for the latest snapshot,
`merged.dmp`), plus the tagged bulk FASTA. See each script's docstring for its exact expected
paths.

1. Download the taxonomy archive snapshots and the RefSeq viral release listed above.
2. Tag the bulk FASTA with current taxids:
   `KDCR_BASE=/path/to/data python scripts/tag_viral_fasta_with_taxid.py`
3. For each snapshot, filter the tagged FASTA to that snapshot's valid taxids:
   `python scripts/filter_fasta_by_taxid_existence.py IN_FASTA NODES_DMP OUT_FASTA`
4. Identify newly-added-taxid candidates between two snapshots:
   `python scripts/find_newly_added_taxon.py OLD_NODES_DMP NEW_NODES_DMP TAGGED_FASTA`
5. Run the cumulative-only merge-aware audit (comparability check only, not Table 2):
   `KDCR_BASE=/path/to/data python scripts/m1_quantify_merge_fix.py`
6. Run the experimental verification on the 56-sequence panel (Figure 1 data):
   `KDCR_BASE=/path/to/data python scripts/m1_artifact_verify.py`
7. Render Figure 1: `python scripts/paperA_figure.py`
8. Run the snapshot-specific merge-aware audit (Table 2): download each snapshot's own
   `nodes.dmp`/`merged.dmp` (not just the latest) and the independently obtained
   population of tagged taxids into `$KDCR_BASE/extracted/<msl>/{nodes.dmp,merged.dmp}`
   and `$KDCR_BASE/all_taxids.txt` (see the script's docstring for the exact commands),
   then: `KDCR_BASE=/path/to/data python scripts/snapshot_aware_validation.py`

## License

Code in this repository is released under the MIT License (see `LICENSE`). The underlying
NCBI taxonomy and RefSeq data are public domain / US government works and are not
redistributed here.

## Citation

If you use this code, please cite the manuscript above. A permanent DOI for this
repository is provided by Zenodo: https://doi.org/10.5281/zenodo.21859349 (concept DOI,
always resolves to the latest archived version).
