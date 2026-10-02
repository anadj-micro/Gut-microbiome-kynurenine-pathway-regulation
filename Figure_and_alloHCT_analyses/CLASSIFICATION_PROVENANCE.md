# SILVA 138.1 classification-generation audit

The existing `inputs/asv_taxonomy_reclassified.csv` remains the authoritative frozen Figure 1 input. `reclassify_human_asvs.R` uses only original public sequences/counts/sample selection to classify the 2,106 ASVs observed in the 92-sample cohort. It preserves IDs and sequences and writes new assignments/bootstrap confidence to a separate `classification_audit/` folder; it does not change abundance counts or the frozen table.

Requires R 4.6.1, DADA2 1.40.0, data.table, digest, and the materialized SILVA 138.1 training FASTA already in `16S_amplicon_sequence_processing/dada2_databases/`. The original mouse DADA2 package version is undocumented; this matches reference and settings, not the exact historical mouse software environment. Settings: seed 100, minBoot 80, tryRC FALSE, one thread, kingdom through genus, no species assignment. Versions and reference SHA-256 are exported.

October 2 rerun: zero changes at kingdom, phylum, class, order, or genus; five family assignment changes relative to the first archived run. Bootstrap scores also vary. This is consistent with the previously documented DADA2 C++ tie-breaking behavior: an R seed alone does not ensure bitwise reproduction. Do not silently overwrite the frozen manuscript input or select a rerun by association strength. Exact downstream figure/statistic reproduction uses the frozen assignments; the new script documents generation and audits its variability.

No new family associations or manuscript figures were adopted from this rerun.
