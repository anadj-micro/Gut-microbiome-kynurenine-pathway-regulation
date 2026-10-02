# Revised Figure S2 RNA-seq

Run `rnaseq.m` from any MATLAB folder. It calls `run_figureS2_rnaseq.m`. Requires MATLAB R2024b, Statistics and Machine Learning Toolbox, and Bioinformatics Toolbox (`nbintest`). `legacy_rnaseq.m` is the superseded original analysis, not the workflow for the revised manuscript.

Inputs are in `current_go_inputs/`: original integer `raw_counts_C_AVN.csv` and `raw_counts_C_LPS.csv` from Zenodo record 23061607's RNAseq.zip; frozen `mgi.gaf.gz`, `go-basic.obo`, and `MRK_ENSEMBL.rpt`. The annotation snapshots were acquired July 16, 2026: GAF generated May 21 (GO version April 27), GO basic release June 15. Original sources are https://current.geneontology.org/annotations/mgi.gaf.gz, https://current.geneontology.org/ontology/go-basic.obo, and https://www.informatics.jax.org/downloads/reports/MRK_ENSEMBL.rpt. Do not replace them with a later current download if exact reproduction is intended.

The workflow aligns original count matrices, verifies shared control counts, runs negative-binomial differential-expression tests with BH correction, and uses median-ratio normalized counts for display. Direct, active Biological Process annotations exclude NOT qualifiers. Enrichment uses the upper hypergeometric tail and BH across active annotated terms separately within organ, contrast, and direction. No GO ancestor propagation is added.

GO:0006569 replaces obsolete GO:0019441 and maps to ten genes. Intestinal AVN reduces Ido1 (q = 2.53369e-5); LPS reduces Tdo2 (q = 0.0118143). Neither downregulated gene set shows significant enrichment of the ten-gene GO set (AVN P = 0.161244; LPS P = 0.0959023; q = 1 for both). These are preliminary n=3/group experiments.

Tables go to `results/`; RNA-seq panels supporting Figure S2C–D go to `panels/` as EPS, PDF, PNG. For continuity, some result CSV filenames retain their original `fig1_rnaseq_` prefix; these are tables, not Figure 1 panels. The script never modifies manuscript figure directories.
