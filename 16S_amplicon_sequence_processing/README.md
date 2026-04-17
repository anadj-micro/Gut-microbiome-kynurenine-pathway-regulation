Scripts and data necessary for 16S rDNA amplicon sequence analysis with dada2.

The structure of the folder needs to be maintained. Fastq files go to "output" folder, database against which ASV alignment will be performed goes to "dada2_databases" folder, and final taxonomy and counts table are generated in "final_tables" folder.

Order:
1. Go to NCBI and download fastq files (PRJNA1435869). Put the fastq file in "output" subfolder.
2. Run dadascript_noprior.R. This takes fastq files, filters them, merges them, removes chimeras, and generates .rds files.
3. Run convertseqtab2newtbls_dada2_taxonomy.R. This takes .rds files, aligns sequences against database (in this case Silva v.138.1), and generates counts and taxonomy tables.
4. Run scParseTaxonomicLevels.m. This separates counts table to generate one abundance table per taxonomic level.