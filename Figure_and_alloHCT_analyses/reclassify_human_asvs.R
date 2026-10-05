#!/usr/bin/env Rscript
# Reclassify observed original sequences without changing counts or frozen inputs.
# Requires DADA2 1.40.0, data.table, digest. Configure R_LIBS_USER if necessary.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)))
here <- dirname(script)
library(dada2)
library(data.table)
library(digest)
stopifnot(as.character(packageVersion("dada2")) == "1.40.0")
input <- file.path(here, "inputs")
output <- file.path(here, "classification_audit")
dir.create(output, showWarnings=FALSE)
metabolites <- fread(file.path(input, "tblMetabolitesUnstackedStoolPlasma.txt"))
embedding <- fread(file.path(input, "current_taxumap_embedding.csv"))
cohort <- intersect(metabolites$SampleID, embedding[[1]])
stopifnot(length(cohort) == 92)
counts <- fread(file.path(input, "tblcounts_asv_melt.csv"))
observed <- sort(unique(counts[SampleID %in% cohort & Count>0, ASV]))
original <- fread(file.path(input, "tblASVtaxonomy_silva132_v4v5_filter.csv"))
stopifnot(!anyDuplicated(original$ASV))
old <- original[match(observed, ASV)]
stopifnot(nrow(old) == 2106, !anyNA(old$ASV), !anyDuplicated(old$Sequence))
reference <- file.path(dirname(here), "16S_amplicon_sequence_processing",
                       "dada2_databases/silva_nr99_v138.1_train_set.fa.gz")
stopifnot(file.exists(reference))
set.seed(100)
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
classified <- assignTaxonomy(old$Sequence, reference, minBoot=80, tryRC=FALSE,
                             outputBootstraps=TRUE, taxLevels=ranks,
                             multithread=FALSE, verbose=TRUE)
stopifnot(identical(rownames(classified$tax), old$Sequence))
revised <- copy(old)
columns <- setdiff(names(old), c("ASV", "Sequence"))
setnames(revised, columns, paste0("Old", columns))
for (rank in ranks) {
    revised[[paste0("New", rank)]] <- classified$tax[,rank]
    revised[[paste0("Bootstrap", rank)]] <- classified$boot[,rank]
}
fwrite(revised, file.path(output, "regenerated_taxonomy.csv"), na="NA")
# Keep the first archived classification as the manuscript input. DADA2's
# C++ tie-breaking can vary even with R's seed; never select a favorable rerun.
frozen <- fread(file.path(input, "asv_taxonomy_reclassified.csv"), na.strings="NA")
frozen <- frozen[match(revised$ASV, ASV)]
stopifnot(identical(frozen$Sequence, revised$Sequence))
differences <- rbindlist(lapply(ranks, function(rank) {
    a <- revised[[paste0("New",rank)]]
    b <- frozen[[paste0("New",rank)]]
    changed <- ifelse(is.na(a) | is.na(b), xor(is.na(a),is.na(b)), a!=b)
    data.table(Rank=rank, ChangedAssignments=sum(changed),
               ChangedBootstraps=sum(revised[[paste0("Bootstrap",rank)]] !=
                                     frozen[[paste0("Bootstrap",rank)]]))
}))
fwrite(differences, file.path(output, "comparison_with_frozen.csv"))
fwrite(data.table(Reference=basename(reference),
                  SHA256=digest(reference, algo="sha256", file=TRUE),
                  DADA2=as.character(packageVersion("dada2")), R=R.version.string,
                  Seed=100, MinBoot=80, TryRC=FALSE),
       file.path(output, "classifier_provenance.csv"))
writeLines(capture.output(sessionInfo()), file.path(output,"session_info.txt"))
print(differences)
