#!/usr/bin/env Rscript
# Script for inferring sequence variants with dada2. Follow with convertseqtab2newtbls_dada2_taxonomy.R. 
# Required: dada2
# Line 23: This command is written for fastq files generated in this project. If re-using
# this script for other projects, update the line according to the names of the new fastq files.
# Setup path to the "output" folder where fastq files are. For example:

# path <- "~/Documents/Analysis/dada2/output"

library(dada2); packageVersion("dada2")

#forward
fnFs <- Sys.glob(file.path(path, "*_R1*.fastq.gz"))

#backward
fnRs <- Sys.glob(file.path(path, "*_R2*.fastq.gz"))

forwardfileinfo <- file.info(fnFs)
forwardsizes <- forwardfileinfo[,"size"]
reversefileinfo <- file.info(fnRs)
reversesizes <- forwardfileinfo[,"size"]

#get sample names (update depending on filename conventions)
sample.names <- sapply(strsplit(sapply(strsplit(basename(fnFs), ".fastq.gz"), `[`, 1), "fastq_"), `[`,2)
#
#temporary files for storing filtered data
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(180,180),
                     maxN=0, maxEE=2, truncQ=2,
                     compress=TRUE, multithread=TRUE)

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

if(!inherits(derepFs,"list")){
  derepFs <- list(derepFs)
  derepRs <- list(derepRs)
}
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

dadaFs <- dada(derepFs, err=errF,multithread=TRUE)

dadaRs <- dada(derepRs, err=errR,multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

seqtabfinal <- seqtab.nochim[,nchar(colnames(seqtab.nochim)) %in% seq(300,360)]

outfile <- file.path(path, 'seqtab_noprior.rds')
saveRDS(seqtabfinal, file=outfile)

outfile <- file.path(path, 'errR.rds')
saveRDS(errR, file=outfile)

outfile <- file.path(path, 'errF.rds')
saveRDS(errF, file=outfile)

outfile <- file.path(path, 'derepFs.rds')
saveRDS(derepFs, file=outfile)

outfile <- file.path(path, 'derepRs.rds')
saveRDS(derepRs, file=outfile)
