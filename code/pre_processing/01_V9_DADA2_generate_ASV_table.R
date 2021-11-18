# Malawi tilapia microbiomes (Newton Aquaculture Project)

# DADA2 short read 18S rRNA V9 amplicon pipeline for generating ASV table 
# from raw data files (*.fq.gz)


# Much of this script is based on DADA2 tutorials:
# https://www.bioconductor.org/packages/release/bioc/vignettes/dada2/inst/doc/dada2-intro.html
# https://benjjneb.github.io/dada2/tutorial.html
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex

library(dada2); packageVersion("dada2")
set.seed(65536)

# Load fastq files in .gz format.
# First, set the path to the data files (using relative paths from project root)
path <- "data/raw/18S"
list.files(path)

fnFs <- sort(list.files(path, pattern="_r1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_r2.fq.gz", full.names = TRUE))

# Extract sample names assuming filenames have format: 2725_SAMPLENAME_rX.fq.gz
# Most sample names include the prefix MLW, which should be removed. Some controls
# from that sequencing run don't have the MLW prefix.
# First, remove all characters after (and including) "_r"
sample.names.temp <- sapply(strsplit(basename(fnFs), "_r"), `[`, 1)

# Then remove the sequencing run number "2725_" and anything preceding it
sample.names <- sapply(strsplit(sample.names.temp, "2725_"), `[`, 2)

# Finally, remove the prefix "MLW"
sample.names <- gsub("MLW_", "", sample.names)

write(sample.names, file="data/process/V9_names.txt")


# Examine the quality profile of the data
plotQualityProfile(fnFs[1:4]) # Forward

plotQualityProfile(fnRs[1:4]) # Reverse


# Trim the forward and reverse reads at position 100
# Filter out all reads with more than  maxN=0 ambiguous nucleotides
# Filter out all reads with more than two expected errors
# Store filtered output files as gzipped fastq files (compress=TRUE)
# Keep only sequences when both F and R reads are kept

# Place filtered files in data/process/V9_filtered directory
filtFs <- file.path("data/process/V9_filtered", paste0(sample.names, 
                                                      "_F_filt.fastq.gz"))
filtRs <- file.path("data/process/V9_filtered", paste0(sample.names, 
                                                      "_R_filt.fastq.gz"))

out <- filterAndTrim(fwd=fnFs, filt=filtFs, rev=fnRs, filt.rev=filtRs,
                     trimLeft=0, truncLen=c(100, 100), maxN=0, maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE,
                     verbose=TRUE)

head(out)


## Learn error rate
# uses a parametric model of the errors introduced by PCR amplification and 
# sequencing

errF <- learnErrors(filtFs, nbases = 1e8, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, nbases = 1e8, multithread=TRUE, randomize=TRUE)


# Plot errors to check that estimated error (black line) is a good fit for 
# observed error (grey dots)
plotErrors(errF[1:9], nominalQ=TRUE)
plotErrors(errR[1:9], nominalQ=TRUE)


## Dereplicate
# derepFastq maintains a summary of the quality information for each 
# dereplicated sequence in $quals
derepF <- derepFastq(filtFs, verbose=TRUE)
derepR <- derepFastq(filtRs, verbose=TRUE)

# The sample names in these objects are initially the file names of the samples, 
# This sets them to the sample names for the rest of the workflow
names(derepF) <- sample.names
names(derepR) <- sample.names


## Infer sample composition
# Using pool=TRUE means sequences from all samples will be pooled prior to ASV
# detection. This is better at detecting rare ASVs but more computationally 
# intensive than the default pool=FALSE
# Can use pool="pseudo" as an intermediate
dadaF <- dada(derepF, err=errF, multithread=TRUE, pool=TRUE)
dadaR <- dada(derepR, err=errR, multithread=TRUE, pool=TRUE)

# Merge pairs
# trimOverhang=TRUE is probably unnecessary but will include just in case
# minOverlap default is only 12, so will set it higher (reads should overlap 
# substantially here)
merger <- mergePairs(dadaF, derepF, dadaR, derepR, trimOverhang=TRUE, 
                      minOverlap=25, verbose=TRUE)

# Remove large temporary files
rm(derepF); rm(derepR)

# Make and save the ASV table
seqtab <- makeSequenceTable(merger)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# This amplicon varies in length, but very short and very long are likely nonspecifics.
# This sets the length range to 90-150
seqtab.trim <- seqtab[,nchar(colnames(seqtab)) %in% seq(90,150)]
table(nchar(getSequences(seqtab.trim)))

saveRDS(seqtab.trim, "data/process/V9_DADA2_seqtab.rds")
# seqtab.trim <- readRDS("data/process/V9_DADA2_seqtab.rds")


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab.trim, method="consensus", multithread=TRUE,
                             verbose=TRUE)

# Save seqtab.nochim to allow re-running of taxonomy without needing to repeat everything!
saveRDS(seqtab.nochim, "data/process/V9_seqtab_nochim.rds")
# seqtab.nochim <- readRDS("data/process/V9_seqtab_nochim.rds")


# Make summary table of no. sequences before, during, after processing
getN <- function(x) sum(getUniques(x))

# 'out' has the input and filtered columns
denoisedF <- sapply(dadaF, getN)
denoisedR <- sapply(dadaR, getN)
merged <- sapply(merger, getN)
no.chim <- rowSums(seqtab.nochim)

track <- cbind(out, denoisedF, denoisedR, merged, no.chim)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
rownames(track) <- sample.names

head(track)
tail(track)

write.csv(track, "data/process/V9_DADA2_summary.csv")



# Assign taxonomy: use both Silva (for identifying bact & arch seqs) and PR2

# Silva taxonomy using naive Bayesian classifier, with a minimum bootstrap
# confidence of 60% to identify bact & arch sequences
taxa.silva <- assignTaxonomy(seqtab.nochim, "data/references/silva_nr_v138_train_set.fa.gz",
                             minBoot = 60, multithread = TRUE)

# Write to disk
saveRDS(taxa.silva, "data/process/V9_tax_final_silva.rds")


# PR2 4.12 database, formatted for DADA2
# https://github.com/pr2database/pr2database/releases/download/v4.12.0/pr2_version_4.12.0_18S_dada2.fasta.gz
# PR2 has different taxlevels so need to specify with taxLevels = xx
# minBoot = 60 because this amplicon is so short
taxa.pr2 <- assignTaxonomy(seqtab.nochim, "data/references/pr2_version_4.12.0_18S_dada2.fasta.gz", 
                           taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),
                           minBoot = 60, multithread = TRUE)

# Write to disk
saveRDS(taxa.pr2, "data/process/V9_tax_final_pr2.rds")


# Next: remove likely contaminants with the decontam package (separate script)