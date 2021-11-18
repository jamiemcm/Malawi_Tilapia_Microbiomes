# Malawi tilapia microbiomes (Newton Aquaculture Project)

# DADA2 short read 16S rRNA V4 amplicon pipeline for generating ASV table 
# from raw data files (*.fq.gz)


# Much of this script is based on DADA2 tutorials:
# https://www.bioconductor.org/packages/release/bioc/vignettes/dada2/inst/doc/dada2-intro.html
# https://benjjneb.github.io/dada2/tutorial.html
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex

library(dada2); packageVersion("dada2")
library(DECIPHER); packageVersion("DECIPHER")
set.seed(65536)

# Load fastq files in .gz format.
# First, set the path to the data files (using relative paths from project root)
path <- "data/raw/"
list.files(path)

fnFs <- sort(list.files(path, pattern="_r1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_r2.fq.gz", full.names = TRUE))

# Extract sample names assuming filenames have format: 2716_SAMPLENAME_rX.fq.gz
# Most sample names include the prefix MLW, which should be removed. Some controls
# from that sequencing run don't have the MLW prefix.
# First, remove all characters after (and including) "_r"
sample.names.temp <- sapply(strsplit(basename(fnFs), "_r"), `[`, 1)

# Then remove the sequencing run number "2716_" and anything preceding it
sample.names <- sapply(strsplit(sample.names.temp, "2716_"), `[`, 2)

# Finally, remove the prefix "MLW"
sample.names2 <- gsub("MLW_", "", sample.names)

write(sample.names2, file="data/process/V4_names.txt")


# Examine the quality profile of the data
plotQualityProfile(fnFs[1:4]) # Forward

plotQualityProfile(fnRs[1:4]) # Reverse


# Trim the forward reads at position 200 and reverse at position 160
# Filter out all reads with more than  maxN=0 ambiguous nucleotides
# Filter out all reads with more than two expected errors
# Store filtered output files as gzipped fastq files (compress=TRUE)
# Keep only sequences when both F and R reads are kept

# Place filtered files in data/process/V4_filtered directory
filtFs <- file.path("data/process/V4_filtered", paste0(sample.names2, 
                                                      "_F_filt.fastq.gz"))
filtRs <- file.path("data/process/V4_filtered", paste0(sample.names2, 
                                                      "_R_filt.fastq.gz"))

out <- filterAndTrim(fwd=fnFs, filt=filtFs, rev=fnRs, filt.rev=filtRs,
                     truncLen=c(200, 160), maxN=0, maxEE=c(2,2), truncQ=2, 
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose=TRUE)

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
derepF2 <- derepFastq(filtFs, verbose=TRUE)
derepR2 <- derepFastq(filtRs, verbose=TRUE)

# The sample names in these objects are initially the file names of the samples, 
# This sets them to the sample names for the rest of the workflow
names(derepF2) <- sample.names2
names(derepR2) <- sample.names2


## Infer sample composition
# Using pool=TRUE means sequences from all samples will be pooled prior to ASV
# detection. This is better at detecting rare ASVs but more computationally 
# intensive than the default pool=FALSE
# Can use pool="pseudo" as an intermediate
dadaF2 <- dada(derepF2, err=errF, multithread=TRUE, pool=TRUE)
dadaR2 <- dada(derepR2, err=errR, multithread=TRUE, pool=TRUE)

# Merge pairs
# trimOverhang=TRUE is probably unnecessary but will include just in case
# minOverlap default is only 12, so will set it higher (reads should overlap 
# substantially here)
merger2 <- mergePairs(dadaF2, derepF2, dadaR2, derepR2, trimOverhang=TRUE, 
                      minOverlap=100, verbose=TRUE)

# Remove large temporary files
rm(derepF2); rm(derepR2)

# Make and save the ASV table
seqtab2 <- makeSequenceTable(merger2)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab2)))

# This amplicon should be 253 bp. Very different lengths are likely 
# nonspecifics so should be removed. This sets the length range to 250-256.
seqtab.trim2 <- seqtab2[,nchar(colnames(seqtab2)) %in% seq(250,256)]
table(nchar(getSequences(seqtab.trim2)))

saveRDS(seqtab.trim2, "data/process/V4_DADA2_seqtab.rds")
# seqtab.trim2 <- readRDS("data/process/V4_DADA2_seqtab.rds")


# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab.trim2, method="consensus", multithread=TRUE,
                             verbose=TRUE)

# Save seqtab.nochim to allow re-running of taxonomy without needing to repeat everything!
saveRDS(seqtab.nochim, "data/process/V4_seqtab_nochim.rds")
# seqtab.nochim <- readRDS("data/process/V4_seqtab_nochim.rds")


# Make summary table of no. sequences before, during, after processing
getN <- function(x) sum(getUniques(x))

# 'out' has the input and filtered columns
denoisedF <- sapply(dadaF2, getN)
denoisedR <- sapply(dadaR2, getN)
merged <- sapply(merger2, getN)
no.chim <- rowSums(seqtab.nochim)

track <- cbind(out, denoisedF, denoisedR, merged, no.chim)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochim")
rownames(track) <- sample.names2

head(track)
tail(track)

write.csv(track, "data/process/V4_DADA2_summary.csv")



# Assign taxonomy using naive Bayesian classifier, with a minimum bootstrap
# confidence of 80%
taxa <- assignTaxonomy(seqtab.nochim, "data/references/silva_nr_v138_train_set.fa.gz",
                       minBoot = 80, multithread = TRUE)

# Save taxa as RDS object in case the next step crashes...
saveRDS(taxa, "data/process/V4_taxa_temp.rds")
# taxa <- readRDS("data/process/V4_taxa_temp.rds")

# Add species assignment to ASVs with 100% match to references:
taxa <- addSpecies(taxa, "data/references/silva_species_assignment_v138.fa.gz")

# How many have a species assignment?
length(taxa[, 7]) - sum(is.na(taxa[, 7]))


# For comparison, assign taxonomy using DECIPHER
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("data/references/SILVA_SSU_r138_2019.RData")
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa.n <- x$taxon[m]
  taxa.n[startsWith(taxa.n, "unclassified_")] <- NA
  taxa.n
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim) 


# Write to disk
saveRDS(seqtab.nochim, "data/process/V4_ASVtable_final.rds")
saveRDS(taxa, "data/process/V4_tax_bayesian_final.rds")
saveRDS(taxid, "data/process/V4_tax_DECIPHER_final.rds")


# Extract standard tables in universal formats rather than R-specific
# Give seq headers more manageable names (ASV_1, ASV_2...) instead of full seq
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "data/process/V4_ASVs.fasta")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "data/process/V4_ASV_counts.tsv", sep="\t", 
            quote=F, col.names=NA)

# Bayesian tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "data/process/V4_taxonomy_Bayesian.tsv", sep="\t", 
            quote=F, col.names=NA)

# DECIPHER tax table:
asv_taxD <- taxid
row.names(asv_taxD) <- sub(">", "", asv_headers)
write.table(asv_taxD, "data/process/V4_taxonomy_DECIPHER.tsv", sep="\t", 
            quote=F, col.names=NA)

# Next: remove likely contaminants with the decontam package (separate script)