# Remove likely contaminants from 18S amplicon data (contaminants in the
# extraction and PCR reagents, not cross-contamination between wells,
# which is trickier to detect)
# This step follows initial processing in DADA2 and generation of ASV table
# found in script 'V9_DADA2_generate_ASV_table.R'

# Will use the decontam package (https://github.com/benjjneb/decontam)

# The 18S V9 primers tend to amplify some bacterial lineages, so the first
# step will be to use the Silva classification to remove non-target domains.

# Load required packages ####
library(tidyverse); packageVersion("tidyverse")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
set.seed(65536)

#settings
theme_set(theme_bw())

# Load data ####
# Combine ASV count table and taxonomy table from V9_DADA2_generate_ASV_table.R
# script, plus metadata table, into a phyloseq object
# Load metadata, with first column (Sample) as row names to match ASV table
metadata <- read.table("data/raw/malawi_tilapia_metadata.txt", header = TRUE,
                       row.names = 1)
seqtab <- readRDS("data/process/V9_seqtab_nochim.rds")
tax <- readRDS("data/process/V9_tax_final_silva.rds") # Silva taxonomy


# Combine into phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(tax))

# Remove large temporary objects
rm(seqtab); rm(tax)

# To get info about phyloseq object:
ps # 2721 taxa and 126 samples
nsamples(ps)
sample_variables(ps)
rank_names(ps)
otu_table(ps)[1:5, 1:5]
tax_table(ps)[1:5, 1:5]
taxa_names(ps)[1:10]

# Library sizes ####
# First, look at library sizes by whether sample is a true sample or a neg
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

ggplot(data = df, aes(x = Index, y = LibrarySize, color = Sample_or_neg)) +
  geom_point()


# Remove bacteria and archaea, but keep NA at domain level (since these might
# have been classified with PR2)
# subset_taxa defaults to removing NAs too, so need to specify that they be kept
# ps.euk <- subset_taxa(ps, (Kingdom == "Eukaryota") | is.na(Kingdom))
# ps.euk # 2435 taxa and 126 samples

# On second thought... Keep only those sequences classified as Eukaryota.
# When I initially kept the NAs, I had a lot of sequences that couldn't
# be classified by PR2 beyond the Domain level.
ps.euk <- subset_taxa(ps, (Kingdom == "Eukaryota"))
ps.euk # 2353 taxa and 126 samples


# Extract this filtered euk-only ASV table from phyloseq object, then make a
# new phyloseq object with the PR2 classification instead of Silva
ASV.euk <- data.frame(as(otu_table(ps.euk), 'matrix'))
tax.pr2 <- readRDS("data/process/V9_tax_final_pr2.rds") # PR2 taxonomy

ps <- phyloseq(otu_table(ASV.euk, taxa_are_rows=FALSE), 
                   sample_data(metadata), 
                   tax_table(tax.pr2))

ps # 2353 taxa and 126 samples



# Separate data for decontam ####
# Need to run decontam separately on fundamentally different sample types:
# 1. water (pond + incubation), 2. swabs

# Make two separate phyloseq objects for these, each with the negative
# controls, but include mock community only with water samples (since the 
# DNA concentration is within that range whereas the other group has much
# lower DNA concentrations).

# Water + mock community + negatives (i.e. everything except swabs)
water <- subset_samples(ps, sample_type != "Swab")
water # 2353 taxa and 74 samples

# Swabs + negatives only
swab <- subset_samples(ps, sample_type == "Swab" | Sample_or_neg == "neg")
swab # 2353 taxa and 58 samples



## WATER ONLY ####
# Remove singletons (there shouldn't be any from DADA2 but will likely be some
# after subsetting)
water.pruned <- prune_taxa(taxa_sums(water) > 1, water)
water.pruned <- prune_samples(sample_sums(water.pruned) > 0, water.pruned)
water.pruned # 2324 taxa and 73 samples


# Extract ASV table, sample metadata and taxonomy table as
# base R data.frames
# get sample sums: each sample's sequencing depth
sample_data(water.pruned)$sample_sums <- sample_sums(water.pruned)
MAP <- data.frame(sample_data(water.pruned))
OTU <- data.frame(as(otu_table(water.pruned), 'matrix'))
TAX <- data.frame(tax_table(water.pruned))

# Apply prevalence filter to filter out taxa present in fewer than 2 samples
# with function(x) sum(x > 0) > 1
# On second thought, don't do this just yet... It could remove taxa that are
# very abundant but in a single sample. Run it with sum(x > 0) > 0, which
# essentially removes taxa not present in any sample (there should be any, but
# best to check anyway...)
water.pruned2 <- filter_taxa(water.pruned, function(x) sum(x > 0) > 0, TRUE)
OTU2 <- OTU[,taxa_names(water.pruned2)]
TAX2 <- TAX[taxa_names(water.pruned2),]
MAP2 <- MAP[sample_names(water.pruned2),]
water.pruned2 # 2324 taxa and 73 samples

# Test three decontam methods on water data ####
# Classify contaminants using the frequency method, the prevalence method, 
# and the combined method, to assess which is most appropriate for these data
conc2 <- MAP2$Conc_18S
neg2 <- MAP2$Sample_or_neg == "neg"

# Frequency-based:
ocf <- isContaminant(as.matrix(OTU2), conc=conc2, threshold=0.1, detailed=TRUE, normalize=TRUE, method='frequency')

# Prevalence-based:
ocp <- isContaminant(as.matrix(OTU2), neg=neg2, threshold=0.5, detailed=TRUE, normalize=TRUE, method='prevalence')

# Combined method:
occ <- isContaminant(as.matrix(OTU2), conc=conc2, neg=neg2, threshold=0.1, detailed=TRUE, normalize=TRUE, method='combined')

# Combine decontam scores from each method into a data.frame.
probcols <- data.frame(row.names=rownames(ocf),prob.f=ocf$p.freq, prob.p=ocp$p.prev, prob.c=occ$p)

# Crude comparison of frequency and prevalence contaminant assignment
table(probcols$prob.f<0.1, probcols$prob.p<0.5)


# Plot score distribution

TAXann <- rbind(cbind(TAX2, Score=probcols$prob.f, Method="Frequency"),
                cbind(TAX2, Score=probcols$prob.p, Method="Prevalence"),
                cbind(TAX2, Score=probcols$prob.c, Method="Combined"))

TAXann$Method <- factor(TAXann$Method, levels=c('Frequency','Prevalence','Combined'))

# Histogram of decontam scores for each method
histo <- ggplot(TAXann, aes(x=Score))
histo <- histo + geom_histogram() + labs(x = 'decontam Score', y ='Number ASVs') + 
  facet_wrap(~ Method, nrow=1) +
  theme(legend.position = "bottom")
histo

# Final water decontam choice ####
# There's no clear bimodal pattern with the frequency/combined approaches,
# making it difficult to select a threshold. Prevalence might be best here, with
# threshold of 0.5, but need to examine the plots of individual ASVs 
# Prevalence method:
water.cont <- isContaminant(as.matrix(OTU2), conc=conc2, neg=neg2, threshold=0.5,
                            detailed=TRUE, normalize=TRUE, method='prevalence')



## SWABS ONLY ####
# Remove singletons (there shouldn't be any from DADA2 but likely some after
# subsetting)
swab.pruned <- prune_taxa(taxa_sums(swab) > 1, swab)
swab.pruned <- prune_samples(sample_sums(swab.pruned) > 0, swab.pruned)
swab # 2353 taxa and 58 samples
swab.pruned # 303 taxa and 57 samples

# Extract ASV table, sample metadata and taxonomy table as
# base R data.frames
# get sample sums: each sample's sequencing depth
sample_data(swab.pruned)$sample_sums <- sample_sums(swab.pruned)
MAP <- data.frame(sample_data(swab.pruned))
OTU <- data.frame(as(otu_table(swab.pruned), 'matrix'))
TAX <- data.frame(tax_table(swab.pruned))

# Apply prevalence filter to filter out taxa present in fewer than 1 sample. 

swab.pruned2 <- filter_taxa(swab.pruned, function(x) sum(x > 0) > 0, TRUE)
OTU2 <- OTU[,taxa_names(swab.pruned2)]
TAX2 <- TAX[taxa_names(swab.pruned2),]
MAP2 <- MAP[sample_names(swab.pruned2),]

swab.pruned2 # 303 taxa and 57 samples

# Test three decontam methods on swabs ####
# Classify contaminants using both methods plus combined
conc2 <- MAP2$Conc_18S
neg2 <- MAP2$Sample_or_neg == "neg"

# Frequency-based:
ocf <- isContaminant(as.matrix(OTU2), conc=conc2, threshold=0.1, detailed=TRUE, normalize=TRUE, method='frequency')

# Prevalence-based:
ocp <- isContaminant(as.matrix(OTU2), neg=neg2, threshold=0.5, detailed=TRUE, normalize=TRUE, method='prevalence')

# Combined method:
occ <- isContaminant(as.matrix(OTU2), conc=conc2, neg=neg2, threshold=0.1, detailed=TRUE, normalize=TRUE, method='combined')

# Combine scores from each method into a data.frame.
probcols <- data.frame(row.names=rownames(ocf),prob.f=ocf$p.freq, prob.p=ocp$p.prev, prob.c=occ$p)

# Crude comparison of prevalence and combined contaminant assignment
table(probcols$prob.f<0.1, probcols$prob.p<0.5)


# Plot score distribution

TAXann <- rbind(cbind(TAX2, Score=probcols$prob.f, Method="Frequency"),
                cbind(TAX2, Score=probcols$prob.p, Method="Prevalence"),
                cbind(TAX2, Score=probcols$prob.c, Method="Combined"))

TAXann$Method <- factor(TAXann$Method, levels=c('Frequency','Prevalence','Combined'))

# Histogram of decontam scores for each method
histo <- ggplot(TAXann, aes(x=Score))
histo <- histo + geom_histogram() + labs(x = 'decontam Score', y ='Number ASVs') + 
  facet_wrap(~ Method, nrow=1) +
  theme(legend.position = "bottom")
histo

# Final swab decontam choice ####
# Prevalence method is most appropriate given the low DNA conc of this group
swab.cont <- isContaminant(as.matrix(OTU2), conc=conc2, neg=neg2, threshold=0.5,
                           detailed=TRUE, normalize=TRUE, method='prevalence')



# Examine ASVs identified as contaminants ####
# Give the total number of ASVs identified as contaminants
table(water.cont$contaminant) # 12
table(swab.cont$contaminant) # 15

# Give the rank of the first few most abundant contaminants
head(which(water.cont$contaminant)) # 4, 35, 176, 513, 722, 746
head(which(swab.cont$contaminant)) # 4, 16, 21, 64, 89, 97



# Frequency of contaminant ASVs against DNA concentration
plot_frequency(water.pruned2, taxa_names(water.pruned2)[which(water.cont$contaminant)[1:12]], conc="Conc_16S") + 
  xlab("DNA Concentration (ng/uL)") + ggtitle("Contaminant ASVs in water samples")


plot_frequency(swab.pruned2, taxa_names(swab.pruned2)[which(swab.cont$contaminant)[1:15]], conc="Conc_16S") + 
  xlab("DNA Concentration (ng/uL)") + ggtitle("Contaminant ASVs in swab samples")



# Prevalence in true samples vs neg controls
# with colours showing identified contaminants
# Water
# Make phyloseq object of presence-absence in negative controls and true samples
wt.pa <- transform_sample_counts(water.pruned2, function(abund) 1*(abund>0))
wt.pa.neg <- prune_samples(sample_data(wt.pa)$Sample_or_neg == "neg", wt.pa)
wt.pa.pos <- prune_samples(sample_data(wt.pa)$Sample_or_neg == "sample", wt.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(wt.pa.pos), pa.neg=taxa_sums(wt.pa.neg),
                    contaminant=water.cont$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point(position = "jitter") +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  ggtitle("Water samples")

# swab
# Make phyloseq object of presence-absence in negative controls and true samples
an.pa <- transform_sample_counts(swab.pruned2, function(abund) 1*(abund>0))
an.pa.neg <- prune_samples(sample_data(an.pa)$Sample_or_neg == "neg", an.pa)
an.pa.pos <- prune_samples(sample_data(an.pa)$Sample_or_neg == "sample", an.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(an.pa.pos), pa.neg=taxa_sums(an.pa.neg),
                    contaminant=swab.cont$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point(position = "jitter") +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  ggtitle("swab samples")



## Remove likely contaminants from phyloseq objects ####
# Before
water.pruned2 # 2324 taxa and 73 samples
swab.pruned2 # 303 taxa and 57 samples

# Remove contaminants from phyloseq objects
water.nocontam <- prune_taxa(!water.cont$contaminant, water.pruned2)
water.nocontam # 2312 taxa and 73 samples

swab.nocontam <- prune_taxa(!swab.cont$contaminant, swab.pruned2)
swab.nocontam # 288 taxa and 57 samples



## Merge the two phyloseq objects back together
# First, need to remove the negatives from swab, as these were included in 
# both groups.

swab.nocontam2 <- subset_samples(swab.nocontam, Sample_or_neg != "neg")


# Merge the two phyloseq objects
ps.nocontam <- merge_phyloseq(water.nocontam, swab.nocontam2)

water.nocontam
swab.nocontam2

ps.nocontam # 2343 taxa and 125 samples

# Write to disk
saveRDS(ps.nocontam, "data/process/V9_phyloseq_obj_after_decontam.rds")


## REMOVE NON-TARGET SEQUENCES ####
# Remaining QC tasks before we can proceed to proper exploration/analysis:
# 1. Removing samples that are not needed in downstream analyses
# 2. Removing non-target sequences (mitochondria, chloroplasts)
# 3. Select Bacteria and Archaea (not Eukaryotes or sequences unclassified
# at the Domain level)
# 4. Pulling out mock community positive controls and formatting them for
# seq.error in mothur to determine overall sequencing error rate


# Load phyloseq object created at the end of the decontam steps
ps <- readRDS("data/process/V9_phyloseq_obj_after_decontam.rds")

# Remove control samples (mock community positives and extraction/PCR negatives)
# that are not needed in downstream analyses
# Sequences in negs have already been addressed in the decontam script
ps.no.controls <- subset_samples(ps, sample_type != "control")
ps.no.controls # 2343 taxa and 116 samples

# Remove ASVs that are unclassified at the Domain level or that don't classify as euk
ps.domain <- subset_taxa(ps.no.controls, Kingdom == "Eukaryota")
ps.domain # 2343 taxa and 116 samples

# There aren't any organelles to remove here, but I'll use this object name
# to preserve downstream code
ps.no.organelles <- ps.domain
ps.no.organelles # 2343 taxa and 116 samples

# Save phyloseq object before filtering out small libraries and caudal/gill swabs
saveRDS(ps.no.organelles, "data/process/V9_phyloseq_obj_final.rds")



# A lot of swab samples have a small number of reads, but I'll keep them for now.
# Will just filter out the smallest ones by choosing a threshold of 200.
ps1 <- prune_samples(sample_sums(ps.no.organelles) >= 200, ps.no.organelles)

# Remove caudal and gill swabs, and any ASVs with fewer than 3 counts in total
ps2 <- ps1 %>% subset_samples(!(sample_type_detailed %in% c("Swab_cd", "Swab_gl")))
ps.no.organelles2 <- prune_taxa(taxa_sums(ps2) > 2, ps2)


# Check again for ASVs that don't occur in any sample (happens when samples are
# removed, especially mock community positive controls, which have different taxa)
# ps.no.organelles2 <- filter_taxa(ps.no.organelles, function(x) sum(x > 0) > 0, TRUE)

# Get sample sums (i.e. sequencing depth) and append
sample_data(ps.no.organelles2)$sample_sums <- sample_sums(ps.no.organelles2)
ps.no.organelles2 # 2293 taxa and 91 samples


# Change ASV names ####
# With the combined phyloseq object, change the ASV names from the full DNA
# sequence to much shorter and more manageable names ASVx, where x is the
# ASV number in decreasing order of overall abundance. The ASV DNA sequences 
# are stored in a separate object called `dna`, kept with the `phyloseq` object.
seqs <- Biostrings::DNAStringSet(taxa_names(ps.no.organelles2))
names(seqs) <- taxa_names(ps.no.organelles2)
ps.no.organelles2 <- merge_phyloseq(ps.no.organelles2, seqs)
taxa_names(ps.no.organelles2) <- paste0("ASV", seq(ntaxa(ps.no.organelles2)))

# Save this object prior to phylogenetic tree step in script V9_generate_tree.R
saveRDS(ps.no.organelles2, "data/process/V9_phyloseq_obj_after_organelle_removal.rds")
