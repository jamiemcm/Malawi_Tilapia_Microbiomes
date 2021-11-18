# Phylogenetic tree ####
# This script won't run on a laptop - it must be run on a cluster.
 
# Constructs a phylogenetic tree of final clean sequences for phylogeny-aware
# distance methods (e.g. unifrac) and for agglomerating by phylogeny instead
# of taxonomy.

# This follows the example from https://f1000research.com/articles/5-1492/v2,
# and uses the DECIPHER and phangorn packages for alignment and tree construction,
# respectively.

# It follows initial processing in DADA2 and generation of ASV table found in 
# script 'V4_DADA2_generate_ASV_table.R', and removal of likely contaminants
# with the decontam package, in script 'V4_decontam_remove_contaminants.R'

# On ISCA, I run this script within my conda R environment that already has
# the required packages

library(phyloseq)
library(DECIPHER)
library(phangorn)

# Load the phyloseq object created at the end of the decontam step
ps <- readRDS("data/process/V4_phyloseq_obj_after_organelle_removal.rds")

# # Test with a small subset of samples/ASVs
# ps1 <- subset_samples(ps, sample_type == "Pond" & Pond == "1")
# ps <- filter_taxa(ps1, function(x) sum(x > 9) > 4, TRUE)

# Phylogenetic tree ####
seqs <- DNAStringSet(refseq(ps))
alignment <- AlignSeqs(seqs, anchor=NA)

# Phangorn code - takes too long to run, even on AWS!
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

dm <- dist.ml(phang.align)

saveRDS(dm, "data/process/phang.ml.dist.rds")

treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

ps.merged <- merge_phyloseq(ps, phy_tree(fitGTR$tree))

# Save the final phyloseq object that now includes the tree
saveRDS(ps.merged, "data/process/V4_phyloseq_obj_with_tree.rds")