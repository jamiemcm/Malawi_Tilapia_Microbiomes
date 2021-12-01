# Malawi_Tilapia_Microbiomes
## Relationships between pond water and tilapia skin microbiomes in aquaculture ponds in Malawi
Jamie McMurtrie, Shayma Alathari, Dominique Chaput, David Bass, Camerson Chambi, Joseph Nagoli, Jérôme Delamare-Deboutteville, Mohan Chadag, Joanne Cable, Ben Temperton, Charles R. Tyler

### Abstract
Intensification of fish farming practices is being driven by the demand for increased food production to support a rapidly growing global human population, particularly in lower-middle income countries. Intensification of production, however, increases the risk of disease outbreaks and thus likelihood for crop losses. The microbial communities that colonise the skin mucosal surface of fish are poorly understood, but are important in maintaining fish health and resistance against disease. This skin microbial community is susceptible to disruption through stressors associated with transport, handling and the environment of intensive practices, and this risks the propagation of disease-causing pathogens. In this study, we characterised the microbial assemblages found on tilapia skin — the most widely farmed finfish globally — and in the surrounding water of seven earthen aquaculture ponds from two pond systems in distinct geographic regions in Malawi. Metabarcoding approaches were used to sequence the prokaryotic and microeukaryotic communities. We found 92% of prokaryotic amplicon sequence variants were common to both skin and water samples. Differentially enriched and core taxa, however, differed between the skin and water samples. In tilapia skin, Cetobacterium, Paucibacter, Pseudomonas and Comamonadaceae were enriched, whereas, the cyanobacteria Cyanobium, Microcystis and/or Synechocystis, and the diatom Cyclotella, were most prevalent in pond water. Ponds that clustered together according to their water prokaryotic communities also had similar microeukaryotic communities indicating strong environmental influences on prokaryotic and microeukaryotic community structures. While strong site-specific clustering was observed in pond water, the grouping of tilapia skin prokaryotes by pond site was less distinct, suggesting fish microbiota have a greater buffering capacity against environmental influences. The characterised diversity, structure and variance of microbial communities associated with tilapia culture in Malawi provides the baseline for studies on how future intensification practices may lead to microbial dysbiosis and disease onset.

### Funding
This work was funded by the BBSRC/Newton Fund project (BB/N00504X/1). JM was supported by the BBSRC/South West Biosciences Doctoral Training Partnership (BB/M009122/1) with CASE partners WorldFish and Cefas. Sequencing was performed at the Exeter Sequencing Service, using equipment funded by the Wellcome Trust Institutional Strategic Support Fund (WT097835MF), Wellcome Trust Multi-User Equipment Award (WT101650MA) and BBSRC LOLA award (BB/K003240/1). This work was further supported by the CGIAR Research Program on Fish Agri-Food Systems (FISH) led by WorldFish. CC and JN salary was supported by WorldFish. 


### How to regenerate this repository
Original sequencing reads are available from the European Nucleotide Archive under the accession number PRJEB46984.
Additionally, final ASV tables and phyloseq objects used for analyses are avaliable under `data/process`
R (v. 3.6.3) should be located in the user's PATH
Main R packages (not exhaustive):
* `knitr`
* `rmarkdown`
* `tidyverse`
* `DADA2`
* `phyloseq`
* `decontam`
* `ggplot2`
* `corncob`
* `microbiome`
* `pheatmap`
* `factoextra`
* `metacoder`
* `here`
