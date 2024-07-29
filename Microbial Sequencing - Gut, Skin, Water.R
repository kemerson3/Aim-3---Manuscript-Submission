##Microbial Sequencing - Aim 3 - gut, skin, and water samples
library (dada2); packageVersion("dada2") 
library(devtools)
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions
library(phangorn)
library(Biostrings)
library(phyloseq)
library(ggplot2)
library(DECIPHER); packageVersion("DECIPHER")
library(tidyverse)
library(ggtext)
library(vegan)
library(dplyr)
library(abdiv)
library(picante)
library(RColorBrewer)
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

##Tutorial from package designers: https://benjjneb.github.io/dada2/tutorial.html
##Youtube tutorial walkthrough with tree construction: https://www.youtube.com/watch?v=wV5_z7rR6yw&t=2574s

#****set working directory/folder to the one that contains the fastq files*******
setwd("C:/BaseSpace/MISEQ792_Kyle-Emerson_45704-394461778/FASTQ_Generation_2023-08-03_04_45_31Z-685575895 - Gut + Skin")
path <- "C:/BaseSpace/MISEQ792_Kyle-Emerson_45704-394461778/FASTQ_Generation_2023-08-03_04_45_31Z-685575895 - Gut + Skin" 
list.files(path) #list the files within that directory/folder within the path

# Forward and reverse fastq filenames have format: 
# Forward = SAMPLENAME_R1_001.fastq and 
# Reverse = SAMPLENAME_R2_001.fastq
#This pulls each FASTQ based on naming pattern and places as Fwd or Rvs
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
#Removes extra values besides the "B1A-2" & "FA-B1A"
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Visualizing quality profiles 
#This will generate a plot for FORWARD READS (R1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnFs[3:4])
plotQualityProfile(fnFs[5:6])
plotQualityProfile(fnFs[7:8])
plotQualityProfile(fnFs[25:26])
plotQualityProfile(fnFs[27:28])
plotQualityProfile(fnFs[29:30])
plotQualityProfile(fnFs[31:32])
#Changing the numbers will show you the quality profile of different sequences

#In gray-scale is a heat map of the frequency of each quality score at each 
#base position. 
#The mean quality score at each position is shown by the green line, and 
#the quartiles of the quality score distribution by the orange lines. 
#The red line shows the scaled proportion of reads that extend to at least 
#that position (this is more useful for other sequencing technologies, 
#as Illumina reads are typically all the same length, hence the flat red line).

#Generally advise to trim the last few nucleotides to avoid less well-controlled errors. 
#Looks like truncate at 270

#Visualize the quality profile of the reverse reads (R2) within a plot:
plotQualityProfile(fnRs[1:2])
plotQualityProfile(fnRs[3:4])
plotQualityProfile(fnRs[5:6])
plotQualityProfile(fnRs[7:8])
plotQualityProfile(fnRs[25:26])
plotQualityProfile(fnRs[27:28])
plotQualityProfile(fnRs[29:30])
plotQualityProfile(fnRs[31:32])
#scores are a touch lower and crash at 220bp 

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#Using the assigned path above, filtered sequences are named either 
# _F_filt or _R_filt whether forward or reserve reads. 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(270, 220), trimLeft = c(19,20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, matchIDs = TRUE) # On Windows set multithread=FALSE
head(out)
# Trim left is my trimming of primers (19 is forward primer (515F) length)
# 20 is my length of the 806R primer length
# Based on quality profile plots, the only trimming required is that of my forward and reverse primers
# On Windows set multithread=FALSE
# The truncQ = 2 code, truncates the reads at the first instance of a quality 
# score less than or equal to 2. 
# The “matchIDs = true” code is used for the paired-read filtering, to enforce 
# matching between the id-line sequence identifiers of the forward and reverse 
# fastq files. Since this is true, ONLY paired reads that share sequence ID 
# fields are shown in the output. 

#Learn the Error Rates

# The DADA2 algorithm makes use of a parametric error model (err) and every 
#amplicon dataset has a different set of error rates. 
#The learnErrors method learns this error model from the data, 
#by alternating estimation of the error rates and inference of sample 
#composition until they converge on a jointly consistent solution. 
#As in many machine-learning problems, the algorithm must begin with an initial guess, 
#for which the maximum possible error rates in this data are used 
#(the error rates if only the most abundant sequence is correct and all the rest are errors).

#Error model, A<>C<>T<>G changes of nucleic bases
errF <- learnErrors(filtFs, multithread=FALSE) 
#130368145 total bases in 519395 reads from 1 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=FALSE)
#103879000 total bases in 519395 reads from 1 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the 
#machine-learning algorithm. 
#The red line shows the error rates expected under the 
#nominal definition of the Q-score. 
#Here the estimated error rates (black line) are a good fit to the 
#observed rates (points), 
#and the error rates drop with increased quality as expected. 
#Everything looks reasonable and we proceed with confidence.

######################################################
#Sample Inference

#We are now ready to apply the core sample inference algorithm to the filtered 
#and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
#Inspecting the returned dada-class object:
dadaFs[[1]]
print(dadaFs[[1]])
print(dadaFs)
#dada-class: object describing DADA2 denoising results
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40
#, BAND_SIZE = 16
#The DADA2 algorithm inferred 95 sequence variants were inferred from 22878 input unique sequences. There is much more to the dada-class 
#return object than this (see help("dada-class") for some info).

#####################################
##Merge paired reads
#We now merge the forward and reverse reads together to obtain the full denoised sequences. 
#Merging is performed by aligning the denoised forward reads with the 
#reverse-complement of the corresponding denoised reverse reads, and then 
#constructing the merged “contig” sequences. 
#By default, merged sequences are only output if the forward and reverse reads 
#overlap by at least 12 bases, and are identical to each other in the overlap 
#region (but these conditions can be changed via function arguments).

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


#Constructing a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) 
# 55 samples, 4999 unique ASVs
# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]
table(nchar(getSequences(seqtab2)))
# Considerations for your own data: Sequences that are much longer or shorter 
# than expected may be the result of non-specific priming. 
# You can remove non-target-length sequences from your sequence table 
# (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). 
# This is analogous to “cutting a band” in-silico to get amplicons of the targeted length

#Remove chimeras
#The core dada method corrects substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of sequence variants after denoising makes identifying 
#chimeric ASVs simpler than when dealing with fuzzy OTUs. 
#Chimeric sequences are identified if they can be exactly reconstructed by 
#combining a left-segment and a right-segment from two more abundant “parent” sequences.
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", 
                                    multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)
#97% of reads not chimeric

####Removal of Singleton reads
##According to tutorial and video, singleton reads are automatically filtered out
##and are not classified as unique sequences
##Our seqtab.nochim is our current ASV table (unrarefied). Currently, I see
##A few sequences at the end of the table with 1 read. Unsure if they are classified as singletons,
##But the link below is from the author and is what I will use to remove these sequences
## https://github.com/benjjneb/dada2/issues/1519

is1 <- colSums(seqtab.nochim) <= 1
seqtab.nochim1 <- seqtab.nochim[,!is1]

#####Track reads through the pipeline
#As a final check of our progress, we’ll look at the number of reads that 
#made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim1))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
list(track)

##### Assign Taxonomy

#It is common at this point, especially in 16S/18S/ITS amplicon sequencing, 
#to assign taxonomy to the sequence variants. The DADA2 package provides a native 
#implementation of the naive Bayesian classifier method for this purpose. 
#The assignTaxonomy function takes as input a set of sequences to be classified 
#and a training set of reference sequences with known taxonomy, 
#and outputs taxonomic assignments with at least minBoot bootstrap confidence.
#We maintain formatted training fastas for the RDP training set, 
#GreenGenes clustered at 97% identity, and the Silva reference database, 
#and additional trainings fastas suitable for protists and certain specific 
#environments have been contributed. For fungal taxonomy, the General Fasta 
#release files from the UNITE ITS database can be used as is. 
#To follow along, download the silva_nr_v132_train_set.fa.gz file, 
#These files can be downloaded here: https://benjjneb.github.io/dada2/training.html
#and place it in the directory with the fastq files.

file.choose()
taxa <- assignTaxonomy(seqtab.nochim1, "C:\\BaseSpace\\MISEQ792_Kyle-Emerson_45704-394461778\\FASTQ_Generation_2023-08-03_04_45_31Z-685575895 - Gut + Skin\\silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)
#Assigning the taxonomy and adding the species level will take a considerable
#amount of time. Best to run the code and check back later. 

write.table(taxa, file = "Emerson Aim 3 - Taxa - Gut Skin Water.txt", sep="\t")
#Creates a txt file in the path listed of the taxa information. Note:
#this goes into your designated R folder (for me: Aim 3 - Stress and the Microbiota)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#####Making a phylogenetic tree
## Note: this is strictly for phylogenetic analysis
## I.e. Faiths. Will not impact ASVs, Shannon, Bray Curtis, etc
## Information found in the youtube tutorial listed at beginning of doc

# Align Sequences for Phylogeny
# Extract sequences from DADA2 output
sequences <- getSequences(seqtab.nochim1)
names(sequences)<-sequences
# Run sequence alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor = NA)
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
# Change sequence alignment output into a phyDat structure

dm <- dist.ml(phang.align)
# Making a distance matrix

treeNJ <- NJ(dm)
# Making a neighbor joining tree

fit = pml(treeNJ, data = phang.align)
# Fitting for internal maximum likelihood

fitGTR <- update(fit, k =4, inv = 0.2)

fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

### Have my tables (unrarefied) and my tree
theme_set(theme_bw())
file.choose()
meta_full <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Emerson Aim 3 - Microbial Map - Gut, Skin, Water.csv")
#This is my map file. This file is necessary to give R a point of reference
#to merge files downstream. The most important consistency is including the
#sample names that are identical across all files

row.names(meta_full) <- meta_full$Sample_ID
#ID needs to be row names for correct merging of files

meta_full$Microbial_Trtmt = factor(meta_full$Microbial_Trtmt)
meta_full$Stressor = factor(meta_full$Stressor)
meta_full$Trtmt_Combo = factor(meta_full$Trtmt_Combo)
meta_full$Source = factor(meta_full$Source)

ps <- phyloseq(otu_table(seqtab.nochim1, taxa_are_rows=FALSE), 
               tax_table(taxa), sample_data(meta_full),
               phy_tree(fitGTR$tree))

###### Rooting my phylogenetic tree

set.seed(711)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))

##### Creating my phyloseq object for analysis

ps <- merge_phyloseq(ps, map)
#Merge ps object with map
#Creates a phyloseq object that uses sequences with chimeras removed from the 
#sequencing run, data from the meta table, and information on taxa from the downloaded
#Silva file 
#had to change sample data to sample names

taxa_names(ps)

wholetax <- do.call(paste, c(as.data.frame(tax_table(ps))
                             [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")],
                             sep = "__"))  # to distinguish from "_" within tax ranks
#generate a vector containing the full taxonomy path for all ASVs

otu_export <- as.data.frame(otu_table(ps))
tmp <- names(otu_export)
# turn the otu_table into a data.frame

for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}
# paste wholetax and OTU_ids together

names(otu_export) <- names(tmp)
head(otu_export)[5]
#overwrite old names

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
list(ps)
#3657 ASVs (55 samples) with all singletons and chimeric sequences removed

refseqPS <-refseq(ps)
#return non-empty slot names of phyloseq object

getslots.phyloseq(ps)
#we have "otu_table" "tax_table" "sam_data" "refseq" "phy_tree"
#our completed phyloseq object

########### Cleaning up my data in R (https://www.youtube.com/watch?v=e3rKYipvdJo)
#write.table(otu_table(ps), file = "Emerson Aim 3 - gut, water & skin otutable.csv", sep = ",")
#write.table(tax_table(ps), file = "Emerson Aim 3 - gut, water & skin taxatable.csv", sep = ",")
#write.table(refseq(ps), file = "Emerson Aim 3 - gut, water & skin refseqtable.csv", sep = ",")

ps_ra = transform_sample_counts(ps,function(x){x/sum(x)})
#write.table(otu_table(ps_ra), file = "Emerson Aim 3 - gut, water & skin relabundtable.csv", sep = ",")
#This is our relative abundance table. 

file.choose()
otu_counts <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Emerson Aim 3 - gut, water & skin otutable.csv") %>%
  pivot_longer(-Sample_ID, names_to="ASV", values_to="count")

file.choose()
taxonomy <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Emerson Aim 3 - gut, water & skin taxatable.csv") %>%
  filter(Kingdom != 'Archaea') %>%
  filter(Kingdom != "Eukaryota")

#This filter step allowed me to remove any non-bacterial taxa from the df
#fileEncoding="UTF-8-BOM" can be used when Excel -> CSV gives weird symbols in column headers

#Now, I want to join my dataframes together
#meta gut & water
#otu_counts gut & water
#taxonomy gut & water
################ Relative Abundance
# https://www.statology.org/dplyr-remove-rows/#:~:text=You%20can%20use%20the%20following%20basic%20syntax%20to,4%29%29%205%205.%20Remove%20rows%20based%20on%20condition

otu_rel_abund <- inner_join(meta_full, otu_counts, by = "Sample_ID") %>%
  inner_join(., taxonomy, by = "ASV") %>%
  group_by(Sample_ID) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV"),
               names_to="level",
               values_to= "taxon") %>%
  filter(level != 'ASV')

##Now, our data frame only includes bacteria

############### Relative Abundance data frame: Phylum
phyla_abundance <- otu_rel_abund %>%
  select(-Microbial_Trtmt, -Stressor, -Trtmt_Combo) %>%
  filter(level =="Phylum") %>%
  group_by(Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  pivot_wider(names_from = "taxon", values_from = "rel_abund")

#write.table(phyla_abundance, file = "Emerson Aim 3 - gut, skin & water phyla abundance.csv", sep = "," )  

############### Relative Abundance data frame: Genus
genus_abundance <- otu_rel_abund %>%
  select(-Microbial_Trtmt, -Stressor, -Trtmt_Combo) %>%
  filter(level =="Genus") %>%
  group_by(Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  pivot_wider(names_from = "taxon", values_from = "rel_abund")

#write.table(genus_abundance, file = "Emerson Aim 3 - gut, skin & water genus abundance.csv", sep = "," )

#Check to make sure our relative abundances add up to 1
otu_rel_abund %>%
  group_by(Sample_ID) %>%
  summarize(total = sum(rel_abund)) %>%
  print(n=57)
#Success, says 6 but that is because each taxa is counted 6 times due to 
#that being the depth of taxonomic assignment

#######Rarefaction for gut, skin water

rarefy_asv_full <- inner_join(meta_full, otu_counts, by = "Sample_ID") %>%
  inner_join(., taxonomy, by = "ASV") %>%
  group_by(Sample_ID) %>%
  select(Sample_ID, Microbial_Trtmt, Stressor, Trtmt_Combo, ASV, Source, count)
# This above command turned our otu_counts table into a new table
# that has archaea removed and is strictly bacteria

sampling_coverage_full <- rarefy_asv_full %>%
  group_by(Sample_ID) %>%
  summarize(n_seqs = sum(count)) 
#gives us the number of sequences in each sample
#Lowest number in an actual sample is 1124 reads, and we lose our blanks and AA301S and NA301S

rarefied_asv_table_full <- rarefy_asv_full %>%
  group_by(Sample_ID) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 1124) %>%
  select(-n)  %>%
  select(-Stressor, -Microbial_Trtmt, -Trtmt_Combo) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  column_to_rownames("Sample_ID")
# Our new, rarefied (1124) ASV counts table for gut and water samples without archaea

#write.table(rarefied_asv_table_full, file = "Emerson Aim 3 - rarefied asv table gut, skin & water.csv", sep = "," )

####### Beta Diversity - Gut, Skin & Water Samples
##Bray Curtis
#https://www.youtube.com/watch?v=G5Qckqq5Erw create PCoA plot
#https://www.youtube.com/watch?v=xyufizOpc5I nmds plot

distance_matrix_full <- rarefy_asv_full %>%
  group_by(Sample_ID) %>%
  select(-Microbial_Trtmt, -Stressor, -Trtmt_Combo, -Source) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 1124) %>%
  ungroup() %>%
  group_by(ASV) %>%
  select(-n) %>%
  mutate(total = sum(count)) %>%
  filter(total !=0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  as.data.frame()

rownames(distance_matrix_full) <- distance_matrix_full$Sample_ID
distance_matrix_full <- distance_matrix_full[,-1]
distance_matrix_full <- as.matrix(distance_matrix_full)

set.seed(19950406)
dist_full <- avgdist(distance_matrix_full, dmethod = "bray", sample = 1124)
set.seed(17)
nmds_full <- metaMDS(dist_full)

metadata_nmds_full <- nmds_full$points %>%
  as_tibble(rownames = "Sample_ID") %>%
  inner_join(., meta_full, by = "Sample_ID") 

Pond_Water <- c("Natural", "Autoclaved")
Stress <- c("Vehicle Control", "Predator Cues", "CORT")

metadata_nmds_full %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Source)) +
  geom_point(aes(shape = Source, size = 5)) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 14))  +
  theme(axis.title = element_text(face = "bold",size = 16)) +
  labs(x = "Bray-Curtis NMDS 1", y = "Bray-Curtis NMDS 2") +
  scale_color_manual(labels = c("Gut", "Skin", "Water"),
                     values = c("pink", "green", "skyblue"),
                     name = "Source") +
  stat_ellipse() +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 16)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  guides(size = "none")

adonis2(dist_full~metadata_nmds_full$Source)

# adonis2(formula = dist_full ~ metadata_nmds_full$Source)
#                           Df SumOfSqs      R2      F  Pr(>F)    
# metadata_nmds_full$Source  2   4.6056 0.29186 9.6855  0.001 ***
# Residual                  47  11.1747 0.70814                  
# Total                     49  15.7803 1.00000 



pairwise.adonis(dist_full, metadata_nmds_full$Source)
#         pairs  Df  SumsOfSqs  F.Model        R2 p.value p.adjusted sig
# 1 Gut vs Skin   1  2.337956 10.339593 0.1975482   0.001      0.003   *
# 2 Gut vs Water  1  2.402528 10.412713 0.2710747   0.001      0.003   *
# 3 Skin vs Water 1  2.128310  7.991117 0.2497918   0.001      0.003   *

bd_gws <- betadisper(dist_full, metadata_nmds_full$Source)
anova(bd_gutmicro)
# Response: Distances
#           Df  Sum Sq    Mean Sq F value Pr(>F)
# Groups     1 0.006862 0.0068622  1.1534 0.2945
# Residuals 22 0.130885 0.0059493

######### Gut vs. Skin - Bray Curtis 

rarefy_asv_gs <- rarefy_asv_full %>%
  filter(Source != 'Water')

# This above command turned our otu_counts table into a new table containing just gut and skin samples
# that has archaea removed and is strictly bacteria

sampling_coverage_gs <- rarefy_asv_gs %>%
  group_by(Sample_ID) %>%
  summarize(n_seqs = sum(count)) 
#gives us the number of sequences in each sample
#Lowest number in an actual sample is 11880 reads, and we lose our blanks and AA301S and NA301S

rarefied_asv_table_gs <- rarefy_asv_gs %>%
  group_by(Sample_ID) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 11880) %>%
  select(-n)  %>%
  select(-Stressor, -Microbial_Trtmt, -Trtmt_Combo) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  column_to_rownames("Sample_ID")
#Our new, rarefied (11880) ASV counts table for gut and skin samples without archaea

#write.table(rarefied_asv_table_gs, file = "Emerson Aim 3 - rarefied asv table gut & skin.csv", sep = "," )

distance_matrix_gs <- rarefy_asv_gs %>%
  group_by(Sample_ID) %>%
  select(-Microbial_Trtmt, -Stressor, -Trtmt_Combo, -Source) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 11880) %>%
  ungroup() %>%
  group_by(ASV) %>%
  select(-n) %>%
  mutate(total = sum(count)) %>%
  filter(total !=0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  as.data.frame()

rownames(distance_matrix_gs) <- distance_matrix_gs$Sample_ID
distance_matrix_gs <- distance_matrix_gs[,-1]
distance_matrix_gs <- as.matrix(distance_matrix_gs)

set.seed(19950406)
dist_gs <- avgdist(distance_matrix_gs, dmethod = "bray", sample = 11880)
set.seed(17)
nmds_gs <- metaMDS(dist_gs)

meta_gs <- meta_full %>%
  filter(Source != 'Water')

metadata_nmds_gs <- nmds_gs$points %>%
  as_tibble(rownames = "Sample_ID") %>%
  inner_join(., meta_gs, by = "Sample_ID") 

Pond_Water <- c("Natural", "Autoclaved")
Stress <- c("Vehicle Control", "Predator Cues", "CORT")

metadata_nmds_gs %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Source)) +
  geom_point(aes(shape = Source, size = 3)) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 14))  +
  theme(axis.title = element_text(face = "bold",size = 16)) +
  labs(x = "Bray-Curtis NMDS 1", y = "Bray-Curtis NMDS 2") +
  scale_color_manual(labels = c("Gut", "Skin"),
                     values = c("pink", "green"),
                     name = "Source") +
  stat_ellipse() +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text( size = 12))

adonis2(dist_gs~metadata_nmds_gs$Source)

# adonis2(formula = dist_gs ~ metadata_nmds_gs$Source)
#                         Df  SumOfSqs      R2    F    Pr(>F)    
# metadata_nmds_gs$Source  1   2.3524 0.20112 10.574  0.001 ***
# Residual                42   9.3439 0.79888                  
# Total                   43  11.6963 1.00000                  


######## Gut Samples: Bray Curtis 

rarefy_asv_gut <- rarefy_asv_full %>%
  filter(Source == 'Gut')

sampling_coverage_gut <- rarefy_asv_gut %>%
  group_by(Sample_ID) %>%
  summarize(n_seqs = sum(count))
#66378 is our rarefaction point

rarefied_asv_table_gut <- rarefy_asv_gut %>%
  group_by(Sample_ID) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 66378) %>%
  select(-n)  %>%
  select(-Stressor, -Microbial_Trtmt, -Trtmt_Combo) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  column_to_rownames("Sample_ID")
#Our new, rarefied (66378) ASV counts table for gut samples without archaea

#write.table(rarefied_asv_table_gut, file = "Emerson Aim 3 - rarefied asv table - gut.csv", sep = "," )

distance_matrix_gut <- rarefy_asv_gut %>%
  group_by(Sample_ID) %>%
  select(-Microbial_Trtmt, -Stressor, -Trtmt_Combo, -Source) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 66378) %>%
  ungroup() %>%
  group_by(ASV) %>%
  select(-n) %>%
  mutate(total = sum(count)) %>%
  filter(total !=0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  as.data.frame()

rownames(distance_matrix_gut) <- distance_matrix_gut$Sample_ID
distance_matrix_gut <- distance_matrix_gut[,-1]
distance_matrix_gut <- as.matrix(distance_matrix_gut)

set.seed(19950406)
dist_gut <- avgdist(distance_matrix_gut, dmethod = "bray", sample = 66378)
set.seed(17)
nmds_gut <- metaMDS(dist_gut)

meta_gut <- meta_full %>%
  filter(Source == 'Gut')

metadata_nmds_gut <- nmds_gut$points %>%
  as_tibble(rownames = "Sample_ID") %>%
  inner_join(., meta_gut, by = "Sample_ID") 

Pond_Water <- c("Natural", "Autoclaved")
Stress <- c("Vehicle Control", "Predator Cues", "CORT")

#Bray - Curtis Beta Diversity  - Gut Samples Plot
metadata_nmds_gut %>%
  ggplot(aes(x = MDS1, y = MDS2, color = factor(Stressor, level=c('B', 'A','C')))) +
  geom_point(aes(shape = Microbial_Trtmt, size = 5)) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 14))  +
  theme(axis.title = element_text(face = "bold",size = 16)) +
  labs(x = "Bray-Curtis NMDS 1", y = "Bray-Curtis NMDS 2") +
  scale_color_manual(labels = c("Vehicle", "Predator", "CORT"),
                     values = c("gray81", "darkred", "deepskyblue1"),
                     name = "Stressor") +
  scale_shape_manual(labels = c("Natural", "Autoclaved"),
                     values = c(16,18),
                     name = "Pond Water") +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 16)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  guides(size = "none")


adonis2(dist_gut~metadata_nmds_gut$Microbial_Trtmt*metadata_nmds_gut$Stressor)
# adonis2(formula = dist_gs ~ metadata_nmds_gs$Source)
#                              Df  SumOfSqs      R2    F    Pr(>F)    
# metadata_nmds_gut$MicroTrtmt 1   1.25          0.266 9.79  0.001 ***
# metadata_nmds_gut$Stressor   2   0.466         0.099 1.83  0.033 *
# Interaction                  2   0.676         0.14  2.65  0.002 **
# Residual                     18  2.29          0.49      
# Total                        23  4.68          1.00

bd_gutmicro <- betadisper(dist_gut, metadata_nmds_gut$Microbial_Trtmt)
anova(bd_gutmicro)
# Response: Distances
#           Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     1 0.006862 0.0068622  1.1534 0.2945
# Residuals 22 0.130885 0.0059493 

bd_gutstress <- betadisper(dist_gut, metadata_nmds_gut$Stressor)
anova(bd_gutstress)

# Response: Distances
#           Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     2 0.000745 0.0003723  0.0745 0.9285
# Residuals 21 0.104965 0.0049983

#https://github.com/pmartinezarbizu/pairwiseAdonis
pairwise.adonis(dist_gut, metadata_nmds_gut$Stressor)
#    pairs Df SumsOfSqs   F.Model        R2 p.value p.adjusted sig
# 1 A vs B  1 0.2033785 0.9925921 0.0662055   0.423       1.00    
# 2 A vs C  1 0.3166640 1.5867374 0.1018005   0.150       0.45    
# 3 B vs C  1 0.1796172 0.9093417 0.0609914   0.454       1.00 


######## Skin Samples: Bray Curtis 

rarefy_asv_skin <- rarefy_asv_full %>%
  filter(Source == 'Skin')

sampling_coverage_skin <- rarefy_asv_skin %>%
  group_by(Sample_ID) %>%
  summarize(n_seqs = sum(count))
#11880 is our rarefaction point

rarefied_asv_table_skin <- rarefy_asv_skin %>%
  group_by(Sample_ID) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 11880) %>%
  select(-n)  %>%
  select(-Stressor, -Microbial_Trtmt, -Trtmt_Combo) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  column_to_rownames("Sample_ID")
#Our new, rarefied (11880) ASV counts table for skin samples without archaea

#write.table(rarefied_asv_table_skin, file = "Emerson Aim 3 - rarefied asv table - skin.csv", sep = "," )

distance_matrix_skin <- rarefy_asv_skin %>%
  group_by(Sample_ID) %>%
  select(-Microbial_Trtmt, -Stressor, -Trtmt_Combo, -Source) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 11880) %>%
  ungroup() %>%
  group_by(ASV) %>%
  select(-n) %>%
  mutate(total = sum(count)) %>%
  filter(total !=0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  as.data.frame()

rownames(distance_matrix_skin) <- distance_matrix_skin$Sample_ID
distance_matrix_skin <- distance_matrix_skin[,-1]
distance_matrix_skin <- as.matrix(distance_matrix_skin)

set.seed(19950406)
dist_skin <- avgdist(distance_matrix_skin, dmethod = "bray", sample = 11880)
set.seed(17)
nmds_skin <- metaMDS(dist_skin)

meta_skin <- meta_full %>%
  filter(Source == 'Skin')

metadata_nmds_skin <- nmds_skin$points %>%
  as_tibble(rownames = "Sample_ID") %>%
  inner_join(., meta_skin, by = "Sample_ID") 

#Bray - Curtis Beta Diversity  - Skin Samples Plot
metadata_nmds_skin %>%
  ggplot(aes(x = MDS1, y = MDS2, color = factor(Stressor, level=c('B', 'A','C')))) +
  geom_point(aes(shape = Microbial_Trtmt, size = 5)) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 14))  +
  theme(axis.title = element_text(face = "bold",size = 16)) +
  labs(x = "Bray-Curtis NMDS 1", y = "Bray-Curtis NMDS 2") +
  scale_color_manual(labels = c("Vehicle", "Predator", "CORT"),
                     values = c("gray81", "darkred", "deepskyblue1"),
                     name = "Stressor") +
  scale_shape_manual(labels = c("Natural", "Autoclaved"),
                     values = c(16,18),
                     name = "Pond Water") +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 16)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  guides(size = "none")

adonis2(dist_skin~metadata_nmds_skin$Microbial_Trtmt*metadata_nmds_skin$Stressor)
# adonis2(formula = dist_gs ~ metadata_nmds_gs$Source)
#                              Df  SumOfSqs         R2    F    Pr(>F)    
# metadata_nmds_gut$MicroTrtmt 1   0.56          0.121  2.47  0.006 **
# metadata_nmds_gut$Stressor   2   0.46          0.099  1.00  0.436 
# Interaction                  2   0.45          0.097  0.98  0.460 
# Residual                     14  3.18          0.68     
# Total                        29  4.65          1.00

bd_skinmicro <- betadisper(dist_skin, metadata_nmds_skin$Microbial_Trtmt)
anova(bd_skinmicro)
# Response: Distances
#           Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     1 0.005895 0.0058954  0.6731 0.4227
# Residuals 18 0.157642 0.0087579  

bd_skinstress <- betadisper(dist_skin, metadata_nmds_skin$Stressor)
anova(bd_skinstress) 
# Response: Distances
#           Df   Sum Sq   Mean Sq F value Pr(>F)
# Groups     2 0.003397 0.0016984  0.1908  0.828
# Residuals 17 0.151313 0.0089008 

pairwise.adonis(dist_skin, metadata_nmds_skin$Stressor)
# pairs   Df  SumsOfSqs   F.Model         R2  p.value  p.adjusted sig
# 1 A vs B  1 0.2239956 0.9151861 0.07680838   0.513          1    
# 2 A vs C  1 0.2520238 1.0257150 0.09302934   0.422          1    
# 3 B vs C  1 0.1977131 0.7864382 0.05704433   0.674          1


############### Relative Abundance Plots - Water Samples
## https://www.youtube.com/watch?v=NVym44SdcaE&t=309s 
## Making stacked barcharts for relative abundance

otu_rel_abund_water <- otu_rel_abund %>%
  filter(Source != 'Skin') %>%
  filter(Source != 'Gut') %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(level =="Phylum") 

otu_rel_abund_water %>%
  group_by(Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Sample_ID, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  ggplot(aes(x = Sample_ID, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_discrete(name = NULL) +
  labs(x = "Pond Water Samples", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(10, "pt"))

taxon_rel_abund_water <- otu_rel_abund_water %>%
  filter(level =="Phylum") %>%
  group_by(Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Sample_ID, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool_water <- taxon_rel_abund_water %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop") 

NPW_names <- c("NPW 1", "NPW 2", "NPW 3", "NPW 4", "NPW 5",
                     "NPW 6")

inner_join(taxon_rel_abund_water, taxon_pool_water, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Sample_ID, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = Sample_ID, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Acidobacteriota", "Actinobacteriota",
                               "Bacteroidota", "Bdellovibrionota",
                               "Chloroflexi", "Desulfobacterota", "Firmicutes",
                               "Myxococcota",
                               "Patescibacteria", "Planctomycetota",
                               "Proteobacteria", "Verrucomicrobiota", "Other"),
                    values = c("darkred", "orangered3","orangered",
                               "lightyellow4","lightyellow3","lightyellow",
                               "orange","yellow","skyblue",
                               "skyblue3","cornflowerblue","darkslateblue", "orchid")) +
  scale_x_discrete(labels = NPW_names) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Pond Water Samples", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(face = "bold", size = 13, angle = 45, hjust = 1))  +
  theme(axis.text.y = element_text(face = "bold", size = 14))  +
  theme(axis.title = element_text(face = "bold", size = 16)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))


############### Relative Abundances Plots - Gut Samples
## https://www.youtube.com/watch?v=NVym44SdcaE&t=309s 
## Making stacked barcharts for relative abundance

Pond_Water_Rel <- c("Natural", "Autoclaved", "Pond Water")

otu_rel_abund_gw <- otu_rel_abund %>%
  filter(Source != 'Skin') %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(level =="Phylum") 

otu_rel_abund_gw %>%
  group_by(Trtmt_Combo, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Trtmt_Combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  ggplot(aes(x = Trtmt_Combo, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_discrete(name = NULL) +
  labs(x = "Treatment Combinations", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(10, "pt"))

taxon_rel_abund_gw <- otu_rel_abund_gw %>%
  filter(level =="Phylum") %>%
  group_by(Trtmt_Combo, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Trtmt_Combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool_gw <- taxon_rel_abund_gw %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop") 



##### Relative Abundance - Treatment Combination (gut samples)

Treatment_Combo <- c("Nat - Ctrl", "Nat - Pred", "Nat - CORT", "AC - Ctrl", "AC - Pred",
                     "AC - CORT", "Pond Water")

a <- ifelse(Treatment_Combo == 7, "blue", "black")

inner_join(taxon_rel_abund_gw, taxon_pool_gw, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Trtmt_Combo, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = factor(Trtmt_Combo, level = c('2', '1', '3', '5', '4', '6', '7')), y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Acidobacteriota", "Actinobacteriota",
                               "Bacteroidota", "Firmicutes",
                               "Fusobacteriota", "Myxococcota",
                               "Patescibacteria", "Planctomycetota",
                               "Proteobacteria", "Verrucomicrobiota", "Other"),
                    values = c(brewer.pal(11, "RdYlBu"))) +
  scale_x_discrete(labels = Treatment_Combo) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Treatment Combination", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(face = "bold", colour = c("black", "black", "black", "black", "black", "black", "blue"), size = 13, angle = 45, hjust = 1))  +
  theme(axis.text.y = element_text(face = "bold", size = 14))  +
  theme(axis.title = element_text(face = "bold", size = 16)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))

#### Relative Abundance - Microbial Treatment (gut samples)

taxon_rel_abund_gutmicro <- otu_rel_abund_gw %>%
  filter(level =="Phylum") %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  group_by(Microbial_Trtmt, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool_gutmicro <- taxon_rel_abund_gutmicro %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop")

inner_join(taxon_rel_abund_gutmicro, taxon_pool_gutmicro, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = Microbial_Trtmt, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Acidobacteriota", "Actinobacteriota",
                               "Bacteroidota", "Firmicutes",
                               "Fusobacteriota", "Myxococcota",
                               "Patescibacteria", "Planctomycetota",
                               "Proteobacteria", "Verrucomicrobiota", "Other"),
                    values = c(brewer.pal(11, "RdYlBu"))) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Microbial Environment", y = "Mean Relative Abundance (%)") +
  scale_x_discrete(labels = Pond_Water_Rel) +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(face = "bold", size = 14))  +
  theme(axis.text.y = element_text(face = "bold", size = 14))  +
  theme(axis.title = element_text(face = "bold", size = 16)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))

#### Relative Abundance - Stressor Trtmt (gut samples)

taxon_rel_abund_gutstressor <- otu_rel_abund_gw %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle061-NPW1') %>%
  filter(Sample_ID != 'Kyle061-NPW2') %>%
  filter(Sample_ID != 'Kyle061-NPW3') %>%
  filter(Sample_ID != 'Kyle061-NPW4') %>%
  filter(Sample_ID != 'Kyle061-NPW5') %>%
  filter(Sample_ID != 'Kyle061-NPW6') %>%
  filter(level =="Phylum") %>%
  group_by(Stressor, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Stressor, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool_gutstressor <- taxon_rel_abund_gutstressor %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop")

Stressor_Rel <- c("Control", "Predator Cue", "CORT", "Pond Water")

inner_join(taxon_rel_abund_gutstressor, taxon_pool_gutstressor, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Stressor, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = factor(Stressor, level = c( 'B', 'A', 'C')), y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Acidobacteriota", "Actinobacteriota",
                               "Bacteroidota", "Desulfobacterota",
                               "Firmicutes",
                               "Fusobacteriota", "Myxococcota",
                               "Patescibacteria", "Planctomycetota",
                               "Proteobacteria", "Verrucomicrobiota", "Other"),
                    values = c(brewer.pal(11, "RdYlBu"), "gray")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = Stressor_Rel) +
  labs(x = "Stressor", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(face = "bold", size = 13))  +
  theme(axis.title = element_text(face = "bold", size = 15)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))

############# Relative Abundances Plots - Skin Samples
## https://www.youtube.com/watch?v=NVym44SdcaE&t=309s 
## Making stacked barcharts for relative abundance

Pond_Water_Rel <- c("Natural", "Autoclaved", "Pond Water")

otu_rel_abund_sw <- otu_rel_abund %>%
  filter(Source != 'Gut') %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle039-AA301S') %>%
  filter(level =="Phylum") 

otu_rel_abund_sw %>%
  group_by(Trtmt_Combo, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Trtmt_Combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  ggplot(aes(x = Trtmt_Combo, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_discrete(name = NULL) +
  labs(x = "Treatment Combinations", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(10, "pt"))

taxon_rel_abund_sw <- otu_rel_abund_sw %>%
  filter(level =="Phylum") %>%
  group_by(Trtmt_Combo, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Trtmt_Combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool_sw <- taxon_rel_abund_sw %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop") 

##### Relative Abundance - Treatment Combination (skin samples)

Treatment_Combo <- c("Nat - Ctrl", "Nat - Pred", "Nat - CORT", "AC - Ctrl", "AC - Pred",
                     "AC - CORT", "Pond Water")

inner_join(taxon_rel_abund_sw, taxon_pool_sw, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Trtmt_Combo, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = factor(Trtmt_Combo, level = c('2', '1', '3', '5', '4', '6', '7')), y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Acidobacteriota", "Actinobacteriota",
                               "Bacteroidota", "Firmicutes",
                               "Fusobacteriota", "Myxococcota",
                               "Patescibacteria", "Planctomycetota",
                               "Proteobacteria", "Verrucomicrobiota", "Other"),
                    values = c(brewer.pal(11, "RdYlBu"))) +
  scale_x_discrete(labels = Treatment_Combo) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Treatment Combination", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(face = "bold", colour = c("black", "black", "black", "black", "black", "black", "blue"), size = 13, angle = 45, hjust = 1))+  
  theme(axis.text.y = element_text(face = "bold", size = 14))  +
  theme(axis.title = element_text(face = "bold", size = 16)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))

#### Relative Abundance - Microbial Treatment (gut samples)

taxon_rel_abund_skinmicro <- otu_rel_abund_sw %>%
  filter(level =="Phylum") %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle039-AA301S') %>%
  group_by(Microbial_Trtmt, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool_skinmicro <- taxon_rel_abund_skinmicro %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop")

inner_join(taxon_rel_abund_skinmicro, taxon_pool_skinmicro, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = Microbial_Trtmt, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Acidobacteriota", "Actinobacteriota",
                               "Bacteroidota", "Firmicutes",
                               "Fusobacteriota", "Myxococcota",
                               "Patescibacteria", "Planctomycetota",
                               "Proteobacteria", "Verrucomicrobiota", "Other"),
                    values = c(brewer.pal(11, "RdYlBu"))) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Microbial Environment", y = "Mean Relative Abundance (%)") +
  scale_x_discrete(labels = Pond_Water_Rel) +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(face = "bold", size = 14))  +
  theme(axis.text.y = element_text(face = "bold", size = 14))  +
  theme(axis.title = element_text(face = "bold", size = 16)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))

#### Relative Abundance - Stressor Trtmt (skin samples)

taxon_rel_abund_skinstressor <- otu_rel_abund_sw %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle061-NPW1') %>%
  filter(Sample_ID != 'Kyle061-NPW2') %>%
  filter(Sample_ID != 'Kyle061-NPW3') %>%
  filter(Sample_ID != 'Kyle061-NPW4') %>%
  filter(Sample_ID != 'Kyle061-NPW5') %>%
  filter(Sample_ID != 'Kyle061-NPW6') %>%
  filter(Sample_ID != 'Kyle039-AA301S') %>%
  filter(level =="Phylum") %>%
  group_by(Stressor, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Stressor, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool_skinstressor <- taxon_rel_abund_skinstressor %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop")

Stressor_Rel <- c("Control", "Predator Cue", "CORT", "Pond Water")

inner_join(taxon_rel_abund_skinstressor, taxon_pool_skinstressor, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Stressor, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = factor(Stressor, level = c( 'B', 'A', 'C')), y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Acidobacteriota", "Actinobacteriota",
                               "Bacteroidota", "Desulfobacterota",
                               "Firmicutes",
                               "Fusobacteriota", "Myxococcota",
                               "Patescibacteria", "Planctomycetota",
                               "Proteobacteria", "Verrucomicrobiota", "Other"),
                    values = c(brewer.pal(11, "RdYlBu"), "gray")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = Stressor_Rel) +
  labs(x = "Stressor", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(face = "bold", size = 13))  +
  theme(axis.title = element_text(face = "bold", size = 15)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))

############# Alpha Diversity Metrics - Full

## No. Observed ASVs - Gut, Skin, & Water

observed_asvs_full <- rarefied_asv_table_full[,-1] %>%
  rarefy(sample = 1124) %>%
  as_tibble(rownames = "Sample_ID") %>%
  select(Sample_ID, Obs_ASVs=value)

#write.table(observed_asvs_full, file = "Emerson Aim 3 - no. observed asvs full.csv", sep = "," )

## Shannon Diversity Index

rarefied_asv_table_full[,-1] %>%
  rrarefy(sample = 1124) %>%
  diversity()

#Now, we have shannon diversity values for our samples
#But, we want to repeat this to get the average
#So we are going to repeat this process, average that output for our rarefied shannon

shannon_iteration <- function(){
  
  rarefied_asv_table_full[,-1] %>%
    rrarefy(sample = 1124) %>%
    diversity()
  
}

rarefied_shannon_full <- replicate (100, shannon_iteration()) %>% 
  as_tibble(rownames = "Sample_ID", .name_repair = "unique" ) %>%
  pivot_longer(-Sample_ID) %>%
  group_by(Sample_ID) %>%
  summarize(shannon = mean(value))

#write.table(rarefied_shannon_full, file = "Emerson Aim 3 - shannon diversity full.csv", sep = "," )

## Faiths Phylogenetic Diversity
# #https://search.r-project.org/CRAN/refmans/abdiv/html/faith_pd.html

faiths.df.full <- rarefied_asv_table_full[,-1]

faiths.df.full <- as.matrix(faiths.df.full)
#Turned our rarefied asv table that we have been using into a matrix

Faith_PD_full <- pd(faiths.df.full, phy_tree(ps), include.root = TRUE) 
#These are our faiths phylogenetic diversity values!
#SR represents species richness values, but we will keep our other ones

#write.table(Faith_PD_full, file = "Emerson Aim 3 - Faiths PD Full.csv", sep = "," )


########## Alpha Diversity Metrics - Gut
## No. Observed ASVs - Gut

observed_asvs_gut <- rarefied_asv_table_gut[,-1] %>%
  rarefy(sample = 66378) %>%
  as_tibble(rownames = "Sample_ID") %>%
  select(Sample_ID, Obs_ASVs=value)

#write.table(observed_asvs_gut, file = "Emerson Aim 3 - no. observed gut asvs.csv", sep = "," )

## Shannon Diversity Index

rarefied_asv_table_gut[,-1] %>%
  rrarefy(sample = 66378) %>%
  diversity()

#Now, we have shannon diversity values for our samples
#But, we want to repeat this to get the average
#So we are going to repeat this process, average that output for our rarefied shannon



shannon_iteration <- function(){
  
  rarefied_asv_table_gut[,-1] %>%
    rrarefy(sample = 66378) %>%
    diversity()
  
}

rarefied_shannon_gut <- replicate (100, shannon_iteration()) %>% 
  as_tibble(rownames = "Sample_ID", .name_repair = "unique" ) %>%
  pivot_longer(-Sample_ID) %>%
  group_by(Sample_ID) %>%
  summarize(shannon = mean(value))

#write.table(rarefied_shannon_gut, file = "Emerson Aim 3 - shannon gut diversity.csv", sep = "," )

## Faiths Phylogenetic Diversity
# #https://search.r-project.org/CRAN/refmans/abdiv/html/faith_pd.html

faiths.df.gut <- rarefied_asv_table_gut[,-1]

faiths.df.gut <- as.matrix(faiths.df.gut)
#Turned our rarefied asv table that we have been using into a matrix

Faith_PD_gut <- pd(faiths.df.gut, phy_tree(ps), include.root = TRUE) 
#These are our faiths phylogenetic diversity values!
#SR represents species richness values, but we will keep our other ones

#write.table(Faith_PD_gut, file = "Emerson Aim 3 - Final FaithsPD.csv", sep = "," )

########## Alpha Diversity Metrics - Skin
## No. Observed ASVs - Skin

observed_asvs_skin <- rarefied_asv_table_skin[,-1] %>%
  rarefy(sample = 11880) %>%
  as_tibble(rownames = "Sample_ID") %>%
  select(Sample_ID, Obs_ASVs=value)

#write.table(observed_asvs_skin, file = "Emerson Aim 3 - no. observed skin asvs.csv", sep = "," )

## Shannon Diversity Index

rarefied_asv_table_skin[,-1] %>%
  rrarefy(sample = 11880) %>%
  diversity()

#Now, we have shannon diversity values for our samples
#But, we want to repeat this to get the average
#So we are going to repeat this process, average that output for our rarefied shannon

shannon_iteration_skin <- function(){
  
  rarefied_asv_table_skin[,-1] %>%
    rrarefy(sample = 11880) %>%
    diversity()
  
}

rarefied_shannon_skin <- replicate (100, shannon_iteration_skin()) %>% 
  as_tibble(rownames = "Sample_ID", .name_repair = "unique" ) %>%
  pivot_longer(-Sample_ID) %>%
  group_by(Sample_ID) %>%
  summarize(shannon = mean(value))

#write.table(rarefied_shannon_skin, file = "Emerson Aim 3 - shannon skin diversity.csv", sep = "," )

## Faiths Phylogenetic Diversity
# #https://search.r-project.org/CRAN/refmans/abdiv/html/faith_pd.html

faiths.df.skin <- rarefied_asv_table_skin[,-1]

faiths.df.skin <- as.matrix(faiths.df.skin)
#Turned our rarefied asv table that we have been using into a matrix

Faith_PD_skin <- pd(faiths.df.skin, phy_tree(ps), include.root = TRUE) 
#These are our faiths phylogenetic diversity values!
#SR represents species richness values, but we will keep our other ones

#write.table(Faith_PD_skin, file = "Emerson Aim 3 - Final skin FaithsPD.csv", sep = "," )
