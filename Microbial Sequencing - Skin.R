##Microbial Sequencing - Aim 3
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

##Tutorial from package designers: https://benjjneb.github.io/dada2/tutorial.html
##Youtube tutorial walkthrough with tree construction: https://www.youtube.com/watch?v=wV5_z7rR6yw&t=2574s

#****set working directory/folder to the one that contains the fastq files*******
setwd("C:/BaseSpace/MISEQ792_Kyle-Emerson_45704-394461778/FASTQ_Generation_2023-08-03_04_45_31Z-685575895 - Skin")
path <- "C:/BaseSpace/MISEQ792_Kyle-Emerson_45704-394461778/FASTQ_Generation_2023-08-03_04_45_31Z-685575895 - Skin" 
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
plotQualityProfile(fnFs[16:17])
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
# Looks like NA301 and AA301 have low reads, will probably remove them

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
#106352716 total bases in 423716 reads from 6 samples will be used for learning the error rates.
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
# 31 samples, 3276 unique ASVs
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
#99.3% of reads not chimeric

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
taxa <- assignTaxonomy(seqtab.nochim1, "C:\\BaseSpace\\MISEQ792_Kyle-Emerson_45704-394461778\\FASTQ_Generation_2023-08-03_04_45_31Z-685575895 - Skin\\silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)
#Assigning the taxonomy and adding the species level will take a considerable
#amount of time. Best to run the code and check back later. 

write.table(taxa, file = "Emerson Aim 3 - Skin Taxa.txt", sep="\t")
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
meta_skinwater <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Emerson Aim 3 - Microbial Map - Skin & Water.csv")
#This is my map file. This file is necessary to give R a point of reference
#to merge files downstream. The most important consistency is including the
#sample names that are identical across all files

row.names(meta_skinwater) <- meta_skinwater$Sample_ID
#ID needs to be row names for correct merging of files

meta_skinwater$Microbial_Trtmt = factor(meta_skinwater$Microbial_Trtmt)
meta_skinwater$Stressor = factor(meta_skinwater$Stressor)
meta_skinwater$Trtmt_Combo = factor(meta_skinwater$Trtmt_Combo)
meta_skinwater$Source = factor(meta_skinwater$Source)


ps <- phyloseq(otu_table(seqtab.nochim1, taxa_are_rows=FALSE), 
               tax_table(taxa), sample_data(meta_skinwater),
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
#3015 ASVs with all singletons and chimeric sequences removed

refseqPS <-refseq(ps)
#return non-empty slot names of phyloseq object

getslots.phyloseq(ps)
#we have "otu_table" "tax_table" "sam_data" "refseq" "phy_tree"
#our completed phyloseq object

########### Cleaning up my data in R (https://www.youtube.com/watch?v=e3rKYipvdJo)
write.table(otu_table(ps), file = "Emerson Aim 3 - skin & water otutable.csv", sep = ",")
write.table(tax_table(ps), file = "Emerson Aim 3 - skin & water taxatable.csv", sep = ",")
write.table(refseq(ps), file = "Emerson Aim 3 - skin & water refseqtable.csv", sep = ",")

ps_ra = transform_sample_counts(ps,function(x){x/sum(x)})
write.table(otu_table(ps_ra), file = "Emerson Aim 3 - skin & water relabundtable.csv", sep = ",")
#This is our relative abundance table. 

file.choose()
otu_counts <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Microbiome- Skin and Water Excel Files\\Emerson Aim 3 - skin & water otutable.csv") %>%
  pivot_longer(-Sample_ID, names_to="ASV", values_to="count")

taxonomy <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Microbiome- Skin and Water Excel Files\\Emerson Aim 3 - skin & water taxatable.csv", fileEncoding="UTF-8-BOM") %>%
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

otu_rel_abund <- inner_join(meta_skinwater, otu_counts, by = "Sample_ID") %>%
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

#write.table(phyla_abundance, file = "Emerson Aim 3 - skin & water phyla abundance.csv", sep = "," )  

############### Relative Abundance data frame: Genus
genus_abundance <- otu_rel_abund %>%
  select(-Microbial_Trtmt, -Stressor, -Trtmt_Combo) %>%
  filter(level =="Genus") %>%
  group_by(Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  pivot_wider(names_from = "taxon", values_from = "rel_abund")

#write.table(genus_abundance, file = "Emerson Aim 3 - skin & water genus abundance.csv", sep = "," )

#Check to make sure our relative abundances add up to 1
otu_rel_abund %>%
  group_by(Sample_ID) %>%
  summarize(total = sum(rel_abund)) %>%
  print(n=40)
#Success, says 6 but that is because each taxa is counted 6 times due to 
#that being the depth of taxonomic assignment

## https://www.youtube.com/watch?v=NVym44SdcaE&t=309s 
## Making stacked barcharts for relative abundance

Pond_Water <- c("Natural", "Autoclaved")

otu_rel_abund %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle039-AA301S') %>%
  filter(level =="Phylum") %>%
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


taxon_rel_abund <- otu_rel_abund %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle039-AA301S') %>%
  filter(level =="Phylum") %>%
  group_by(Trtmt_Combo, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Trtmt_Combo, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop") 


##### Relative Abundance - Treatment Combination
Treatment_Combo <- c("Nat - Ctrl", "Nat - Pred", "Nat - CORT", "AC - Ctrl", "AC - Pred",
                     "AC - CORT", "Pond Water")

inner_join(taxon_rel_abund, taxon_pool, by = "taxon") %>%
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
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels = Treatment_Combo) +
  labs(x = "Treatment Combination", y = "Mean Relative Abundance (%)") +
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

#### Relative Abundance - Mircobial Treatment

Pond_Water_Rel <- c("Natural", "Autoclaved", "Pond Water")

taxon_rel_abund_micro <- otu_rel_abund %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle039-AA301S') %>%
  filter(level =="Phylum") %>%
  group_by(Microbial_Trtmt, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool_micro <- taxon_rel_abund_micro %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop")

inner_join(taxon_rel_abund_micro, taxon_pool_micro, by = "taxon") %>%
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

#### Relative Abundance - Stressor Trtmt

Stressor_Rel <- c("Control", "Predator Cues", "CORT", "Pond Water")

taxon_rel_abund_stressor <- otu_rel_abund %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle039-AA301S') %>%
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

taxon_pool_stressor <- taxon_rel_abund_stressor %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop")

inner_join(taxon_rel_abund_stressor, taxon_pool_stressor, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Stressor, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = Stressor, y=mean_rel_abund, fill=taxon)) +
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


#######Rarefaction for skin and water

rarefy_asv_skinwater <- inner_join(meta_skinwater, otu_counts, by = "Sample_ID") %>%
  inner_join(., taxonomy, by = "ASV") %>%
  group_by(Sample_ID) %>%
  select(Sample_ID, Microbial_Trtmt, Stressor, Trtmt_Combo, ASV, Source, count)
# This above command turned our otu_counts table into a new table
# that has archaea removed and is strictly bacteria

sampling_coverage_sw <- rarefy_asv_skinwater %>%
  group_by(Sample_ID) %>%
  summarize(n_seqs = sum(count)) 
#gives us the number of sequences in each sample
#Lowest number in an actual sample is 14020 reads

sampling_coverage_sw %>%
  ggplot(aes(x = n_seqs)) +
  geom_histogram(binwidth = 5000) +
  coord_cartesian(xlim = c(0, 90000))
# Helps visualize if there are any breaks in my data that
# will allow me to rarefy my sequences

sampling_coverage_sw %>%
  ggplot(aes(x = 1, y = n_seqs)) +
  geom_jitter()

sampling_coverage_sw %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1: nrow (.), y = n_seqs)) +
  geom_line()

sampling_coverage_sw %>%
  arrange(n_seqs) %>%
  print(n = 40)
# 2 samples (AA301S and NA301S have no sequences) will be removed
# To properly compare to water samples, will rarefy to the lowest water sequence
# That looks to be NPW5 at 1120 sequences


rarefied_asv_table_sw <- rarefy_asv_skinwater %>%
  group_by(Sample_ID) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 1120) %>%
  select(-n)  %>%
  select(-Stressor, -Microbial_Trtmt, -Trtmt_Combo) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  column_to_rownames("Sample_ID")
# Our new, rarefied (1120) ASV counts table for gut and water samples without archaea

write.table(rarefied_asv_table_sw, file = "Emerson Aim 3 - rarefied asv table skin + water.csv", sep = "," )

####### Beta Diversity - Gut and Water Samples
##Bray Curtis
#https://www.youtube.com/watch?v=G5Qckqq5Erw create PCoA plot
#https://www.youtube.com/watch?v=xyufizOpc5I nmds plot

distance_matrix_skinwater <- rarefy_asv_skinwater %>%
  group_by(Sample_ID) %>%
  select(-Microbial_Trtmt, -Stressor, -Trtmt_Combo, -Source) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 1120) %>%
  ungroup() %>%
  group_by(ASV) %>%
  select(-n) %>%
  mutate(total = sum(count)) %>%
  filter(total !=0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  as.data.frame()

rownames(distance_matrix_skinwater) <- distance_matrix_skinwater$Sample_ID
distance_matrix_skinwater <- distance_matrix_skinwater[,-1]
distance_matrix_skinwater <- as.matrix(distance_matrix_skinwater)

set.seed(19950406)
dist_sw <- avgdist(distance_matrix_skinwater, dmethod = "bray", sample = 1120)
set.seed(17)
nmds_sw <- metaMDS(dist_sw)

metadata_nmds_sw <- nmds_sw$points %>%
  as_tibble(rownames = "Sample_ID") %>%
  inner_join(., meta_skinwater, by = "Sample_ID") 

Pond_Water <- c("Natural", "Autoclaved")
Stress <- c("Vehicle Control", "Predator Cues", "CORT")

###Bray Curtis for Gut vs Water

metadata_nmds_sw %>%
  ggplot(aes(x = MDS1, y = MDS2, color = Source)) +
  geom_point(aes(shape = Source, size = 3)) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 14))  +
  theme(axis.title = element_text(face = "bold",size = 16)) +
  labs(x = "Bray-Curtis NMDS 1", y = "Bray-Curtis NMDS 2") +
  scale_color_manual(labels = c("Skin", "Water"),
                     values = c("pink", "skyblue"),
                     name = "Source") +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text( size = 12))

#### Bray Curtis for Skin Samples

rarefy_asv_skin <- inner_join(meta_skinwater, otu_counts, by = "Sample_ID") %>%
  inner_join(., taxonomy, by = "ASV") %>%
  group_by(Sample_ID) %>%
  select(Sample_ID, Microbial_Trtmt, Stressor, Trtmt_Combo, ASV, Source, count)
# This above command turned our otu_counts table into a new table
# that has archaea removed and is strictly bacteria

sampling_coverage_skin <- rarefy_asv_skin %>%
  group_by(Sample_ID) %>%
  summarize(n_seqs = sum(count)) 
#gives us the number of sequences in each sample


sampling_coverage_skin %>%
  ggplot(aes(x = n_seqs)) +
  geom_histogram(binwidth = 5000) +
  coord_cartesian(xlim = c(0, 90000))
# Helps visualize if there are any breaks in my data that
# will allow me to rarefy my sequences

sampling_coverage_skin %>%
  ggplot(aes(x = 1, y = n_seqs)) +
  geom_jitter()

sampling_coverage_skin %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1: nrow (.), y = n_seqs)) +
  geom_line()

sampling_coverage_skin %>%
  arrange(n_seqs) %>%
  print(n = 40)
# So I will rarefy to 11882 (lowest number of sequences in a sample), which will cause me to drop off the blanks 
# and I will also filter out our pond water samples

rarefied_asv_table_skin <- rarefy_asv_skin %>%
  group_by(Sample_ID) %>%
  filter(Sample_ID != 'Kyle061-NPW1') %>%
  filter(Sample_ID != 'Kyle062-NPW2') %>%
  filter(Sample_ID != 'Kyle063-NPW3') %>%
  filter(Sample_ID != 'Kyle064-NPW4') %>%
  filter(Sample_ID != 'Kyle065-NPW5') %>%
  filter(Sample_ID != 'Kyle066-NPW6') %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 11880) %>%
  select(-n)  %>%
  select(-Stressor, -Microbial_Trtmt, -Trtmt_Combo, -Source) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  column_to_rownames("Sample_ID")
# Our new, rarefied (11880) ASV counts table for gut and water samples without archaea

write.table(rarefied_asv_table_skin, file = "Emerson Aim 3 - rarefied asv table skin.csv", sep = "," )

####### Beta Diversity - Skin Samples
##Bray Curtis
#https://www.youtube.com/watch?v=G5Qckqq5Erw create PCoA plot
#https://www.youtube.com/watch?v=xyufizOpc5I nmds plot

distance_matrix_skin <- rarefy_asv_skin %>%
  group_by(Sample_ID) %>%
  filter(Sample_ID != 'Kyle061-NPW1') %>%
  filter(Sample_ID != 'Kyle062-NPW2') %>%
  filter(Sample_ID != 'Kyle063-NPW3') %>%
  filter(Sample_ID != 'Kyle064-NPW4') %>%
  filter(Sample_ID != 'Kyle065-NPW5') %>%
  filter(Sample_ID != 'Kyle066-NPW6') %>%
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

metadata_nmds_skin <- nmds_skin$points %>%
  as_tibble(rownames = "Sample_ID") %>%
  inner_join(., meta_skin, by = "Sample_ID") 

Pond_Water <- c("Natural", "Autoclaved")
Stress <- c("Vehicle Control", "Predator Cues", "CORT")

###Bray Curtis for Skin Samples 

metadata_nmds_skin %>%
  ggplot(aes(x = MDS1, y = MDS2, color = factor(Stressor, level=c('B', 'A','C')))) +
  geom_point(aes(shape = Microbial_Trtmt, size = 3)) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 14))  +
  theme(axis.title = element_text(face = "bold",size = 16)) +
  labs(x = "Bray-Curtis NMDS 1", y = "Bray-Curtis NMDS 2") +
  scale_color_manual(labels = c("Vehicle", "Predator", "CORT"),
                     values = c("gray81", "darkred", "deepskyblue1"),
                     name = "Pond Water") +
  scale_shape_manual(values = c(16,18)) +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text( size = 12))


##Plot of our Bray Curtis NMDS (rarefied)
## Bray Curtis Analysis: https://www.youtube.com/watch?v=oLf0EpMJ4yA
adonis2(dist_skin~metadata_nmds_skin$Microbial_Trtmt*metadata_nmds_skin$Stressor)
#Essentially our PERMANOVA. Significant on treatment
#                                                                Df SumOfSqs      R2
# metadata_nmds_skin$Microbial_Trtmt                              1   0.5625 0.12091
# metadata_nmds_skin$Stressor                                     2   0.4591 0.09868
# metadata_nmds_skin$Microbial_Trtmt:metadata_nmds_skin$Stressor  2   0.4491 0.09653
# Residual                                                       14   3.1816 0.68388
# Total                                                          19   4.6522 1.00000
#                                                                   F     Pr(>F)   
# metadata_nmds_skin$Microbial_Trtmt                             2.4753  0.004 **
# metadata_nmds_skin$Stressor                                    1.0100  0.464   
# metadata_nmds_skin$Microbial_Trtmt:metadata_nmds_skin$Stressor 0.9880  0.483 

bd_micro <- betadisper(dist_skin, metadata_nmds_skin$Microbial_Trtmt)
anova(bd_micro)
#Our Permdisp, nonsignificant for microbial treatment
#            Df   Sum Sq   Mean Sq   F value   Pr(>F)
# Groups     1  0.005594  0.0055944  0.6444   0.4326

bd_stressor <- betadisper(dist_skin, metadata_nmds_skin$Stressor)
anova(bd_stressor)
#Nonsignificant for stessor 
#            Df   Sum Sq   Mean Sq  F value   Pr(>F)
# Groups     2  0.003783 0.0018917   0.2168   0.8073


########## Alpha Diversity Metrics
## No. Observed ASVs

observed_asvs <- rarefied_asv_table_skin %>%
  rarefy(sample = 11880) %>%
  as_tibble(rownames = "Sample_ID") %>%
  select(Sample_ID, Obs_ASVs=value)

write.table(observed_asvs, file = "Emerson Aim 3 - no. observed skin asvs.csv", sep = "," )

## Shannon Diversity Index

rarefied_asv_table_skin %>%
  rrarefy(sample = 11880) %>%
  diversity()

#Now, we have shannon diversity values for our samples
#But, we want to repeat this to get the average
#So we are going to repeat this process, average that output for our rarefied shannon

shannon_iteration <- function(){
  
  rarefied_asv_table_skin %>%
    rrarefy(sample = 11880) %>%
    diversity()
  
}

rarefied_shannon <- replicate (100, shannon_iteration()) %>% 
  as_tibble(rownames = "Sample_ID", .name_repair = "unique" ) %>%
  pivot_longer(-Sample_ID) %>%
  group_by(Sample_ID) %>%
  summarize(shannon = mean(value))

write.table(rarefied_shannon, file = "Emerson Aim 3 - shannon skin diversity.csv", sep = "," )

## Faiths Phylogenetic Diversity
# #https://search.r-project.org/CRAN/refmans/abdiv/html/faith_pd.html

faiths.df <- rarefied_asv_table_skin

faiths.df <- as.matrix(faiths.df)
#Turned our rarefied asv table that we have been using into a matrix

Faith_PD <- pd(faiths.df, phy_tree(ps), include.root = TRUE) 
#These are our faiths phylogenetic diversity values!
#SR represents species richness values, but we will keep our other ones

write.table(Faith_PD, file = "Emerson Aim 3 - Final skin FaithsPD.csv", sep = "," )


####### Updated Code for revisions(KE - 2-19-25)
####### Beta Diversity - Skin Samples
##Jaccard - starting with the same distance matrix used for Bray Curtis

set.seed(19950406)
dist_skin_jaccard <- avgdist(distance_matrix_skin, dmethod = "jaccard", binary = TRUE, sample = 11880)
set.seed(17)
nmds_skin_jaccard <- metaMDS(dist_skin_jaccard)
##have to specify "binary = TRUE" for jaccard analysis to make sure that the species input is binary
##https://stats.stackexchange.com/questions/242110/nmds-from-jaccard-and-bray-curtis-identical-is-that-a-bad-thing

metadata_nmds_skin_jaccard <- nmds_skin_jaccard$points %>%
  as_tibble(rownames = "Sample_ID") %>%
  inner_join(., meta_skin, by = "Sample_ID") 

Pond_Water <- c("Natural", "Autoclaved")
Stress <- c("Vehicle Control", "Predator Cues", "CORT")

###Jaccard Ordination plot for skin Samples 

metadata_nmds_skin_jaccard %>%
  ggplot(aes(x = MDS1, y = MDS2, color = factor(Stressor, level=c('B', 'A','C')))) +
  geom_point(aes(shape = Microbial_Trtmt, size = 3)) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 14))  +
  theme(axis.title = element_text(face = "bold",size = 16)) +
  labs(x = "Jaccard NMDS 1", y = "Jaccard NMDS 2") +
  scale_color_manual(labels = c("Vehicle", "Predator", "CORT"),
                     values = c("gray81", "darkred", "deepskyblue1"),
                     name = "Pond Water") +
  scale_shape_manual(labels = c("Natural", "Autoclaved"),
                     values = c(16,18),
                     name = "Pond Water") +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size = 16)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  guides(size = "none")


##Plot of our Jaccard NMDS (rarefied)
adonis2(dist_skin_jaccard~metadata_nmds_skin_jaccard$Microbial_Trtmt*metadata_nmds_skin_jaccard$Stressor)
#PERMANOVA. Microbial Trtmt Significant
#                                                              
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_skin_jaccard ~ metadata_nmds_skin_jaccard$Microbial_Trtmt * metadata_nmds_skin_jaccard$Stressor)
#                                                                                Df SumOfSqs      R2      F   Pr(>F)    
# metadata_nmds_skin_jaccard$Microbial_Trtmt                                      1   0.7083 0.11496 2.3971  0.001 ***
# metadata_nmds_skin_jaccard$Stressor                                             2   0.6575 0.10672 1.1127  0.148    
# metadata_nmds_skin_jaccard$Microbial_Trtmt:metadata_nmds_skin_jaccard$Stressor  2   0.6587 0.10691 1.1147  0.152    
# Residual                                                                       14   4.1365 0.67140                  
# Total                                                                          19   6.1609 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

bd_micro_skin_jac <- betadisper(dist_skin_jaccard, metadata_nmds_skin_jaccard$Microbial_Trtmt)
anova(bd_micro_skin_jac)
#Our Permdisp, significant for microbial treatment

# #Analysis of Variance Table
# 
# Response: Distances
#           Df   Sum Sq   Mean Sq  F value    Pr(>F)   
# Groups     1 0.014856 0.0148558  14.118 0.001442 **
# Residuals 18 0.018940 0.0010522  

bd_stressor_skin_jac <- betadisper(dist_skin_jaccard, metadata_nmds_skin_jaccard$Stressor)
anova(bd_stressor_skin_jac)
#Nonsignificant for stessor

# Analysis of Variance Table
# 
# Response: Distances
#           Df    Sum Sq   Mean Sq F value Pr(>F)
# Groups     2 0.0028984 0.0014492  1.0357 0.3763
# Residuals 17 0.0237875 0.0013993 