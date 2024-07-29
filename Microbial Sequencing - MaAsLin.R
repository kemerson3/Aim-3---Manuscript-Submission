##Microbial Sequencing - MaAsLin
library(lme4)
library(car)
library(readr)
library(moments)
library(psych)
library(pastecs)
library(ggplot2)
library(Maaslin2)

##import metadata and phyla/genus abundance data

file.choose()
Maslin_meta <- read.table("C:\\R\\Aim 3 - Stress and the Microbiota\\Emerson Aim 3 - Microbial Map - Gut & Water.txt", sep = '\t',
                          header=T) %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle061-NPW1') %>%
  filter(Sample_ID != 'Kyle062-NPW2') %>%
  filter(Sample_ID != 'Kyle063-NPW3') %>%
  filter(Sample_ID != 'Kyle064-NPW4') %>%
  filter(Sample_ID != 'Kyle065-NPW5') %>%
  filter(Sample_ID != 'Kyle066-NPW6') 


Maslin_meta$Microbial_Trtmt = factor(Maslin_meta$Microbial_Trtmt)
Maslin_meta$Stressor = factor(Maslin_meta$Stressor)
Maslin_meta$Trtmt_Combo = factor(Maslin_meta$Trtmt_Combo)


file.choose()
Genera <- read.table("C:\\R\\Aim 3 - Stress and the Microbiota\\Microbiome - Gut and Water Excel Files\\Emerson Aim 3 - gut & water genus abundance.txt", sep = '\t',
                     header=T) %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle061-NPW1') %>%
  filter(Sample_ID != 'Kyle062-NPW2') %>%
  filter(Sample_ID != 'Kyle063-NPW3') %>%
  filter(Sample_ID != 'Kyle064-NPW4') %>%
  filter(Sample_ID != 'Kyle065-NPW5') %>%
  filter(Sample_ID != 'Kyle066-NPW6') 

Phyla <- read.table("C:\\R\\Aim 3 - Stress and the Microbiota\\Microbiome - Gut and Water Excel Files\\Emerson Aim 3 - gut & water phyla abundance.txt", sep = '\t',
                    header=T) %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle061-NPW1') %>%
  filter(Sample_ID != 'Kyle062-NPW2') %>%
  filter(Sample_ID != 'Kyle063-NPW3') %>%
  filter(Sample_ID != 'Kyle064-NPW4') %>%
  filter(Sample_ID != 'Kyle065-NPW5') %>%
  filter(Sample_ID != 'Kyle066-NPW6') 

##create tmp files to work around bug: "https://github.com/biobakery/Maaslin2/issues/1"
library(tidyverse)
tmp_data_phyla = tempfile(pattern = "data")
write_delim(Phyla, tmp_data_phyla, delim = "\t")

tmp_data_gen = tempfile(pattern = "data")
write_delim(Genera, tmp_data_gen, delim = "\t")

tmp_metadata = tempfile(pattern = "metadata")
write_delim(Maslin_meta, tmp_metadata, delim = "\t")

##Next step will not work if the Sample_IDs do not match
##Need to go into excel files for phyla and genera, make sure Sample_IDs
##Match in all 3 files. Save new txt files if need be.

##Run at phyla level - pond
Maaslin_phyla_pond <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_pond_output', 
  transform = "AST",
  fixed_effects = c('Microbial_Trtmt'),
  standardize = FALSE)

##Run at genus level - pond
Maaslin_gen_pond <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_pond_output', 
  transform = "AST",
  fixed_effects = c('Microbial_Trtmt'),
  standardize = FALSE)

##Run at phyla level - stressors
Maaslin_phyla_stressor <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_stressor_output', 
  transform = "AST",
  fixed_effects = c('Stressor'),
  reference = c('Stressor','B'),
  standardize = FALSE)

##Run at genus level - stressor
Maaslin_gen_stressor <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_stressor_output', 
  transform = "AST",
  fixed_effects = c('Stressor'),
  reference = c('Stressor','B'),
  standardize = FALSE)

##Run at phyla level - both
Maaslin_phyla_trtmtcmbo <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_trtmt_output', 
  transform = "AST",
  fixed_effects = c('Stressor', 'Microbial_Trtmt'),
  reference = c('Stressor','B'),
  standardize = FALSE)

##Run at genus level - both
Maaslin_gen_trtmtcmbo <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_trtmt_output', 
  transform = "AST",
  fixed_effects = c('Stressor', 'Microbial_Trtmt'),
  reference = c('Stressor','B'),
  standardize = FALSE)

##### Now, we want to find associations between specific bacterial taxa and some aspects of host physiology
file.choose()
Maslin_meta <- read.table("C:\\R\\Aim 3 - Stress and the Microbiota\\Emerson Aim 3 - Microbial Map - Gut & Water.txt", sep = '\t',
                          header=T) %>%
  filter(Sample_ID != 'Kyle073-KE-Blank') %>%
  filter(Sample_ID != 'Kyle074-PCR1Blank1') %>%
  filter(Sample_ID != 'Kyle075-PCR1Blank2') %>%
  filter(Sample_ID != 'Kyle061-NPW1') %>%
  filter(Sample_ID != 'Kyle062-NPW2') %>%
  filter(Sample_ID != 'Kyle063-NPW3') %>%
  filter(Sample_ID != 'Kyle064-NPW4') %>%
  filter(Sample_ID != 'Kyle065-NPW5') %>%
  filter(Sample_ID != 'Kyle066-NPW6') 


Maslin_meta$Microbial_Trtmt = factor(Maslin_meta$Microbial_Trtmt)
Maslin_meta$Stressor = factor(Maslin_meta$Stressor)
Maslin_meta$Trtmt_Combo = factor(Maslin_meta$Trtmt_Combo)
Maslin_meta$Gosner_Stage = as.numeric(Maslin_meta$Gosner_Stage)

##Relative Brain Mass
Maaslin_phyla_brainmass <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_brainmass_output', 
  transform = "AST",
  fixed_effects = c('Relative_Brain_Mass'),
  standardize = FALSE)

Maaslin_genus_brainmass <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_brainmass_output', 
  transform = "AST",
  fixed_effects = c('Relative_Brain_Mass'),
  standardize = FALSE)

##Relative Brain Shape - PC1
Maaslin_phyla_brainshape1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_brainshapePC1_output', 
  transform = "AST",
  fixed_effects = c('BrainShape_RC1'),
  standardize = FALSE)

Maaslin_genus_brainshape1 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_brainshapePC1_output', 
  transform = "AST",
  fixed_effects = c('BrainShape_RC1'),
  standardize = FALSE)

##Relative Brain Shape - PC2
Maaslin_phyla_brainshape2 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_brainshapePC2_output', 
  transform = "AST",
  fixed_effects = c('BrainShape_RC2'),
  standardize = FALSE)

Maaslin_genus_brainshape2 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_brainshapePC2_output', 
  transform = "AST",
  fixed_effects = c('BrainShape_RC2'),
  standardize = FALSE)

##Relative Brain Shape - PC3
Maaslin_phyla_brainshape3 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_brainshapePC3_output', 
  transform = "AST",
  fixed_effects = c('BrainShape_RC3'),
  standardize = FALSE)

Maaslin_genus_brainshape3 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_brainshapePC3_output', 
  transform = "AST",
  fixed_effects = c('BrainShape_RC3'),
  standardize = FALSE)

##Baseline Behavior - PC1
Maaslin_phyla_baseline1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_baseline1_output', 
  transform = "AST",
  fixed_effects = c('Baseline_RC1'),
  standardize = FALSE)

Maaslin_genus_baseline1 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_baseline1_output', 
  transform = "AST",
  fixed_effects = c('Baseline_RC1'),
  standardize = FALSE)

##Baseline Behavior - PC2
Maaslin_phyla_baseline2 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_baseline2_output', 
  transform = "AST",
  fixed_effects = c('Baseline_RC2'),
  standardize = FALSE)

Maaslin_genus_baseline2 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_baseline2_output', 
  transform = "AST",
  fixed_effects = c('Baseline_RC2'),
  standardize = FALSE)

##Visual Empty Behavior - PC1
Maaslin_phyla_empty1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_empty1_output', 
  transform = "AST",
  fixed_effects = c('Empty_RC1'),
  standardize = FALSE)

Maaslin_genus_empty1 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_empty1_output', 
  transform = "AST",
  fixed_effects = c('Empty_RC1'),
  standardize = FALSE)

##Visual Empty Behavior - PC2
Maaslin_phyla_empty2 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_empty2_output', 
  transform = "AST",
  fixed_effects = c('Empty_RC2'),
  standardize = FALSE)

Maaslin_genus_empty2 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_empty2_output', 
  transform = "AST",
  fixed_effects = c('Empty_RC2'),
  standardize = FALSE)

##Visual Empty Behavior - PC3
Maaslin_phyla_empty3 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_empty3_output', 
  transform = "AST",
  fixed_effects = c('Empty_RC3'),
  standardize = FALSE)

Maaslin_genus_empty3 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_empty3_output', 
  transform = "AST",
  fixed_effects = c('Empty_RC3'),
  standardize = FALSE)

##Visual Food Behavior - PC1
Maaslin_phyla_food1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_food1_output', 
  transform = "AST",
  fixed_effects = c('Food_RC1'),
  standardize = FALSE)

Maaslin_genus_food1 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_food1_output', 
  transform = "AST",
  fixed_effects = c('Food_RC1'),
  standardize = FALSE)

##Visual Food Behavior - PC2
Maaslin_phyla_food2 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_food2_output', 
  transform = "AST",
  fixed_effects = c('Food_RC2'),
  standardize = FALSE)

Maaslin_genus_food2 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_food2_output', 
  transform = "AST",
  fixed_effects = c('Food_RC2'),
  standardize = FALSE)

##Visual Food Behavior - PC3
Maaslin_phyla_food3 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_food3_output', 
  transform = "AST",
  fixed_effects = c('Food_RC3'),
  standardize = FALSE)

Maaslin_genus_food3 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_food3_output', 
  transform = "AST",
  fixed_effects = c('Food_RC3'),
  standardize = FALSE)

##Olfactory Behavior - PC1
Maaslin_phyla_olfactory1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_olfactory1_output', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC1'),
  standardize = FALSE)

Maaslin_genus_olfactory1 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_olfactory1_output', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC1'),
  standardize = FALSE)

##Olfactory Behavior - PC2
Maaslin_phyla_olfactory2 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_olfactory2_output', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC2'),
  standardize = FALSE)

Maaslin_genus_olfactory2 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_olfactory2_output', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC2'),
  standardize = FALSE)

##Body Mass
Maaslin_phyla_bodymass <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_bodymass_output', 
  transform = "AST",
  fixed_effects = c('Body_Mass'),
  standardize = FALSE)

Maaslin_genus_bodymass <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_bodymass_output', 
  transform = "AST",
  fixed_effects = c('Body_Mass'),
  standardize = FALSE)

##Gosner Stage
Maaslin_phyla_gosner <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_phyla_gosner_output', 
  transform = "AST",
  fixed_effects = c('Gosner_Stage'),
  standardize = FALSE)

Maaslin_genus_gosner <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Aim 3 - Stress and the Microbiota/EM3_Maaslin_genus_gosner_output', 
  transform = "AST",
  fixed_effects = c('Gosner_Stage'),
  standardize = FALSE)
