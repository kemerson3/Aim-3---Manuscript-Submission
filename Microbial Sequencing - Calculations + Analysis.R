##Microbial Sequencing - Calculations and Analysis - FULL ANALYSIS
library(devtools)
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
library(lme4)
library(car)
library(emmeans)

##### Alpha Diversity

file.choose()
df_alpha <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Emerson Aim 3 - alpha diversity full.csv")
df_alpha$Microbial_Trtmt = factor(df_alpha$Microbial_Trtmt)
df_alpha$Stressor = factor(df_alpha$Stressor)
df_alpha$Trtmt_Combo = factor(df_alpha$Trtmt_Combo)
df_alpha$Source = factor(df_alpha$Source)

##################### Alpha Diversity - Source
### ASV
asv_source_glmm <- glm(Obs_ASVs~Source, data = df_alpha, family = "gaussian")
Anova(asv_source_glmm)

# Response: Obs_ASVs
# LR        Chisq Df Pr(>Chisq)    
# Source    22.35  2  1.402e-05 ***

pairwise.t.test(df_alpha$Obs_ASVs, df_alpha$Source, p.adj = 'bonferroni')
#          Gut    Skin   
# Skin  0.1521    -      
# Water 0.0033**  6.6e-05***
# P value adjustment method: bonferroni
# NOTE: Going to use Tukey from now on

emmeans(asv_source_glmm, list(pairwise ~ Source), adjust = "tukey")
# 1            estimate   SE df t.ratio p.value
# Gut - Skin         47 23.5 47   2.005  0.1221
# Gut - Water      -123 35.4 47  -3.477  0.0031
# Skin - Water     -170 36.1 47  -4.714  0.0001
# P value adjustment: tukey method for comparing a family of 3 estimates 

ggplot(df_alpha, aes(x= Source, y = Obs_ASVs, fill = Source)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  geom_jitter(width = .2, size = 2.5, shape = 21, color = "black") +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("pink", "green", "skyblue"),
                    name = "Source",
                    labels = c("Gut", "Skin", "Water")) +
  labs(x = "Source", y = "No. Observed ASVs") +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.12, .92)) +
  theme(legend.title = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "none")

### Shannon

shannon_source_glmm <- glm(shannon~Source, data = df_alpha, family = "gaussian")
Anova(shannon_source_glmm)

# Response: Obs_ASVs
# LR        Chisq Df Pr(>Chisq)    
# Source   39.841  2  2.232e-09 ***

emmeans(shannon_source_glmm, list(pairwise ~ Source), adjust = "tukey")
# $`pairwise differences of Source`
# 1            estimate    SE df t.ratio p.value
# Gut - Skin     0.0533 0.161 47   0.330  0.9417
# Gut - Water   -1.4380 0.243 47  -5.910  <.0001
# Skin - Water  -1.4913 0.248 47  -6.010  <.0001
# 
# P value adjustment: tukey method for comparing a family of 3 estimates 

ggplot(df_alpha, aes(x= Source, y = shannon, fill = Source)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  geom_jitter(width = .2, size = 2.5, shape = 21, color = "black") +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("pink", "green", "skyblue"),
                    name = "Source",
                    labels = c("Gut", "Skin", "Water")) +
  labs(x = "Source", y = "Shannon Diversity Index") +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.12, .92)) +
  theme(legend.title = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "none")

### Faith

faith_source_glmm <- glm(PD~Source, data = df_alpha, family = "gaussian")
Anova(faith_source_glmm)

# Response: Obs_ASVs
# LR        Chisq Df Pr(>Chisq)    
# Source   49.897  2  1.462e-11 ***

emmeans(faith_source_glmm, list(pairwise ~ Source), adjust = "tukey")
# $`pairwise differences of Source`
# 1            estimate    SE df t.ratio p.value
# Gut - Skin       6.75 5.77 47   1.171  0.4759
# Gut - Water    -54.66 8.69 47  -6.289  <.0001
# Skin - Water   -61.41 8.86 47  -6.928  <.0001
# P value adjustment: tukey method for comparing a family of 3 estimates 

ggplot(df_alpha, aes(x= Source, y = PD, fill = Source)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  geom_jitter(width = .2, size = 2.5, shape = 21, color = "black") +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("pink", "green", "skyblue"),
                    name = "Source",
                    labels = c("Gut", "Skin", "Water")) +
  labs(x = "Source", y = "Faith's Phylogenetic Diversity") +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.12, .92)) +
  theme(legend.title = element_text(face = "bold", size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "none")
 
########## Alpha Diversity - Gut Samples
### ASV
df_alpha_gut <- df_alpha %>%
  filter(Source == "Gut")
  
asv_gut_glmm <- glm(Obs_ASVs~Microbial_Trtmt*Stressor, data = df_alpha_gut, family = "gaussian")
Anova(asv_gut_glmm)

# Microbial_Trtmt            41.109  1   1.44e-10 ***
# Stressor                    0.699  2     0.7050    
# Microbial_Trtmt:Stressor    3.693  2     0.1578

emmeans(asv_gut_glmm, list(pairwise ~ Stressor), adjust = "tukey")
# $`pairwise differences of Stressor`
# 1     estimate   SE df t.ratio p.value
# A - B   15.441 20.8 18   0.743  0.7417
# A - C    0.811 20.8 18   0.039  0.9992
# B - C  -14.630 20.8 18  -0.704  0.7643
# 
# Results are averaged over the levels of: Microbial_Trtmt 
# P value adjustment: tukey method for comparing a family of 3 estimates

pairwise.t.test(df_alpha_gut$Obs_ASVs, df_alpha_gut$Stressor, p.adj = 'bonferroni')
#     A B
# B   1 -
# C   1 1
# P value adjustment method: bonferroni

Stress <- c("Vehicle Control", "Predator Cues", "CORT")
#For labeling X axis

ggplot(df_alpha_gut, aes(x= factor(Stressor, level=c('B', 'A','C')), y = Obs_ASVs, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Pond Water",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "No. Observed ASVs") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16, colour = "black"))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.83, .92)) +
  theme(legend.title = element_text(face = "bold", size = 18)) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.position= c(.16, .9))

###Shannon

shannon_gut_glmm <- glm(shannon~Microbial_Trtmt*Stressor, data = df_alpha_gut, family = "gaussian")
Anova(shannon_gut_glmm)

#                         LR Chisq Df Pr(>Chisq)    
# Microbial_Trtmt            31.888  1  1.633e-08 ***
# Stressor                    4.394  2     0.1111    
# Microbial_Trtmt:Stressor    3.996  2     0.1356 

emmeans(shannon_gut_glmm, list(pairwise ~ Stressor), adjust = "tukey")
# $`pairwise differences of Stressor`
# 1     estimate   SE df t.ratio p.value
# A - B   0.0644 0.14 18   0.460  0.8906
# A - C  -0.2158 0.14 18  -1.541  0.2961
# B - C  -0.2803 0.14 18  -2.001  0.1406
# 
# Results are averaged over the levels of: Microbial_Trtmt 
# P value adjustment: tukey method for comparing a family of 3 estimates 


ggplot(df_alpha_gut, aes(x= factor(Stressor, level=c('B', 'A','C')), y = shannon, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "Shannon Diversity Index") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16, colour = "black"))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.83, .92)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position="none")

### Faith's

faith_gut_glmm <- glm(PD~Microbial_Trtmt*Stressor, data = df_alpha_gut, family = "gaussian")
Anova(faith_gut_glmm)

#                         LR Chisq Df Pr(>Chisq)    
# Microbial_Trtmt           16.9500  1  3.838e-05 ***
# Stressor                   1.2882  2     0.5251    
# Microbial_Trtmt:Stressor   3.0101  2     0.2220 

emmeans(faith_gut_glmm, list(pairwise ~ Stressor), adjust = "tukey")
# $`pairwise differences of Stressor`
# 1     estimate   SE df t.ratio p.value
# A - B     4.91 4.55 18   1.079  0.5385
# A - C     1.07 4.55 18   0.236  0.9699
# B - C    -3.84 4.55 18  -0.844  0.6815
# 
# Results are averaged over the levels of: Microbial_Trtmt 
# P value adjustment: tukey method for comparing a family of 3 estimates 


ggplot(df_alpha_gut, aes(x= factor(Stressor, level=c('B', 'A','C')), y = PD, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "Faith's Phylogenetic Diversity") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16, colour = "black"))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.83, .92)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position="none")

########## Alpha Diversity - Skin Samples
### ASV
df_alpha_skin <- df_alpha %>%
  filter(Source == "Skin")

asv_skin_glmm <- glm(Obs_ASVs~Microbial_Trtmt*Stressor, data = df_alpha_skin, family = "gaussian")
Anova(asv_skin_glmm)

#         LR                   Chisq Df Pr(>Chisq)    
# Microbial_Trtmt           21.3478  1  3.831e-06 ***
# Stressor                   5.3882  2     0.0676 .  
# Microbial_Trtmt:Stressor   1.1700  2     0.5571 

emmeans(asv_skin_glmm, list(pairwise ~ Stressor), adjust = "tukey")
# $`pairwise differences of Stressor`
# 1     estimate   SE df t.ratio p.value
# A - B   -27.60 13.5 14  -2.037  0.1399
# A - C   -26.57 14.0 14  -1.903  0.1745
# B - C     1.03 12.2 14   0.084  0.9961
# 
# Results are averaged over the levels of: Microbial_Trtmt 
# P value adjustment: tukey method for comparing a family of 3 estimates


ggplot(df_alpha_skin, aes(x= factor(Stressor, level=c('B', 'A','C')), y = Obs_ASVs, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Pond Water",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "No. Observed ASVs") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.83, .92)) +
  theme(legend.title = element_text(face = "bold", size = 18)) +
  theme(legend.text = element_text(size = 16)) +
  theme(legend.position= "none")

###Shannon

shannon_skin_glmm <- glm(shannon~Microbial_Trtmt*Stressor, data = df_alpha_skin, family = "gaussian")
Anova(shannon_skin_glmm)

#                         LR Chisq Df Pr(>Chisq)    
# Microbial_Trtmt           1.53228  1     0.2158
# Stressor                  1.46444  2     0.4808
# Microbial_Trtmt:Stressor  0.71329  2     0.7000

emmeans(shannon_skin_glmm, list(pairwise ~ Stressor), adjust = "tukey")
# $`pairwise differences of Stressor`
# 1     estimate   SE df t.ratio p.value
# A - B   -0.217 0.270 14  -0.803  0.7074
# A - C   -0.328 0.279 14  -1.175  0.4863
# B - C   -0.111 0.244 14  -0.453  0.8937
# 
# Results are averaged over the levels of: Microbial_Trtmt 
# P value adjustment: tukey method for comparing a family of 3 estimates 


ggplot(df_alpha_skin, aes(x= factor(Stressor, level=c('B', 'A','C')), y = shannon, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "Shannon Diversity Index") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.83, .92)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position="none")

### Faith's

faith_skin_glmm <- glm(PD~Microbial_Trtmt*Stressor, data = df_alpha_skin, family = "gaussian")
Anova(faith_skin_glmm)

#                         LR Chisq Df Pr(>Chisq)    
# Microbial_Trtmt           13.6122  1  0.0002247 ***
# Stressor                   4.2135  2  0.1216343    
# Microbial_Trtmt:Stressor   1.2318  2  0.5401481 

emmeans(faith_skin_glmm, list(pairwise ~ Stressor), adjust = "tukey")
# $`pairwise differences of Stressor`
# 1     estimate   SE df t.ratio p.value
# A - B   -5.198 3.26 14  -1.595  0.2799
# A - C   -5.625 3.36 14  -1.675  0.2489
# B - C   -0.427 2.94 14  -0.145  0.9884
# 
# Results are averaged over the levels of: Microbial_Trtmt 
# P value adjustment: tukey method for comparing a family of 3 estimates 


ggplot(df_alpha_skin, aes(x= factor(Stressor, level=c('B', 'A','C')), y = PD, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "Faith's Phylogenetic Diversity") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.83, .92)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position="none")


####Gut Microbiome - Physiology Associations
file.choose()
df <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota\\Aim 3 (Restarted)\\Aim 3 metadata.csv")
df$Water_Trtmt = factor(df$Water_Trtmt)
df$Stressor = factor(df$Stressor)
df$Trtmt_Combo = factor(df$Trtmt_Combo)
df$Replicate = factor(df$Replicate)
df$Gosner_Stage = as.numeric(df$Gosner_Stage)

####Body Mass ~ Alpha Diversity
asv_bodymass_glmm <- glm(Body_Mass~Water_Trtmt*Stressor*Gut_ASV, data = df, family = "gaussian")
Anova(asv_bodymass_glmm)
# Response: Body_Mass
#                               LR Chisq Df Pr(>Chisq)
# Water_Trtmt                    1.3858  1     0.2391
# Stressor                       2.3144  2     0.3144
# Gut_ASV                        0.4003  1     0.5269
# Water_Trtmt:Stressor           0.3528  2     0.8383
# Water_Trtmt:Gut_ASV            0.0257  1     0.8727
# Stressor:Gut_ASV               0.4922  2     0.7819
# Water_Trtmt:Stressor:Gut_ASV   3.6051  2     0.1649

shannon_bodymass_glmm <- glm(Body_Mass~Water_Trtmt*Stressor*Gut_Shannon, data = df, family = "gaussian")
Anova(shannon_bodymass_glmm)
# Response: Body_Mass
#                                   LR Chisq Df Pr(>Chisq)
# Water_Trtmt                       0.06454  1     0.7995
# Stressor                          2.38556  2     0.3034
# Gut_Shannon                       0.54912  1     0.4587
# Water_Trtmt:Stressor              0.25131  2     0.8819
# Water_Trtmt:Gut_Shannon           0.10744  1     0.7431
# Stressor:Gut_Shannon              0.10291  2     0.9498
# Water_Trtmt:Stressor:Gut_Shannon  1.81555  2     0.4034

faith_bodymass_glmm <- glm(Body_Mass~Water_Trtmt*Stressor*Gut_PD, data = df, family = "gaussian")
Anova(faith_bodymass_glmm)
# Response: Body_Mass
#                              LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                   4.1035  1    0.04279 *
# Stressor                      3.2641  2    0.19553  
# Gut_PD                        1.3071  1    0.25292  
# Water_Trtmt:Stressor          1.8115  2    0.40425  
# Water_Trtmt:Gut_PD            0.0595  1    0.80736  
# Stressor:Gut_PD               2.4417  2    0.29498  
# Water_Trtmt:Stressor:Gut_PD   4.0308  2    0.13327  


####Gosner Stage ~ Alpha Diversity
asv_gosner_glmm <- glm(Gosner_Stage~Water_Trtmt*Stressor*Gut_ASV, data = df, family = "gaussian")
Anova(asv_gosner_glmm)
# Response: Gosner_Stage
#                               LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                    1.7970  1    0.18008  
# Stressor                       3.6848  2    0.15843  
# Gut_ASV                        0.8494  1    0.35673  
# Water_Trtmt:Stressor           7.6130  2    0.02223 *
# Water_Trtmt:Gut_ASV            0.2879  1    0.59159  
# Stressor:Gut_ASV               5.0396  2    0.08048 .
# Water_Trtmt:Stressor:Gut_ASV   0.4070  2    0.81587  

shannon_gosner_glmm <- glm(Gosner_Stage~Water_Trtmt*Stressor*Gut_Shannon, data = df, family = "gaussian")
Anova(shannon_gosner_glmm)
# Response: Gosner_Stage
#                                   LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                        0.3558  1    0.55087  
# Stressor                           4.7281  2    0.09404 .
# Gut_Shannon                        1.1385  1    0.28596  
# Water_Trtmt:Stressor               4.3419  2    0.11407  
# Water_Trtmt:Gut_Shannon            0.0370  1    0.84754  
# Stressor:Gut_Shannon               2.2944  2    0.31753  
# Water_Trtmt:Stressor:Gut_Shannon   1.4699  2    0.47952 

faith_gosner_glmm <- glm(Gosner_Stage~Water_Trtmt*Stressor*Gut_PD, data = df, family = "gaussian")
Anova(faith_gosner_glmm)
# Response: Gosner_Stage
#                             LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                   3.0480  1    0.08083 .
# Stressor                      4.0393  2    0.13270  
# Gut_PD                        0.9097  1    0.34019  
# Water_Trtmt:Stressor          7.0427  2    0.02956 *
# Water_Trtmt:Gut_PD            0.7097  1    0.39956  
# Stressor:Gut_PD               4.1871  2    0.12325  
# Water_Trtmt:Stressor:Gut_PD   0.9619  2    0.61820  

####Relative Brain Mass ~ Alpha Diversity
asv_brainmass_glmm <- glm(MAV_Brain_Mass~Water_Trtmt*Stressor*Gut_ASV, data = df, family = "gaussian")
Anova(asv_brainmass_glmm)
# Response: MAV_Brain_Mass
#                              LR Chisq Df Pr(>Chisq)   
# Water_Trtmt                    7.2674  1   0.007022 **
# Stressor                       5.0216  2   0.081203 . 
# Gut_ASV                        3.5942  1   0.057980 . 
# Water_Trtmt:Stressor           8.4257  2   0.014804 * 
# Water_Trtmt:Gut_ASV            0.9026  1   0.342082   
# Stressor:Gut_ASV               9.8364  2   0.007312 **
# Water_Trtmt:Stressor:Gut_ASV   0.5931  2   0.743367   

shannon_brainmass_glmm <- glm(MAV_Brain_Mass~Water_Trtmt*Stressor*Gut_Shannon, data = df, family = "gaussian")
Anova(shannon_brainmass_glmm)
# Response: MAV_Brain_Mass
#                                  LR Chisq Df Pr(>Chisq)
# Water_Trtmt                       0.23648  1     0.6268
# Stressor                          3.05677  2     0.2169
# Gut_Shannon                       0.01051  1     0.9183
# Water_Trtmt:Stressor              0.30899  2     0.8568
# Water_Trtmt:Gut_Shannon           0.07394  1     0.7857
# Stressor:Gut_Shannon              0.03609  2     0.9821
# Water_Trtmt:Stressor:Gut_Shannon  0.93224  2     0.6274

faith_brainmass_glmm <- glm(MAV_Brain_Mass~Water_Trtmt*Stressor*Gut_PD, data = df, family = "gaussian")
Anova(faith_brainmass_glmm)
# Response: MAV_Brain_Mass
#                             LR Chisq Df Pr(>Chisq)   
# Water_Trtmt                   7.5099  1   0.006136 **
# Stressor                      6.5484  2   0.037847 * 
# Gut_PD                        2.4340  1   0.118731   
# Water_Trtmt:Stressor          8.7626  2   0.012509 * 
# Water_Trtmt:Gut_PD            2.4019  1   0.121189   
# Stressor:Gut_PD              12.1041  2   0.002353 **
# Water_Trtmt:Stressor:Gut_PD   0.7016  2   0.704109   

##Alpha Diversity - Brain Shape (PC1)
asv_brainshape1_glmm <- glm(BrainShape_RC1~Water_Trtmt*Stressor*Gut_ASV, data = df, family = "gaussian")
Anova(asv_brainshape1_glmm)
# Response: BrainShape_RC1
#                              LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                    0.5486  1    0.45891  
# Stressor                       0.6603  2    0.71882  
# Gut_ASV                        0.1055  1    0.74528  
# Water_Trtmt:Stressor           6.1694  2    0.04574 *
# Water_Trtmt:Gut_ASV            3.2528  1    0.07130 .
# Stressor:Gut_ASV               6.9960  2    0.03026 *
# Water_Trtmt:Stressor:Gut_ASV   0.1688  2    0.91908  

shannon_brainshape1_glmm <- glm(BrainShape_RC1~Water_Trtmt*Stressor*Gut_Shannon, data = df, family = "gaussian")
Anova(shannon_brainshape1_glmm)
# Response: BrainShape_RC1
#                                  LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                        0.0846  1    0.77116  
# Stressor                           0.6769  2    0.71287  
# Gut_Shannon                        0.2717  1    0.60222  
# Water_Trtmt:Stressor               3.2336  2    0.19853  
# Water_Trtmt:Gut_Shannon            2.0473  1    0.15247  
# Stressor:Gut_Shannon               2.4470  2    0.29420  
# Water_Trtmt:Stressor:Gut_Shannon   5.6272  2    0.05999 .

faith_brainshape1_glmm <- glm(BrainShape_RC1~Water_Trtmt*Stressor*Gut_PD, data = df, family = "gaussian")
Anova(faith_brainshape1_glmm)
# Response: BrainShape_RC1
# LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                   0.5898  1    0.44249  
# Stressor                      0.7011  2    0.70430  
# Gut_PD                        0.0115  1    0.91468  
# Water_Trtmt:Stressor          2.6281  2    0.26873  
# Water_Trtmt:Gut_PD            1.7815  1    0.18197  
# Stressor:Gut_PD               5.2250  2    0.07335 .
# Water_Trtmt:Stressor:Gut_PD   0.1372  2    0.93371  

##Alpha Diversity - Brain Shape (PC2)
asv_brainshape2_glmm <- glm(BrainShape_RC2~Water_Trtmt*Stressor*Gut_ASV, data = df, family = "gaussian")
Anova(asv_brainshape2_glmm)
# Response: BrainShape_RC2
#                              LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                    0.3292  1    0.56613  
# Stressor                       5.0684  2    0.07932 .
# Gut_ASV                        0.0008  1    0.97800  
# Water_Trtmt:Stressor           0.2862  2    0.86665  
# Water_Trtmt:Gut_ASV            0.2183  1    0.64035  
# Stressor:Gut_ASV               1.0960  2    0.57810  
# Water_Trtmt:Stressor:Gut_ASV   0.1047  2    0.94900 

shannon_brainshape2_glmm <- glm(BrainShape_RC2~Water_Trtmt*Stressor*Gut_Shannon, data = df, family = "gaussian")
Anova(shannon_brainshape2_glmm)
# Response: BrainShape_RC2
#                                   LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                        0.7220  1    0.39548  
# Stressor                           5.6619  2    0.05896 .
# Gut_Shannon                        0.1672  1    0.68261  
# Water_Trtmt:Stressor               0.6535  2    0.72127  
# Water_Trtmt:Gut_Shannon            0.2358  1    0.62727  
# Stressor:Gut_Shannon               0.1372  2    0.93368  
# Water_Trtmt:Stressor:Gut_Shannon   0.5555  2    0.75750  

faith_brainshape2_glmm <- glm(BrainShape_RC2~Water_Trtmt*Stressor*Gut_PD, data = df, family = "gaussian")
Anova(faith_brainshape2_glmm)
# Response: BrainShape_RC2
#                              LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                   0.0241  1    0.87665  
# Stressor                      5.8988  2    0.05237 .
# Gut_PD                        0.1025  1    0.74889  
# Water_Trtmt:Stressor          0.1357  2    0.93438  
# Water_Trtmt:Gut_PD            0.1058  1    0.74502  
# Stressor:Gut_PD               1.3127  2    0.51875  
# Water_Trtmt:Stressor:Gut_PD   0.3326  2    0.84681 

##Alpha Diversity - Brain Shape (PC3)
asv_brainshape3_glmm <- glm(BrainShape_RC3~Water_Trtmt*Stressor*Gut_ASV, data = df, family = "gaussian")
Anova(asv_brainshape3_glmm)
# Response: BrainShape_RC3
#                              LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                    0.2006  1    0.65423  
# Stressor                       1.3414  2    0.51136  
# Gut_ASV                        0.0644  1    0.79974  
# Water_Trtmt:Stressor           6.4560  2    0.03964 *
# Water_Trtmt:Gut_ASV            2.7536  1    0.09703 .
# Stressor:Gut_ASV               5.1987  2    0.07432 .
# Water_Trtmt:Stressor:Gut_ASV   0.0571  2    0.97185  

shannon_brainshape3_glmm <- glm(BrainShape_RC3~Water_Trtmt*Stressor*Gut_Shannon, data = df, family = "gaussian")
Anova(shannon_brainshape3_glmm)
# Response: BrainShape_RC3
#                                  LR Chisq Df Pr(>Chisq)
# Water_Trtmt                        0.0564  1     0.8122
# Stressor                           1.7955  2     0.4075
# Gut_Shannon                        0.0101  1     0.9200
# Water_Trtmt:Stressor               1.7571  2     0.4154
# Water_Trtmt:Gut_Shannon            0.0001  1     0.9935
# Stressor:Gut_Shannon               0.5668  2     0.7532
# Water_Trtmt:Stressor:Gut_Shannon   3.7798  2     0.1511

faith_brainshape3_glmm <- glm(BrainShape_RC3~Water_Trtmt*Stressor*Gut_PD, data = df, family = "gaussian")
Anova(faith_brainshape3_glmm)
# Response: BrainShape_RC3
#                             LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                   0.0160  1    0.89934  
# Stressor                      1.4955  2    0.47343  
# Gut_PD                        0.0073  1    0.93191  
# Water_Trtmt:Stressor          4.9688  2    0.08338 .
# Water_Trtmt:Gut_PD            2.8468  1    0.09156 .
# Stressor:Gut_PD               4.9439  2    0.08442 .
# Water_Trtmt:Stressor:Gut_PD   0.1663  2    0.92021  

####Alpha Diversity - Baseline Behavior
file.choose()
df_baseline <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota Behavior\\Behavior Analysis\\Combined data sheets\\Data Sheets with Factor Scores - No Change from Baseline\\Aim3baselinefactorscores.csv")
df_baseline$Micro_Trtmt = factor(df_baseline$Micro_Trtmt)
df_baseline$Stressor = factor(df_baseline$Stressor)

##PC1
asv_baseline_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_ASV, data = df_baseline, family = "gaussian")
Anova(asv_baseline_glmm)
# Response: RC1
#                              LR Chisq Df Pr(>Chisq)   
# Micro_Trtmt                    3.7136  1   0.053972 . 
# Stressor                       0.2400  2   0.886912   
# Gut_ASV                        8.7689  1   0.003064 **
# Micro_Trtmt:Stressor           3.5280  2   0.171358   
# Micro_Trtmt:Gut_ASV            6.7432  1   0.009410 **
# Stressor:Gut_ASV               2.1519  2   0.340970   
# Micro_Trtmt:Stressor:Gut_ASV  10.7376  2   0.004660 **

shannon_baseline_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_Shannon, data = df_baseline, family = "gaussian")
Anova(shannon_baseline_glmm)
# Response: RC1
#                                  LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                       0.38977  1    0.53242  
# Stressor                          0.39332  2    0.82147  
# Gut_Shannon                       2.75485  1    0.09696 .
# Micro_Trtmt:Stressor              0.73528  2    0.69237  
# Micro_Trtmt:Gut_Shannon           1.69626  1    0.19278  
# Stressor:Gut_Shannon              1.01867  2    0.60090  
# Micro_Trtmt:Stressor:Gut_Shannon  1.91961  2    0.38297  

faith_baseline_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_PD, data = df_baseline, family = "gaussian")
Anova(faith_baseline_glmm)
# Response: RC1
#                             LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                   0.9926  1    0.31911  
# Stressor                      0.2170  2    0.89719  
# Gut_PD                        4.0978  1    0.04294 *
# Micro_Trtmt:Stressor          1.9755  2    0.37242  
# Micro_Trtmt:Gut_PD            2.1171  1    0.14567  
# Stressor:Gut_PD               1.3914  2    0.49874  
# Micro_Trtmt:Stressor:Gut_PD   2.7600  2    0.25157 

##PC2
asv_baseline2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_ASV, data = df_baseline, family = "gaussian")
Anova(asv_baseline2_glmm)
# Response: RC2
#                               LR Chisq Df Pr(>Chisq)   
# Micro_Trtmt                    0.5513  1   0.457789   
# Stressor                       2.5217  2   0.283411   
# Gut_ASV                        6.9500  1   0.008382 **
# Micro_Trtmt:Stressor           0.9481  2   0.622464   
# Micro_Trtmt:Gut_ASV            8.0984  1   0.004430 **
# Stressor:Gut_ASV               2.6112  2   0.271016   
# Micro_Trtmt:Stressor:Gut_ASV   1.5971  2   0.449980 

shannon_baseline2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_Shannon, data = df_baseline, family = "gaussian")
Anova(shannon_baseline2_glmm)
# Response: RC2
#                                  LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                       0.52250  1     0.4698
# Stressor                          1.71491  2     0.4242
# Gut_Shannon                       0.05866  1     0.8086
# Micro_Trtmt:Stressor              1.16021  2     0.5598
# Micro_Trtmt:Gut_Shannon           2.05123  1     0.1521
# Stressor:Gut_Shannon              0.91992  2     0.6313
# Micro_Trtmt:Stressor:Gut_Shannon  0.70152  2     0.7042

faith_baseline2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_PD, data = df_baseline, family = "gaussian")
Anova(faith_baseline2_glmm)
# Response: RC2
#                             LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                   0.0726  1    0.78758  
# Stressor                      1.1110  2    0.57380  
# Gut_PD                        4.6212  1    0.03158 *
# Micro_Trtmt:Stressor          0.0252  2    0.98746  
# Micro_Trtmt:Gut_PD            5.4264  1    0.01983 *
# Stressor:Gut_PD               2.1115  2    0.34793  
# Micro_Trtmt:Stressor:Gut_PD   1.7568  2    0.41544  

####Alpha Diversity - Visual Empty
file.choose()
df_empty <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota Behavior\\Behavior Analysis\\Combined data sheets\\Data Sheets with Factor Scores - No Change from Baseline\\Aim3visualemptynochangefactorscores.csv")
df_empty$Micro_Trtmt = factor(df_empty$Micro_Trtmt)
df_empty$Stressor = factor(df_empty$Stressor)

##PC1
asv_empty1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_ASV, data = df_empty, family = "gaussian")
Anova(asv_empty1_glmm)
# Response: RC1
#                             LR Chisq Df Pr(>Chisq)   
# Micro_Trtmt                    0.0283  1   0.866317   
# Stressor                       1.6077  2   0.447608   
# Gut_ASV                        0.7704  1   0.380100   
# Micro_Trtmt:Stressor           2.1430  2   0.342496   
# Micro_Trtmt:Gut_ASV            7.8846  1   0.004986 **
# Stressor:Gut_ASV               2.5565  2   0.278523   
# Micro_Trtmt:Stressor:Gut_ASV   3.5076  2   0.173114 

shannon_empty1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_Shannon, data = df_empty, family = "gaussian")
Anova(shannon_empty1_glmm)
# Response: RC1
#                                 LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                        0.0478  1    0.82689  
# Stressor                           1.1294  2    0.56852  
# Gut_Shannon                        1.9490  1    0.16269  
# Micro_Trtmt:Stressor               1.3834  2    0.50071  
# Micro_Trtmt:Gut_Shannon            3.6612  1    0.05569 .
# Stressor:Gut_Shannon               2.4313  2    0.29652  
# Micro_Trtmt:Stressor:Gut_Shannon   1.5112  2    0.46973  

faith_empty1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_PD, data = df_empty, family = "gaussian")
Anova(faith_empty1_glmm)
# Response: RC1
#                             LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                   0.1317  1     0.7167  
# Stressor                      0.8593  2     0.6507  
# Gut_PD                        0.6158  1     0.4326  
# Micro_Trtmt:Stressor          0.4794  2     0.7869  
# Micro_Trtmt:Gut_PD            3.6911  1     0.0547 .
# Stressor:Gut_PD               1.0069  2     0.6045  
# Micro_Trtmt:Stressor:Gut_PD   1.4859  2     0.4757 

##PC2
asv_empty2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_ASV, data = df_empty, family = "gaussian")
Anova(asv_empty2_glmm)
# Response: RC2
#                              LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                   0.10035  1     0.7514
# Stressor                      0.75472  2     0.6857
# Gut_ASV                       0.13833  1     0.7099
# Micro_Trtmt:Stressor          0.81461  2     0.6654
# Micro_Trtmt:Gut_ASV           0.47149  1     0.4923
# Stressor:Gut_ASV              0.84947  2     0.6539
# Micro_Trtmt:Stressor:Gut_ASV  0.03209  2     0.9841


shannon_empty2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_Shannon, data = df_empty, family = "gaussian")
Anova(shannon_empty2_glmm)
# Response: RC2
#                                   LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                        1.9645  1     0.1610
# Stressor                           1.2916  2     0.5242
# Gut_Shannon                        0.9571  1     0.3279
# Micro_Trtmt:Stressor               3.4963  2     0.1741
# Micro_Trtmt:Gut_Shannon            0.2124  1     0.6449
# Stressor:Gut_Shannon               2.2944  2     0.3175
# Micro_Trtmt:Stressor:Gut_Shannon   0.9671  2     0.6166

faith_empty2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_PD, data = df_empty, family = "gaussian")
Anova(faith_empty2_glmm)
# Response: RC2
#                             LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                  0.27789  1     0.5981
# Stressor                     0.80824  2     0.6676
# Gut_PD                       0.35819  1     0.5495
# Micro_Trtmt:Stressor         0.24758  2     0.8836
# Micro_Trtmt:Gut_PD           0.17164  1     0.6787
# Stressor:Gut_PD              0.14884  2     0.9283
# Micro_Trtmt:Stressor:Gut_PD  0.27278  2     0.8725

##PC3
asv_empty3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Gut_ASV, data = df_empty, family = "gaussian")
Anova(asv_empty3_glmm)
# Response: RC3
#                              LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                    0.9996  1    0.31741  
# Stressor                       4.6727  2    0.09668 .
# Gut_ASV                        2.0957  1    0.14771  
# Micro_Trtmt:Stressor           0.7289  2    0.69457  
# Micro_Trtmt:Gut_ASV            1.0551  1    0.30434  
# Stressor:Gut_ASV               0.4158  2    0.81229  
# Micro_Trtmt:Stressor:Gut_ASV   1.8536  2    0.39583 

shannon_empty3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Gut_Shannon, data = df_empty, family = "gaussian")
Anova(shannon_empty3_glmm)
# Response: RC3
#                                  LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                        0.7417  1    0.38913  
# Stressor                           6.7397  2    0.03439 *
# Gut_Shannon                        0.2199  1    0.63912  
# Micro_Trtmt:Stressor               1.5990  2    0.44955  
# Micro_Trtmt:Gut_Shannon            4.8615  1    0.02746 *
# Stressor:Gut_Shannon               3.5280  2    0.17136  
# Micro_Trtmt:Stressor:Gut_Shannon   1.2976  2    0.52266  

faith_empty3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Gut_PD, data = df_empty, family = "gaussian")
Anova(faith_empty3_glmm)
# Response: RC3
#                              LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                  0.01745  1     0.8949
# Stressor                     2.84466  2     0.2412
# Gut_PD                       0.32719  1     0.5673
# Micro_Trtmt:Stressor         0.98078  2     0.6124
# Micro_Trtmt:Gut_PD           0.50039  1     0.4793
# Stressor:Gut_PD              0.10924  2     0.9468
# Micro_Trtmt:Stressor:Gut_PD  0.94208  2     0.6244

####Alpha Diversity - Visual Food
file.choose()
df_food <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota Behavior\\Behavior Analysis\\Combined data sheets\\Data Sheets with Factor Scores - No Change from Baseline\\Aim3visualfoodnochangefactorscores.csv")
df_food$Micro_Trtmt = factor(df_food$Micro_Trtmt)
df_food$Stressor = factor(df_food$Stressor)

##PC1
asv_food1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_ASV, data = df_food, family = "gaussian")
Anova(asv_food1_glmm)
# Response: RC1
#                               LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                    0.0002  1    0.98758  
# Stressor                       3.9035  2    0.14202  
# Gut_ASV                        1.6543  1    0.19837  
# Micro_Trtmt:Stressor           5.4407  2    0.06585 .
# Micro_Trtmt:Gut_ASV            5.2762  1    0.02162 *
# Stressor:Gut_ASV               5.0353  2    0.08065 .
# Micro_Trtmt:Stressor:Gut_ASV   2.2103  2    0.33115  

shannon_food1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_Shannon, data = df_food, family = "gaussian")
Anova(shannon_food1_glmm)
# Response: RC1
#                                  LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                        0.0948  1    0.75822  
# Stressor                           5.2054  2    0.07407 .
# Gut_Shannon                        1.1077  1    0.29259  
# Micro_Trtmt:Stressor               0.5546  2    0.75784  
# Micro_Trtmt:Gut_Shannon            3.7504  1    0.05280 .
# Stressor:Gut_Shannon               2.3382  2    0.31064  
# Micro_Trtmt:Stressor:Gut_Shannon   2.7678  2    0.25059 

faith_food1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_PD, data = df_food, family = "gaussian")
Anova(faith_food1_glmm)
# Response: RC1
#                             LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                   0.0204  1    0.88634  
# Stressor                      4.4354  2    0.10886  
# Gut_PD                        1.4017  1    0.23644  
# Micro_Trtmt:Stressor          7.8153  2    0.02009 *
# Micro_Trtmt:Gut_PD            4.2624  1    0.03896 *
# Stressor:Gut_PD               8.5849  2    0.01367 *
# Micro_Trtmt:Stressor:Gut_PD   0.1724  2    0.91741  

##PC2
asv_food2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_ASV, data = df_food, family = "gaussian")
Anova(asv_food2_glmm)
# Response: RC2
#                               LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                    0.2093  1    0.64735  
# Stressor                       0.0085  2    0.99574  
# Gut_ASV                        0.0869  1    0.76820  
# Micro_Trtmt:Stressor           7.0384  2    0.02962 *
# Micro_Trtmt:Gut_ASV            0.8435  1    0.35840  
# Stressor:Gut_ASV               7.8933  2    0.01932 *
# Micro_Trtmt:Stressor:Gut_ASV   1.3339  2    0.51327 

shannon_food2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_Shannon, data = df_food, family = "gaussian")
Anova(shannon_food2_glmm)
# Response: RC2
#                                  LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                        1.7257  1     0.1890  
# Stressor                           0.1193  2     0.9421  
# Gut_Shannon                        0.8425  1     0.3587  
# Micro_Trtmt:Stressor               3.7304  2     0.1549  
# Micro_Trtmt:Gut_Shannon            1.6165  1     0.2036  
# Stressor:Gut_Shannon               3.6151  2     0.1641  
# Micro_Trtmt:Stressor:Gut_Shannon   6.6153  2     0.0366 *

faith_food2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_PD, data = df_food, family = "gaussian")
Anova(faith_food2_glmm)
# Response: RC2
#                             LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                  0.11051  1     0.7396
# Stressor                     0.03587  2     0.9822
# Gut_PD                       0.83292  1     0.3614
# Micro_Trtmt:Stressor         1.03488  2     0.5960
# Micro_Trtmt:Gut_PD           0.83478  1     0.3609
# Stressor:Gut_PD              1.42871  2     0.4895
# Micro_Trtmt:Stressor:Gut_PD  0.50866  2     0.7754

##PC3
asv_food3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Gut_ASV, data = df_food, family = "gaussian")
Anova(asv_food3_glmm)
# Response: RC3
#                               LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                   0.08109  1     0.7758
# Stressor                      0.55191  2     0.7588
# Gut_ASV                       0.68177  1     0.4090
# Micro_Trtmt:Stressor          3.09326  2     0.2130
# Micro_Trtmt:Gut_ASV           0.95905  1     0.3274
# Stressor:Gut_ASV              2.33871  2     0.3106
# Micro_Trtmt:Stressor:Gut_ASV  0.43562  2     0.8043

shannon_food3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Gut_Shannon, data = df_food, family = "gaussian")
Anova(shannon_food3_glmm)
# Response: RC3
#                                   LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                        4.0576  1    0.04397 *
# Stressor                           0.0921  2    0.95498  
# Gut_Shannon                        4.3017  1    0.03807 *
# Micro_Trtmt:Stressor               1.4708  2    0.47931  
# Micro_Trtmt:Gut_Shannon            0.0653  1    0.79830  
# Stressor:Gut_Shannon               0.4251  2    0.80853  
# Micro_Trtmt:Stressor:Gut_Shannon   1.7357  2    0.41986  

faith_food3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Gut_PD, data = df_food, family = "gaussian")
Anova(faith_food3_glmm)
# Response: RC3
#                             LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                   0.1765  1    0.67436  
# Stressor                      0.4236  2    0.80914  
# Gut_PD                        1.0062  1    0.31582  
# Micro_Trtmt:Stressor          6.8765  2    0.03212 *
# Micro_Trtmt:Gut_PD            2.2687  1    0.13201  
# Stressor:Gut_PD               5.3186  2    0.07000 .
# Micro_Trtmt:Stressor:Gut_PD   0.0423  2    0.97909  

####Alpha Diversity - Olfactory
file.choose()
df_olfactory <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota Behavior\\Behavior Analysis\\Combined data sheets\\Data Sheets with Factor Scores - No Change from Baseline\\Aim3olfactorynochangefactorscores.csv")
df_olfactory$Micro_Trtmt = factor(df_olfactory$Micro_Trtmt)
df_olfactory$Stressor = factor(df_olfactory$Stressor)

##PC1
asv_olfactory1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_ASV, data = df_olfactory, family = "gaussian")
Anova(asv_olfactory1_glmm)
# Response: RC1
#                               LR Chisq Df Pr(>Chisq)   
# Micro_Trtmt                    0.4293  1   0.512344   
# Stressor                       5.2652  2   0.071893 . 
# Gut_ASV                        3.4237  1   0.064267 . 
# Micro_Trtmt:Stressor           5.2921  2   0.070930 . 
# Micro_Trtmt:Gut_ASV            8.1156  1   0.004389 **
# Stressor:Gut_ASV               5.5342  2   0.062843 . 
# Micro_Trtmt:Stressor:Gut_ASV   4.1028  2   0.128556   

shannon_olfactory1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_Shannon, data = df_olfactory, family = "gaussian")
Anova(shannon_olfactory1_glmm)
# Response: RC1
#                                  LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                        1.8061  1    0.17897  
# Stressor                           4.8609  2    0.08800 .
# Gut_Shannon                        0.4146  1    0.51963  
# Micro_Trtmt:Stressor               2.8500  2    0.24051  
# Micro_Trtmt:Gut_Shannon            4.1825  1    0.04084 *
# Stressor:Gut_Shannon               3.1367  2    0.20839  
# Micro_Trtmt:Stressor:Gut_Shannon   0.6199  2    0.73349 

faith_olfactory1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Gut_PD, data = df_olfactory, family = "gaussian")
Anova(faith_olfactory1_glmm)
# Response: RC1
#                             LR Chisq Df Pr(>Chisq)   
# Micro_Trtmt                   2.0626  1   0.150956   
# Stressor                      8.8220  2   0.012143 * 
# Gut_PD                        6.2186  1   0.012641 * 
# Micro_Trtmt:Stressor          8.5160  2   0.014151 * 
# Micro_Trtmt:Gut_PD           10.2804  1   0.001345 **
# Stressor:Gut_PD              12.7095  2   0.001738 **
# Micro_Trtmt:Stressor:Gut_PD   6.3429  2   0.041943 *

model_ORC1 <- lm(RC1~Gut_PD, data=df_olfactory)
summary(model_ORC1)
#https://www.statology.org/multiple-linear-regression-r/


ggplot(df_olfactory, aes(x= Gut_PD, y = RC1, color = Micro_Trtmt)) +
  geom_point (aes(shape = Micro_Trtmt, size = 3))+
  geom_smooth(method=lm)+
  coord_cartesian(xlim = c(0, 65)) +
  theme_classic() +
  labs(x = "Faith's Phylogenetic Diversity", y = "Change in Locomotory Activity (PC1): Olfactory") +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 17)) +
  scale_color_manual(values = c("seagreen3", "skyblue3"),
                     labels = c("Natural" , "Autoclaved")) +
  scale_shape_manual(values = c(16,17)) +
  theme(legend.text = element_text(size = 14))


##PC2
asv_olfactory2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_ASV, data = df_olfactory, family = "gaussian")
Anova(asv_olfactory2_glmm)
# Response: RC2
#                               LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                   0.00953  1     0.9222
# Stressor                      0.41485  2     0.8127
# Gut_ASV                       0.09284  1     0.7606
# Micro_Trtmt:Stressor          2.73008  2     0.2554
# Micro_Trtmt:Gut_ASV           0.00377  1     0.9511
# Stressor:Gut_ASV              0.89178  2     0.6403
# Micro_Trtmt:Stressor:Gut_ASV  0.95779  2     0.6195

shannon_olfactory2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_Shannon, data = df_olfactory, family = "gaussian")
Anova(shannon_olfactory2_glmm)
# Response: RC2
#                                   LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                       0.21874  1     0.6400
# Stressor                          0.10987  2     0.9465
# Gut_Shannon                       0.30453  1     0.5811
# Micro_Trtmt:Stressor              1.85596  2     0.3954
# Micro_Trtmt:Gut_Shannon           0.78044  1     0.3770
# Stressor:Gut_Shannon              0.20895  2     0.9008
# Micro_Trtmt:Stressor:Gut_Shannon  1.91445  2     0.3840

faith_olfactory2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Gut_PD, data = df_olfactory, family = "gaussian")
Anova(faith_olfactory2_glmm)
# Response: RC2
#                              LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                   0.0785  1     0.7793
# Stressor                      0.6386  2     0.7267
# Gut_PD                        0.0650  1     0.7987
# Micro_Trtmt:Stressor          3.6846  2     0.1585
# Micro_Trtmt:Gut_PD            0.5690  1     0.4507
# Stressor:Gut_PD               0.9544  2     0.6205
# Micro_Trtmt:Stressor:Gut_PD   0.9597  2     0.6189


####Skin Microbiome - Physiology Associations

####Body Mass ~ Alpha Diversity
asv_skinbodymass_glmm <- glm(Body_Mass~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asv_skinbodymass_glmm)

#Analysis of Deviance Table (Type II tests)

# Response: Body_Mass
# LR Chisq Df Pr(>Chisq)   
# Water_Trtmt                     2.8554  1   0.091070 . 
# Stressor                        2.5725  2   0.276311   
# Skin_ASV                        0.0050  1   0.943629   
# Water_Trtmt:Stressor            4.4601  2   0.107522   
# Water_Trtmt:Skin_ASV            6.3612  1   0.011664 * 
# Stressor:Skin_ASV               9.8320  2   0.007328 **
# Water_Trtmt:Stressor:Skin_ASV   1.3786  2   0.501937   

shannonskin_bodymass_glmm <- glm(Body_Mass~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannonskin_bodymass_glmm)

# Response: Body_Mass
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                        1.76966  1     0.1834
# Stressor                           1.28335  2     0.5264
# Skin_Shannon                       0.44685  1     0.5038
# Water_Trtmt:Stressor               0.06008  2     0.9704
# Water_Trtmt:Skin_Shannon           0.77353  1     0.3791
# Stressor:Skin_Shannon              1.99624  2     0.3686
# Water_Trtmt:Stressor:Skin_Shannon  1.87140  2     0.3923

faithskin_bodymass_glmm <- glm(Body_Mass~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faithskin_bodymass_glmm)

# Response: Body_Mass
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                    0.5674  1     0.4513
# Stressor                       1.6159  2     0.4458
# Skin_PD                        0.0798  1     0.7775
# Water_Trtmt:Stressor           0.6103  2     0.7370
# Water_Trtmt:Skin_PD            1.7107  1     0.1909
# Stressor:Skin_PD               1.7122  2     0.4248
# Water_Trtmt:Stressor:Skin_PD   3.4150  2     0.1813

####Gosner Stage ~ Alpha Diversity
asvskin_gosner_glmm <- glm(Gosner_Stage~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asvskin_gosner_glmm)

# Response: Gosner_Stage
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                    1.24981  1     0.2636
# Stressor                       0.61812  2     0.7341
# Skin_ASV                       0.10791  1     0.7425
# Water_Trtmt:Stressor           0.89404  2     0.6395
# Water_Trtmt:Skin_ASV           0.00378  1     0.9510
# Stressor:Skin_ASV              0.73488  2     0.6925
# Water_Trtmt:Stressor:Skin_ASV  0.03207  2     0.9841

shannonskin_gosner_glmm <- glm(Gosner_Stage~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannonskin_gosner_glmm)

# Response: Gosner_Stage
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                        1.07494  1     0.2998
# Stressor                           0.97778  2     0.6133
# Skin_Shannon                       2.70297  1     0.1002
# Water_Trtmt:Stressor               1.90499  2     0.3858
# Water_Trtmt:Skin_Shannon           0.00700  1     0.9333
# Stressor:Skin_Shannon              0.55446  2     0.7579
# Water_Trtmt:Stressor:Skin_Shannon  0.61433  2     0.7355

faithskin_gosner_glmm <- glm(Gosner_Stage~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faithskin_gosner_glmm)

# Response: Gosner_Stage
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                   0.27202  1     0.6020
# Stressor                      1.00317  2     0.6056
# Skin_PD                       0.00655  1     0.9355
# Water_Trtmt:Stressor          1.61009  2     0.4471
# Water_Trtmt:Skin_PD           0.63459  1     0.4257
# Stressor:Skin_PD              0.53592  2     0.7649
# Water_Trtmt:Stressor:Skin_PD  0.64032  2     0.7260

####Relative Brain Mass ~ Alpha Diversity
asvskin_brainmass_glmm <- glm(MAV_Brain_Mass~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asvskin_brainmass_glmm)

# Response: MAV_Brain_Mass
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                    0.41828  1     0.5178
# Stressor                       0.58908  2     0.7449
# Skin_ASV                       0.04525  1     0.8316
# Water_Trtmt:Stressor           1.66273  2     0.4355
# Water_Trtmt:Skin_ASV           0.03726  1     0.8469
# Stressor:Skin_ASV              0.96476  2     0.6173
# Water_Trtmt:Stressor:Skin_ASV  1.02395  2     0.5993

shannonskin_brainmass_glmm <- glm(MAV_Brain_Mass~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannonskin_brainmass_glmm)

# Response: MAV_Brain_Mass
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                        2.06516  1     0.1507
# Stressor                           0.13231  2     0.9360
# Skin_Shannon                       2.47195  1     0.1159
# Water_Trtmt:Stressor               0.59576  2     0.7424
# Water_Trtmt:Skin_Shannon           0.54286  1     0.4612
# Stressor:Skin_Shannon              0.83098  2     0.6600
# Water_Trtmt:Stressor:Skin_Shannon  1.76066  2     0.4146

faithskin_brainmass_glmm <- glm(MAV_Brain_Mass~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faithskin_brainmass_glmm)

# Response: MAV_Brain_Mass
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                    0.7820  1     0.3765
# Stressor                       0.9169  2     0.6323
# Skin_PD                        0.2584  1     0.6112
# Water_Trtmt:Stressor           0.5337  2     0.7658
# Water_Trtmt:Skin_PD            0.2400  1     0.6242
# Stressor:Skin_PD               0.3144  2     0.8545
# Water_Trtmt:Stressor:Skin_PD   3.2468  2     0.1972

##Alpha Diversity - Brain Shape (PC1)
asvskin_brainshape1_glmm <- glm(BrainShape_RC1~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asvskin_brainshape1_glmm)

# Response: BrainShape_RC1
# LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                     2.0670  1    0.15052  
# Stressor                        4.3520  2    0.11349  
# Skin_ASV                        0.6372  1    0.42473  
# Water_Trtmt:Stressor            4.4577  2    0.10765  
# Water_Trtmt:Skin_ASV            0.0895  1    0.76486  
# Stressor:Skin_ASV               4.8178  2    0.08991 .
# Water_Trtmt:Stressor:Skin_ASV   2.5301  2    0.28223  

shannonskin_brainshape1_glmm <- glm(BrainShape_RC1~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannonskin_brainshape1_glmm)

# Response: BrainShape_RC1
# LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                         2.2195  1    0.13628  
# Stressor                            3.9433  2    0.13923  
# Skin_Shannon                        0.3654  1    0.54554  
# Water_Trtmt:Stressor                1.7862  2    0.40939  
# Water_Trtmt:Skin_Shannon            0.2026  1    0.65264  
# Stressor:Skin_Shannon               5.9459  2    0.05115 .
# Water_Trtmt:Stressor:Skin_Shannon   0.0917  2    0.95519 

faithskin_brainshape1_glmm <- glm(BrainShape_RC1~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faithskin_brainshape1_glmm)

# Response: BrainShape_RC1
# LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                    3.2505  1    0.07140 .
# Stressor                       6.8675  2    0.03227 *
# Skin_PD                        0.7413  1    0.38926  
# Water_Trtmt:Stressor           5.4852  2    0.06440 .
# Water_Trtmt:Skin_PD            0.1347  1    0.71357  
# Stressor:Skin_PD               6.0977  2    0.04741 *
# Water_Trtmt:Stressor:Skin_PD   3.4642  2    0.17691  

##Alpha Diversity - Brain Shape (PC2)
asvskin_brainshape2_glmm <- glm(BrainShape_RC2~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asvskin_brainshape2_glmm)

# Response: BrainShape_RC2
# LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                     0.2377  1    0.62585  
# Stressor                        1.3725  2    0.50346  
# Skin_ASV                        0.3068  1    0.57966  
# Water_Trtmt:Stressor            0.6824  2    0.71093  
# Water_Trtmt:Skin_ASV            2.8018  1    0.09416 .
# Stressor:Skin_ASV               0.9447  2    0.62355  
# Water_Trtmt:Stressor:Skin_ASV   4.6679  2    0.09691 .

shannonskin_brainshape2_glmm <- glm(BrainShape_RC2~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannonskin_brainshape2_glmm)

# Response: BrainShape_RC2
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                        0.06611  1     0.7971
# Stressor                           1.92257  2     0.3824
# Skin_Shannon                       0.00050  1     0.9821
# Water_Trtmt:Stressor               0.11410  2     0.9445
# Water_Trtmt:Skin_Shannon           0.18965  1     0.6632
# Stressor:Skin_Shannon              0.82946  2     0.6605
# Water_Trtmt:Stressor:Skin_Shannon  0.65546  2     0.7206

faithskin_brainshape2_glmm <- glm(BrainShape_RC2~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faithskin_brainshape2_glmm)

# Response: BrainShape_RC2
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                   0.05823  1     0.8093
# Stressor                      1.67896  2     0.4319
# Skin_PD                       0.01374  1     0.9067
# Water_Trtmt:Stressor          0.19491  2     0.9071
# Water_Trtmt:Skin_PD           0.54776  1     0.4592
# Stressor:Skin_PD              0.17937  2     0.9142
# Water_Trtmt:Stressor:Skin_PD  0.50122  2     0.7783

##Alpha Diversity - Brain Shape (PC3)
asvskin_brainshape3_glmm <- glm(BrainShape_RC3~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asvskin_brainshape3_glmm)

# Response: BrainShape_RC3
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                     0.3731  1     0.5413
# Stressor                        3.3036  2     0.1917
# Skin_ASV                        0.0562  1     0.8126
# Water_Trtmt:Stressor            0.7650  2     0.6821
# Water_Trtmt:Skin_ASV            0.0210  1     0.8849
# Stressor:Skin_ASV               0.6495  2     0.7227
# Water_Trtmt:Stressor:Skin_ASV   0.1186  2     0.9424

shannonskin_brainshape3_glmm <- glm(BrainShape_RC3~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannonskin_brainshape3_glmm)

# Response: BrainShape_RC3
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                         1.5102  1    0.21910  
# Stressor                            8.6185  2    0.01344 *
# Skin_Shannon                        1.9887  1    0.15848  
# Water_Trtmt:Stressor                0.7627  2    0.68293  
# Water_Trtmt:Skin_Shannon            0.0276  1    0.86800  
# Stressor:Skin_Shannon               6.7516  2    0.03419 *
# Water_Trtmt:Stressor:Skin_Shannon   1.1005  2    0.57680  

faithskin_brainshape3_glmm <- glm(BrainShape_RC3~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faithskin_brainshape3_glmm)

# Response: BrainShape_RC3
# LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                    1.2627  1    0.26115  
# Stressor                       4.7068  2    0.09505 .
# Skin_PD                        0.5810  1    0.44593  
# Water_Trtmt:Stressor           0.4361  2    0.80407  
# Water_Trtmt:Skin_PD            0.4085  1    0.52274  
# Stressor:Skin_PD               0.0304  2    0.98492  
# Water_Trtmt:Stressor:Skin_PD   1.9670  2    0.37400  

####Alpha Diversity - Baseline Behavior
file.choose()
df_baseline <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota Behavior\\Behavior Analysis\\Combined data sheets\\Data Sheets with Factor Scores - No Change from Baseline\\Aim3baselinefactorscores.csv")
df_baseline$Micro_Trtmt = factor(df_baseline$Micro_Trtmt)
df_baseline$Stressor = factor(df_baseline$Stressor)

##PC1
asvskin_baseline_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_ASV, data = df_baseline, family = "gaussian")
Anova(asvskin_baseline_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                    0.03134  1     0.8595
# Stressor                       0.19921  2     0.9052
# Skin_ASV                       0.36847  1     0.5438
# Micro_Trtmt:Stressor           1.96474  2     0.3744
# Micro_Trtmt:Skin_ASV           1.70828  1     0.1912
# Stressor:Skin_ASV              0.98762  2     0.6103
# Micro_Trtmt:Stressor:Skin_ASV  2.37197  2     0.3054

shannonskin_baseline_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_Shannon, data = df_baseline, family = "gaussian")
Anova(shannonskin_baseline_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                         1.0024  1    0.31673  
# Stressor                            0.1394  2    0.93266  
# Skin_Shannon                        2.0596  1    0.15125  
# Micro_Trtmt:Stressor                1.4473  2    0.48499  
# Micro_Trtmt:Skin_Shannon            6.5080  1    0.01074 *
# Stressor:Skin_Shannon               0.1041  2    0.94930  
# Micro_Trtmt:Stressor:Skin_Shannon   3.3535  2    0.18698 

faithskin_baseline_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_PD, data = df_baseline, family = "gaussian")
Anova(faithskin_baseline_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                    0.0043  1     0.9475
# Stressor                       0.9013  2     0.6372
# Skin_PD                        1.8169  1     0.1777
# Micro_Trtmt:Stressor           4.1910  2     0.1230
# Micro_Trtmt:Skin_PD            1.5660  1     0.2108
# Stressor:Skin_PD               1.3021  2     0.5215
# Micro_Trtmt:Stressor:Skin_PD   2.0731  2     0.3547

asvskin_baseline2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_ASV, data = df_baseline, family = "gaussian")
Anova(asvskin_baseline2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     0.2756  1     0.5996  
# Stressor                        1.7902  2     0.4086  
# Skin_ASV                        3.7188  1     0.0538 .
# Micro_Trtmt:Stressor            1.5541  2     0.4598  
# Micro_Trtmt:Skin_ASV            0.0719  1     0.7886  
# Stressor:Skin_ASV               2.8438  2     0.2413  
# Micro_Trtmt:Stressor:Skin_ASV   0.3268  2     0.8493  

shannonskin_baseline2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_Shannon, data = df_baseline, family = "gaussian")
Anova(shannonskin_baseline2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                         1.6721  1   0.195971   
# Stressor                            1.7224  2   0.422648   
# Skin_Shannon                        0.7017  1   0.402204   
# Micro_Trtmt:Stressor                0.8230  2   0.662651   
# Micro_Trtmt:Skin_Shannon            7.2842  1   0.006956 **
# Stressor:Skin_Shannon               4.2826  2   0.117500   
# Micro_Trtmt:Stressor:Skin_Shannon   1.4941  2   0.473766 

faithskin_baseline2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_PD, data = df_baseline, family = "gaussian")
Anova(faithskin_baseline2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                   0.01910  1     0.8901
# Stressor                      0.84250  2     0.6562
# Skin_PD                       2.55092  1     0.1102
# Micro_Trtmt:Stressor          0.20612  2     0.9021
# Micro_Trtmt:Skin_PD           0.18458  1     0.6675
# Stressor:Skin_PD              0.07305  2     0.9641
# Micro_Trtmt:Stressor:Skin_PD  0.73410  2     0.6928

####Alpha Diversity - Visual Empty
file.choose()
df_empty <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota Behavior\\Behavior Analysis\\Combined data sheets\\Data Sheets with Factor Scores - No Change from Baseline\\Aim3visualemptynochangefactorscores.csv")
df_empty$Micro_Trtmt = factor(df_empty$Micro_Trtmt)
df_empty$Stressor = factor(df_empty$Stressor)

##PC1
asvskin_empty1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_ASV, data = df_empty, family = "gaussian")
Anova(asvskin_empty1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     0.4704  1    0.49279  
# Stressor                        1.4813  2    0.47680  
# Skin_ASV                        0.1236  1    0.72515  
# Micro_Trtmt:Stressor            2.9217  2    0.23204  
# Micro_Trtmt:Skin_ASV            1.1065  1    0.29284  
# Stressor:Skin_ASV               2.7781  2    0.24931  
# Micro_Trtmt:Stressor:Skin_ASV   5.9728  2    0.05047 .

shannonskin_empty1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_Shannon, data = df_empty, family = "gaussian")
Anova(shannonskin_empty1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                         2.3599  1    0.12449  
# Stressor                            3.0803  2    0.21435  
# Skin_Shannon                        0.0000  1    0.99630  
# Micro_Trtmt:Stressor                1.6219  2    0.44443  
# Micro_Trtmt:Skin_Shannon            3.1641  1    0.07527 .
# Stressor:Skin_Shannon               2.7698  2    0.25035  
# Micro_Trtmt:Stressor:Skin_Shannon   1.7445  2    0.41800  

faithskin_empty1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_PD, data = df_empty, family = "gaussian")
Anova(faithskin_empty1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                    1.9466  1    0.16295  
# Stressor                       2.8322  2    0.24266  
# Skin_PD                        0.1147  1    0.73490  
# Micro_Trtmt:Stressor           6.0991  2    0.04738 *
# Micro_Trtmt:Skin_PD            3.4692  1    0.06252 .
# Stressor:Skin_PD               6.2477  2    0.04399 *
# Micro_Trtmt:Stressor:Skin_PD   1.2404  2    0.53783  

##PC2
asvskin_empty2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_ASV, data = df_empty, family = "gaussian")
Anova(asvskin_empty2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                    0.94335  1     0.3314
# Stressor                       0.62786  2     0.7306
# Skin_ASV                       0.19043  1     0.6626
# Micro_Trtmt:Stressor           0.51421  2     0.7733
# Micro_Trtmt:Skin_ASV           0.02107  1     0.8846
# Stressor:Skin_ASV              0.53357  2     0.7658
# Micro_Trtmt:Stressor:Skin_ASV  0.32718  2     0.8491

shannonskin_empty2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_Shannon, data = df_empty, family = "gaussian")
Anova(shannonskin_empty2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                        0.70962  1     0.3996
# Stressor                           0.45892  2     0.7950
# Skin_Shannon                       0.02006  1     0.8874
# Micro_Trtmt:Stressor               0.33048  2     0.8477
# Micro_Trtmt:Skin_Shannon           0.07272  1     0.7874
# Stressor:Skin_Shannon              0.19705  2     0.9062
# Micro_Trtmt:Stressor:Skin_Shannon  0.00592  2     0.9970

faithskin_empty2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_PD, data = df_empty, family = "gaussian")
Anova(faithskin_empty2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                    0.1161  1    0.73335  
# Stressor                       0.1198  2    0.94186  
# Skin_PD                        1.3154  1    0.25141  
# Micro_Trtmt:Stressor           5.8045  2    0.05490 .
# Micro_Trtmt:Skin_PD            0.6265  1    0.42862  
# Stressor:Skin_PD               8.5344  2    0.01402 *
# Micro_Trtmt:Stressor:Skin_PD   0.2788  2    0.86987  

##PC3
asvskin_empty3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_ASV, data = df_empty, family = "gaussian")
Anova(asvskin_empty3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     0.1345  1    0.71380  
# Stressor                        5.8944  2    0.05249 .
# Skin_ASV                        0.0266  1    0.87034  
# Micro_Trtmt:Stressor            0.5083  2    0.77556  
# Micro_Trtmt:Skin_ASV            1.5392  1    0.21473  
# Stressor:Skin_ASV               0.4946  2    0.78091  
# Micro_Trtmt:Stressor:Skin_ASV   4.1365  2    0.12641  

shannonskin_empty3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_Shannon, data = df_empty, family = "gaussian")
Anova(shannonskin_empty3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                         0.0758  1    0.78303  
# Stressor                            7.1604  2    0.02787 *
# Skin_Shannon                        1.5354  1    0.21530  
# Micro_Trtmt:Stressor                0.8223  2    0.66290  
# Micro_Trtmt:Skin_Shannon            0.3364  1    0.56191  
# Stressor:Skin_Shannon               7.0218  2    0.02987 *
# Micro_Trtmt:Stressor:Skin_Shannon   2.7862  2    0.24831  

faithskin_empty3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_PD, data = df_empty, family = "gaussian")
Anova(faithskin_empty3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                    0.0749  1    0.78430  
# Stressor                       4.0579  2    0.13147  
# Skin_PD                        0.3599  1    0.54859  
# Micro_Trtmt:Stressor           5.3525  2    0.06882 .
# Micro_Trtmt:Skin_PD            4.9764  1    0.02570 *
# Stressor:Skin_PD               4.4663  2    0.10719  
# Micro_Trtmt:Stressor:Skin_PD   0.5885  2    0.74510 

####Alpha Diversity - Visual Food
file.choose()
df_food <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota Behavior\\Behavior Analysis\\Combined data sheets\\Data Sheets with Factor Scores - No Change from Baseline\\Aim3visualfoodnochangefactorscores.csv")
df_food$Micro_Trtmt = factor(df_food$Micro_Trtmt)
df_food$Stressor = factor(df_food$Stressor)

##PC1
asvskin_food1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_ASV, data = df_food, family = "gaussian")
Anova(asvskin_food1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     0.4086  1    0.52267  
# Stressor                        2.1620  2    0.33926  
# Skin_ASV                        0.0177  1    0.89420  
# Micro_Trtmt:Stressor            0.8226  2    0.66278  
# Micro_Trtmt:Skin_ASV            3.8959  1    0.04841 *
# Stressor:Skin_ASV               0.4182  2    0.81133  
# Micro_Trtmt:Stressor:Skin_ASV   3.1650  2    0.20546  

shannonskin_food1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_Shannon, data = df_food, family = "gaussian")
Anova(shannonskin_food1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                         2.4372  1    0.11849  
# Stressor                            3.9376  2    0.13962  
# Skin_Shannon                        0.1662  1    0.68354  
# Micro_Trtmt:Stressor                1.7759  2    0.41149  
# Micro_Trtmt:Skin_Shannon            3.0014  1    0.08319 .
# Stressor:Skin_Shannon               2.8354  2    0.24227  
# Micro_Trtmt:Stressor:Skin_Shannon   1.2501  2    0.53523 

faithskin_food1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_PD, data = df_food, family = "gaussian")
Anova(faithskin_food1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                    0.8844  1    0.34700  
# Stressor                       3.6063  2    0.16478  
# Skin_PD                        0.3488  1    0.55479  
# Micro_Trtmt:Stressor           5.5865  2    0.06122 .
# Micro_Trtmt:Skin_PD            4.9640  1    0.02588 *
# Stressor:Skin_PD               3.9133  2    0.14133  
# Micro_Trtmt:Stressor:Skin_PD   0.6630  2    0.71784  

##PC2
asvskin_food2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_ASV, data = df_food, family = "gaussian")
Anova(asvskin_food2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                    0.06032  1     0.8060
# Stressor                       0.06016  2     0.9704
# Skin_ASV                       0.07382  1     0.7859
# Micro_Trtmt:Stressor           0.55585  2     0.7574
# Micro_Trtmt:Skin_ASV           0.34102  1     0.5592
# Stressor:Skin_ASV              0.76826  2     0.6810
# Micro_Trtmt:Stressor:Skin_ASV  1.66629  2     0.4347

shannonskin_food2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_Shannon, data = df_food, family = "gaussian")
Anova(shannonskin_food2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                        0.06886  1     0.7930
# Stressor                           0.01899  2     0.9905
# Skin_Shannon                       0.41776  1     0.5181
# Micro_Trtmt:Stressor               0.33998  2     0.8437
# Micro_Trtmt:Skin_Shannon           0.00276  1     0.9581
# Stressor:Skin_Shannon              0.49299  2     0.7815
# Micro_Trtmt:Stressor:Skin_Shannon  0.44851  2     0.7991

faithskin_food2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_PD, data = df_food, family = "gaussian")
Anova(faithskin_food2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                    0.1382  1     0.7101
# Stressor                       0.0733  2     0.9640
# Skin_PD                        0.0021  1     0.9637
# Micro_Trtmt:Stressor           3.5948  2     0.1657
# Micro_Trtmt:Skin_PD            1.9320  1     0.1645
# Stressor:Skin_PD               3.1529  2     0.2067
# Micro_Trtmt:Stressor:Skin_PD   1.3908  2     0.4989

##PC3
asvskin_food3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_ASV, data = df_food, family = "gaussian")
Anova(asvskin_food3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                     0.2755  1     0.5997
# Stressor                        0.3033  2     0.8593
# Skin_ASV                        0.0287  1     0.8655
# Micro_Trtmt:Stressor            0.1847  2     0.9118
# Micro_Trtmt:Skin_ASV            0.1194  1     0.7297
# Stressor:Skin_ASV               0.3982  2     0.8195
# Micro_Trtmt:Stressor:Skin_ASV   4.1955  2     0.1227

shannonskin_food3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_Shannon, data = df_food, family = "gaussian")
Anova(shannonskin_food3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                         0.0016  1    0.96781  
# Stressor                            0.1631  2    0.92167  
# Skin_Shannon                        4.4667  1    0.03456 *
# Micro_Trtmt:Stressor                4.6266  2    0.09893 .
# Micro_Trtmt:Skin_Shannon            1.0875  1    0.29703  
# Stressor:Skin_Shannon               5.8714  2    0.05309 .
# Micro_Trtmt:Stressor:Skin_Shannon   1.0369  2    0.59543

faithskin_food3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_PD, data = df_food, family = "gaussian")
Anova(faithskin_food3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                   0.00153  1     0.9688
# Stressor                      0.13279  2     0.9358
# Skin_PD                       0.00174  1     0.9667
# Micro_Trtmt:Stressor          2.21339  2     0.3306
# Micro_Trtmt:Skin_PD           0.15871  1     0.6903
# Stressor:Skin_PD              1.50677  2     0.4708
# Micro_Trtmt:Stressor:Skin_PD  1.33870  2     0.5120

####Alpha Diversity - Olfactory
file.choose()
df_olfactory <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota Behavior\\Behavior Analysis\\Combined data sheets\\Data Sheets with Factor Scores - No Change from Baseline\\Aim3olfactorynochangefactorscores.csv")
df_olfactory$Micro_Trtmt = factor(df_olfactory$Micro_Trtmt)
df_olfactory$Stressor = factor(df_olfactory$Stressor)

##PC1
asvskin_olfactory1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_ASV, data = df_olfactory, family = "gaussian")
Anova(asvskin_olfactory1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     5.9176  1    0.01499 *
# Stressor                        2.6726  2    0.26282  
# Skin_ASV                        0.0651  1    0.79855  
# Micro_Trtmt:Stressor            2.8808  2    0.23684  
# Micro_Trtmt:Skin_ASV            1.9169  1    0.16620  
# Stressor:Skin_ASV               2.0452  2    0.35967  
# Micro_Trtmt:Stressor:Skin_ASV   2.9102  2    0.23338  

shannonskin_olfactory1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_Shannon, data = df_olfactory, family = "gaussian")
Anova(shannonskin_olfactory1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)   
# Micro_Trtmt                        10.8063  1   0.001012 **
# Stressor                            3.1009  2   0.212157   
# Skin_Shannon                        0.0365  1   0.848462   
# Micro_Trtmt:Stressor                0.1395  2   0.932614   
# Micro_Trtmt:Skin_Shannon            0.3924  1   0.531030   
# Stressor:Skin_Shannon               1.3531  2   0.508376   
# Micro_Trtmt:Stressor:Skin_Shannon   4.1753  2   0.123977   

faithskin_olfactory1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_PD, data = df_olfactory, family = "gaussian")
Anova(faithskin_olfactory1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)   
# Micro_Trtmt                    8.6722  1   0.003231 **
# Stressor                       3.8335  2   0.147081   
# Skin_PD                        0.0017  1   0.966731   
# Micro_Trtmt:Stressor           1.2436  2   0.536977   
# Micro_Trtmt:Skin_PD            1.2651  1   0.260684   
# Stressor:Skin_PD               2.0806  2   0.353356   
# Micro_Trtmt:Stressor:Skin_PD   3.8261  2   0.147633 

##PC2
asvskin_olfactory2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_ASV, data = df_olfactory, family = "gaussian")
Anova(asvskin_olfactory2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     3.3155  1    0.06863 .
# Stressor                        0.7177  2    0.69848  
# Skin_ASV                        1.1924  1    0.27485  
# Micro_Trtmt:Stressor            0.6061  2    0.73856  
# Micro_Trtmt:Skin_ASV            0.1071  1    0.74343  
# Stressor:Skin_ASV               0.3978  2    0.81964  
# Micro_Trtmt:Stressor:Skin_ASV   2.1578  2    0.33996  

shannonskin_olfactory2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_Shannon, data = df_olfactory, family = "gaussian")
Anova(shannonskin_olfactory2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                         1.0671  1    0.30161  
# Stressor                            0.0757  2    0.96287  
# Skin_Shannon                        2.1300  1    0.14444  
# Micro_Trtmt:Stressor                6.6718  2    0.03558 *
# Micro_Trtmt:Skin_Shannon            1.9144  1    0.16648  
# Stressor:Skin_Shannon               5.0685  2    0.07932 .
# Micro_Trtmt:Stressor:Skin_Shannon   6.8215  2    0.03302 *

faithskin_olfactory2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_PD, data = df_olfactory, family = "gaussian")
Anova(faithskin_olfactory2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                    4.4213  1    0.03549 *
# Stressor                       0.9044  2    0.63622  
# Skin_PD                        3.6546  1    0.05591 .
# Micro_Trtmt:Stressor           6.3140  2    0.04255 *
# Micro_Trtmt:Skin_PD            3.2137  1    0.07302 .
# Stressor:Skin_PD               2.9939  2    0.22381  
# Micro_Trtmt:Stressor:Skin_PD   4.4761  2    0.10666  
