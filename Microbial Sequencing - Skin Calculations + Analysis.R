##Microbial Sequencing - Calculations and Analysis
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

##### Alpha Diversity
## ASVs

file.choose()
df_asv <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Microbiome- Skin and Water Excel Files\\Emerson Aim 3 - no. observed skin asvs.csv")
df_asv$Microbial_Trtmt = factor(df_asv$Microbial_Trtmt)
df_asv$Stressor = factor(df_asv$Stressor)
df_asv$Trtmt_Combo = factor(df_asv$Trtmt_Combo)

asvs_glmm <- glm(Obs_ASVs~Microbial_Trtmt*Stressor, data = df_asv, family = "gaussian")
Anova(asvs_glmm)
# Response: Obs_ASVs
#                          LR Chisq Df Pr(>Chisq)    
# Microbial_Trtmt           21.3478  1  3.831e-06 ***
# Stressor                   5.3882  2     0.0676 .  
# Microbial_Trtmt:Stressor   1.1700  2     0.5571   

summary(df_asv)

Stress <- c("Vehicle Control", "Predator Cues", "CORT")
#For labeling X axis

ggplot(df_asv, aes(x= factor(Stressor, level=c('B', 'A','C')), y = Obs_ASVs, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "No. Observed ASVs") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = c(.83, .92)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position="none")

## Shannon Diversity Index
file.choose()
df_shannon <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Microbiome- Skin and Water Excel Files\\Emerson Aim 3 - shannon skin diversity.csv")
df_shannon$Microbial_Trtmt = factor(df_shannon$Microbial_Trtmt)
df_shannon$Stressor = factor(df_shannon$Stressor)
df_shannon$Trtmt_Combo = factor(df_shannon$Trtmt_Combo)

shannon_glmm <- glm(shannon~Microbial_Trtmt*Stressor, data = df_shannon, family = "gaussian")
Anova(shannon_glmm)
# Response: shannon
#                          LR Chisq Df Pr(>Chisq)
# Microbial_Trtmt            0.1515  1     0.6971
# Stressor                   1.1763  2     0.5554
# Microbial_Trtmt:Stressor   3.4738  2     0.1761


ggplot(df_shannon, aes(x= factor(Stressor, level=c('B', 'A','C')), y = shannon, fill = Microbial_Trtmt)) +
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


## Faiths Phylogenetic Diversity
file.choose()
df_faith <- read.csv("C:\\R\\Aim 3 - Stress and the Microbiota\\Microbiome- Skin and Water Excel Files\\Emerson Aim 3 - Final skin FaithsPD.csv")
df_faith$Microbial_Trtmt = factor(df_faith$Microbial_Trtmt)
df_faith$Stressor = factor(df_faith$Stressor)
df_faith$Trtmt_Combo = factor(df_faith$Trtmt_Combo)

faith_glmm <- glm(PD~Microbial_Trtmt*Stressor, data = df_faith, family = "gaussian")
Anova(faith_glmm)
# Response: PD
#                         LR Chisq   Df Pr(>Chisq)    
# Microbial_Trtmt           13.6122  1  0.0002247 ***
# Stressor                   4.2135  2  0.1216343    
# Microbial_Trtmt:Stressor   1.2318  2  0.5401481   

ggplot(df_faith, aes(x= factor(Stressor, level=c('B', 'A','C')), y = PD, fill = Microbial_Trtmt)) +
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

####Skin Microbiome - Physiology Associations
file.choose()
df <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota\\Aim 3 (Restarted)\\Aim 3 metadata.csv")
df$Water_Trtmt = factor(df$Water_Trtmt)
df$Stressor = factor(df$Stressor)
df$Trtmt_Combo = factor(df$Trtmt_Combo)
df$Replicate = factor(df$Replicate)
df$Gosner_Stage = as.numeric(df$Gosner_Stage)

####Body Mass ~ Alpha Diversity
asv_bodymass_glmm <- glm(Body_Mass~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asv_bodymass_glmm)

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

shannon_bodymass_glmm <- glm(Body_Mass~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannon_bodymass_glmm)

# Response: Body_Mass
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                        1.89189  1     0.1690
# Stressor                           0.43961  2     0.8027
# Skin_Shannon                       0.02688  1     0.8698
# Water_Trtmt:Stressor               0.76922  2     0.6807
# Water_Trtmt:Skin_Shannon           0.75974  1     0.3834
# Stressor:Skin_Shannon              0.70376  2     0.7034
# Water_Trtmt:Stressor:Skin_Shannon  0.32484  2     0.8501

faith_bodymass_glmm <- glm(Body_Mass~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faith_bodymass_glmm)

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
asv_gosner_glmm <- glm(Gosner_Stage~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asv_gosner_glmm)

# Response: Gosner_Stage
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                    1.24981  1     0.2636
# Stressor                       0.61812  2     0.7341
# Skin_ASV                       0.10791  1     0.7425
# Water_Trtmt:Stressor           0.89404  2     0.6395
# Water_Trtmt:Skin_ASV           0.00378  1     0.9510
# Stressor:Skin_ASV              0.73488  2     0.6925
# Water_Trtmt:Stressor:Skin_ASV  0.03207  2     0.9841

shannon_gosner_glmm <- glm(Gosner_Stage~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannon_gosner_glmm)

# Response: Gosner_Stage
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                        0.51538  1     0.4728
# Stressor                           0.50848  2     0.7755
# Skin_Shannon                       0.51236  1     0.4741
# Water_Trtmt:Stressor               2.66237  2     0.2642
# Water_Trtmt:Skin_Shannon           0.00201  1     0.9643
# Stressor:Skin_Shannon              0.14257  2     0.9312
# Water_Trtmt:Stressor:Skin_Shannon  2.64639  2     0.2663

faith_gosner_glmm <- glm(Gosner_Stage~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faith_gosner_glmm)

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
asv_brainmass_glmm <- glm(MAV_Brain_Mass~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asv_brainmass_glmm)

# Response: MAV_Brain_Mass
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                    0.41828  1     0.5178
# Stressor                       0.58908  2     0.7449
# Skin_ASV                       0.04525  1     0.8316
# Water_Trtmt:Stressor           1.66273  2     0.4355
# Water_Trtmt:Skin_ASV           0.03726  1     0.8469
# Stressor:Skin_ASV              0.96476  2     0.6173
# Water_Trtmt:Stressor:Skin_ASV  1.02395  2     0.5993

shannon_brainmass_glmm <- glm(MAV_Brain_Mass~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannon_brainmass_glmm)

# Response: MAV_Brain_Mass
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                         0.0226  1     0.8805
# Stressor                            2.1781  2     0.3365
# Skin_Shannon                        1.3810  1     0.2399
# Water_Trtmt:Stressor                2.0131  2     0.3655
# Water_Trtmt:Skin_Shannon            0.3903  1     0.5321
# Stressor:Skin_Shannon               4.4522  2     0.1079
# Water_Trtmt:Stressor:Skin_Shannon   0.2179  2     0.8968

faith_brainmass_glmm <- glm(MAV_Brain_Mass~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faith_brainmass_glmm)

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
asv_brainshape1_glmm <- glm(BrainShape_RC1~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asv_brainshape1_glmm)

# Response: BrainShape_RC1
# LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                     2.0670  1    0.15052  
# Stressor                        4.3520  2    0.11349  
# Skin_ASV                        0.6372  1    0.42473  
# Water_Trtmt:Stressor            4.4577  2    0.10765  
# Water_Trtmt:Skin_ASV            0.0895  1    0.76486  
# Stressor:Skin_ASV               4.8178  2    0.08991 .
# Water_Trtmt:Stressor:Skin_ASV   2.5301  2    0.28223  

shannon_brainshape1_glmm <- glm(BrainShape_RC1~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannon_brainshape1_glmm)

# Response: BrainShape_RC1
# LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                         2.6867  1    0.10119  
# Stressor                            3.5832  2    0.16669  
# Skin_Shannon                        0.0821  1    0.77449  
# Water_Trtmt:Stressor                3.0407  2    0.21864  
# Water_Trtmt:Skin_Shannon            0.2192  1    0.63964  
# Stressor:Skin_Shannon               3.5469  2    0.16974  
# Water_Trtmt:Stressor:Skin_Shannon   5.2689  2    0.07176 .

faith_brainshape1_glmm <- glm(BrainShape_RC1~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faith_brainshape1_glmm)

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
asv_brainshape2_glmm <- glm(BrainShape_RC2~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asv_brainshape2_glmm)

# Response: BrainShape_RC2
# LR Chisq Df Pr(>Chisq)  
# Water_Trtmt                     0.2377  1    0.62585  
# Stressor                        1.3725  2    0.50346  
# Skin_ASV                        0.3068  1    0.57966  
# Water_Trtmt:Stressor            0.6824  2    0.71093  
# Water_Trtmt:Skin_ASV            2.8018  1    0.09416 .
# Stressor:Skin_ASV               0.9447  2    0.62355  
# Water_Trtmt:Stressor:Skin_ASV   4.6679  2    0.09691 .

shannon_brainshape2_glmm <- glm(BrainShape_RC2~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannon_brainshape2_glmm)

# Response: BrainShape_RC2
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                        0.09577  1     0.7570
# Stressor                           1.75188  2     0.4165
# Skin_Shannon                       0.63689  1     0.4248
# Water_Trtmt:Stressor               0.02515  2     0.9875
# Water_Trtmt:Skin_Shannon           1.07237  1     0.3004
# Stressor:Skin_Shannon              0.17974  2     0.9140
# Water_Trtmt:Stressor:Skin_Shannon  0.19212  2     0.9084

faith_brainshape2_glmm <- glm(BrainShape_RC2~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faith_brainshape2_glmm)

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
asv_brainshape3_glmm <- glm(BrainShape_RC3~Water_Trtmt*Stressor*Skin_ASV, data = df, family = "gaussian")
Anova(asv_brainshape3_glmm)

# Response: BrainShape_RC3
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                     0.3731  1     0.5413
# Stressor                        3.3036  2     0.1917
# Skin_ASV                        0.0562  1     0.8126
# Water_Trtmt:Stressor            0.7650  2     0.6821
# Water_Trtmt:Skin_ASV            0.0210  1     0.8849
# Stressor:Skin_ASV               0.6495  2     0.7227
# Water_Trtmt:Stressor:Skin_ASV   0.1186  2     0.9424

shannon_brainshape3_glmm <- glm(BrainShape_RC3~Water_Trtmt*Stressor*Skin_Shannon, data = df, family = "gaussian")
Anova(shannon_brainshape3_glmm)

# Response: BrainShape_RC3
# LR Chisq Df Pr(>Chisq)
# Water_Trtmt                         0.0271  1     0.8692
# Stressor                            4.0475  2     0.1322
# Skin_Shannon                        0.1332  1     0.7152
# Water_Trtmt:Stressor                0.1177  2     0.9429
# Water_Trtmt:Skin_Shannon            0.2138  1     0.6438
# Stressor:Skin_Shannon               1.5031  2     0.4716
# Water_Trtmt:Stressor:Skin_Shannon   0.4970  2     0.7800

faith_brainshape3_glmm <- glm(BrainShape_RC3~Water_Trtmt*Stressor*Skin_PD, data = df, family = "gaussian")
Anova(faith_brainshape3_glmm)

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
asv_baseline_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_ASV, data = df_baseline, family = "gaussian")
Anova(asv_baseline_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                    0.03134  1     0.8595
# Stressor                       0.19921  2     0.9052
# Skin_ASV                       0.36847  1     0.5438
# Micro_Trtmt:Stressor           1.96474  2     0.3744
# Micro_Trtmt:Skin_ASV           1.70828  1     0.1912
# Stressor:Skin_ASV              0.98762  2     0.6103
# Micro_Trtmt:Stressor:Skin_ASV  2.37197  2     0.3054

shannon_baseline_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_Shannon, data = df_baseline, family = "gaussian")
Anova(shannon_baseline_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                        0.79268  1     0.3733
# Stressor                           1.18701  2     0.5524
# Skin_Shannon                       1.01126  1     0.3146
# Micro_Trtmt:Stressor               2.43222  2     0.2964
# Micro_Trtmt:Skin_Shannon           0.53222  1     0.4657
# Stressor:Skin_Shannon              0.22951  2     0.8916
# Micro_Trtmt:Stressor:Skin_Shannon  0.01860  2     0.9907

faith_baseline_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_PD, data = df_baseline, family = "gaussian")
Anova(faith_baseline_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                    0.0043  1     0.9475
# Stressor                       0.9013  2     0.6372
# Skin_PD                        1.8169  1     0.1777
# Micro_Trtmt:Stressor           4.1910  2     0.1230
# Micro_Trtmt:Skin_PD            1.5660  1     0.2108
# Stressor:Skin_PD               1.3021  2     0.5215
# Micro_Trtmt:Stressor:Skin_PD   2.0731  2     0.3547

asv_baseline2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_ASV, data = df_baseline, family = "gaussian")
Anova(asv_baseline2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     0.2756  1     0.5996  
# Stressor                        1.7902  2     0.4086  
# Skin_ASV                        3.7188  1     0.0538 .
# Micro_Trtmt:Stressor            1.5541  2     0.4598  
# Micro_Trtmt:Skin_ASV            0.0719  1     0.7886  
# Stressor:Skin_ASV               2.8438  2     0.2413  
# Micro_Trtmt:Stressor:Skin_ASV   0.3268  2     0.8493  

shannon_baseline2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_Shannon, data = df_baseline, family = "gaussian")
Anova(shannon_baseline2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                        0.71382  1     0.3982
# Stressor                           0.30272  2     0.8595
# Skin_Shannon                       0.24531  1     0.6204
# Micro_Trtmt:Stressor               0.60841  2     0.7377
# Micro_Trtmt:Skin_Shannon           0.00166  1     0.9675
# Stressor:Skin_Shannon              0.32818  2     0.8487
# Micro_Trtmt:Stressor:Skin_Shannon  1.84209  2     0.3981

faith_baseline2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_PD, data = df_baseline, family = "gaussian")
Anova(faith_baseline2_glmm)

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
asv_empty1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_ASV, data = df_empty, family = "gaussian")
Anova(asv_empty1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     0.4704  1    0.49279  
# Stressor                        1.4813  2    0.47680  
# Skin_ASV                        0.1236  1    0.72515  
# Micro_Trtmt:Stressor            2.9217  2    0.23204  
# Micro_Trtmt:Skin_ASV            1.1065  1    0.29284  
# Stressor:Skin_ASV               2.7781  2    0.24931  
# Micro_Trtmt:Stressor:Skin_ASV   5.9728  2    0.05047 .

shannon_empty1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_Shannon, data = df_empty, family = "gaussian")
Anova(shannon_empty1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                         3.5645  1    0.05903 .
# Stressor                            1.0548  2    0.59014  
# Skin_Shannon                        1.3154  1    0.25142  
# Micro_Trtmt:Stressor                3.0364  2    0.21911  
# Micro_Trtmt:Skin_Shannon            0.7005  1    0.40260  
# Stressor:Skin_Shannon               1.6205  2    0.44474  
# Micro_Trtmt:Stressor:Skin_Shannon   0.5807  2    0.74800  

faith_empty1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_PD, data = df_empty, family = "gaussian")
Anova(faith_empty1_glmm)

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
asv_empty2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_ASV, data = df_empty, family = "gaussian")
Anova(asv_empty2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                    0.94335  1     0.3314
# Stressor                       0.62786  2     0.7306
# Skin_ASV                       0.19043  1     0.6626
# Micro_Trtmt:Stressor           0.51421  2     0.7733
# Micro_Trtmt:Skin_ASV           0.02107  1     0.8846
# Stressor:Skin_ASV              0.53357  2     0.7658
# Micro_Trtmt:Stressor:Skin_ASV  0.32718  2     0.8491

shannon_empty2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_Shannon, data = df_empty, family = "gaussian")
Anova(shannon_empty2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                         1.4230  1     0.2329  
# Stressor                            0.6417  2     0.7255  
# Skin_Shannon                        0.0093  1     0.9233  
# Micro_Trtmt:Stressor                1.5732  2     0.4554  
# Micro_Trtmt:Skin_Shannon            0.5143  1     0.4733  
# Stressor:Skin_Shannon               2.9360  2     0.2304  
# Micro_Trtmt:Stressor:Skin_Shannon   9.1133  2     0.0105 *

faith_empty2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_PD, data = df_empty, family = "gaussian")
Anova(faith_empty2_glmm)

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
asv_empty3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_ASV, data = df_empty, family = "gaussian")
Anova(asv_empty3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     0.1345  1    0.71380  
# Stressor                        5.8944  2    0.05249 .
# Skin_ASV                        0.0266  1    0.87034  
# Micro_Trtmt:Stressor            0.5083  2    0.77556  
# Micro_Trtmt:Skin_ASV            1.5392  1    0.21473  
# Stressor:Skin_ASV               0.4946  2    0.78091  
# Micro_Trtmt:Stressor:Skin_ASV   4.1365  2    0.12641  

shannon_empty3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_Shannon, data = df_empty, family = "gaussian")
Anova(shannon_empty3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                         0.1761  1     0.6747
# Stressor                            4.4926  2     0.1058
# Skin_Shannon                        2.5544  1     0.1100
# Micro_Trtmt:Stressor                1.4636  2     0.4810
# Micro_Trtmt:Skin_Shannon            0.1007  1     0.7510
# Stressor:Skin_Shannon               0.7819  2     0.6764
# Micro_Trtmt:Stressor:Skin_Shannon   0.5283  2     0.7679

faith_empty3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_PD, data = df_empty, family = "gaussian")
Anova(faith_empty3_glmm)

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
asv_food1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_ASV, data = df_food, family = "gaussian")
Anova(asv_food1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     0.4086  1    0.52267  
# Stressor                        2.1620  2    0.33926  
# Skin_ASV                        0.0177  1    0.89420  
# Micro_Trtmt:Stressor            0.8226  2    0.66278  
# Micro_Trtmt:Skin_ASV            3.8959  1    0.04841 *
# Stressor:Skin_ASV               0.4182  2    0.81133  
# Micro_Trtmt:Stressor:Skin_ASV   3.1650  2    0.20546  

shannon_food1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_Shannon, data = df_food, family = "gaussian")
Anova(shannon_food1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                        2.32600  1     0.1272
# Stressor                           1.56448  2     0.4574
# Skin_Shannon                       0.09823  1     0.7540
# Micro_Trtmt:Stressor               2.13779  2     0.3434
# Micro_Trtmt:Skin_Shannon           1.35454  1     0.2445
# Stressor:Skin_Shannon              0.27231  2     0.8727
# Micro_Trtmt:Stressor:Skin_Shannon  0.49354  2     0.7813

faith_food1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_PD, data = df_food, family = "gaussian")
Anova(faith_food1_glmm)

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
asv_food2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_ASV, data = df_food, family = "gaussian")
Anova(asv_food2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                    0.06032  1     0.8060
# Stressor                       0.06016  2     0.9704
# Skin_ASV                       0.07382  1     0.7859
# Micro_Trtmt:Stressor           0.55585  2     0.7574
# Micro_Trtmt:Skin_ASV           0.34102  1     0.5592
# Stressor:Skin_ASV              0.76826  2     0.6810
# Micro_Trtmt:Stressor:Skin_ASV  1.66629  2     0.4347

shannon_food2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_Shannon, data = df_food, family = "gaussian")
Anova(shannon_food2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                        0.63192  1     0.4267
# Stressor                           0.08292  2     0.9594
# Skin_Shannon                       0.18488  1     0.6672
# Micro_Trtmt:Stressor               1.67604  2     0.4326
# Micro_Trtmt:Skin_Shannon           0.38174  1     0.5367
# Stressor:Skin_Shannon              1.69799  2     0.4278
# Micro_Trtmt:Stressor:Skin_Shannon  0.33148  2     0.8473

faith_food2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_PD, data = df_food, family = "gaussian")
Anova(faith_food2_glmm)

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
asv_food3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_ASV, data = df_food, family = "gaussian")
Anova(asv_food3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                     0.2755  1     0.5997
# Stressor                        0.3033  2     0.8593
# Skin_ASV                        0.0287  1     0.8655
# Micro_Trtmt:Stressor            0.1847  2     0.9118
# Micro_Trtmt:Skin_ASV            0.1194  1     0.7297
# Stressor:Skin_ASV               0.3982  2     0.8195
# Micro_Trtmt:Stressor:Skin_ASV   4.1955  2     0.1227

shannon_food3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_Shannon, data = df_food, family = "gaussian")
Anova(shannon_food3_glmm)

# Response: RC3
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                        0.14738  1     0.7011
# Stressor                           0.52144  2     0.7705
# Skin_Shannon                       0.12192  1     0.7270
# Micro_Trtmt:Stressor               0.68403  2     0.7103
# Micro_Trtmt:Skin_Shannon           0.39309  1     0.5307
# Stressor:Skin_Shannon              0.68751  2     0.7091
# Micro_Trtmt:Stressor:Skin_Shannon  0.15849  2     0.9238

faith_food3_glmm <- glm(RC3~Micro_Trtmt*Stressor*Skin_PD, data = df_food, family = "gaussian")
Anova(faith_food3_glmm)

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
asv_olfactory1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_ASV, data = df_olfactory, family = "gaussian")
Anova(asv_olfactory1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     5.9176  1    0.01499 *
# Stressor                        2.6726  2    0.26282  
# Skin_ASV                        0.0651  1    0.79855  
# Micro_Trtmt:Stressor            2.8808  2    0.23684  
# Micro_Trtmt:Skin_ASV            1.9169  1    0.16620  
# Stressor:Skin_ASV               2.0452  2    0.35967  
# Micro_Trtmt:Stressor:Skin_ASV   2.9102  2    0.23338  

shannon_olfactory1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_Shannon, data = df_olfactory, family = "gaussian")
Anova(shannon_olfactory1_glmm)

# Response: RC1
# LR Chisq Df Pr(>Chisq)   
# Micro_Trtmt                         8.6260  1   0.003314 **
# Stressor                            2.6560  2   0.265005   
# Skin_Shannon                        0.2792  1   0.597251   
# Micro_Trtmt:Stressor                0.2303  2   0.891223   
# Micro_Trtmt:Skin_Shannon            0.6876  1   0.406965   
# Stressor:Skin_Shannon               0.5872  2   0.745587   
# Micro_Trtmt:Stressor:Skin_Shannon   3.9926  2   0.135836   

faith_olfactory1_glmm <- glm(RC1~Micro_Trtmt*Stressor*Skin_PD, data = df_olfactory, family = "gaussian")
Anova(faith_olfactory1_glmm)

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
asv_olfactory2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_ASV, data = df_olfactory, family = "gaussian")
Anova(asv_olfactory2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                     3.3155  1    0.06863 .
# Stressor                        0.7177  2    0.69848  
# Skin_ASV                        1.1924  1    0.27485  
# Micro_Trtmt:Stressor            0.6061  2    0.73856  
# Micro_Trtmt:Skin_ASV            0.1071  1    0.74343  
# Stressor:Skin_ASV               0.3978  2    0.81964  
# Micro_Trtmt:Stressor:Skin_ASV   2.1578  2    0.33996  

shannon_olfactory2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_Shannon, data = df_olfactory, family = "gaussian")
Anova(shannon_olfactory2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)
# Micro_Trtmt                        1.58672  1     0.2078
# Stressor                           0.73653  2     0.6919
# Skin_Shannon                       0.69973  1     0.4029
# Micro_Trtmt:Stressor               1.31561  2     0.5180
# Micro_Trtmt:Skin_Shannon           1.82401  1     0.1768
# Stressor:Skin_Shannon              1.42888  2     0.4895
# Micro_Trtmt:Stressor:Skin_Shannon  2.75905  2     0.2517

faith_olfactory2_glmm <- glm(RC2~Micro_Trtmt*Stressor*Skin_PD, data = df_olfactory, family = "gaussian")
Anova(faith_olfactory2_glmm)

# Response: RC2
# LR Chisq Df Pr(>Chisq)  
# Micro_Trtmt                    4.4213  1    0.03549 *
# Stressor                       0.9044  2    0.63622  
# Skin_PD                        3.6546  1    0.05591 .
# Micro_Trtmt:Stressor           6.3140  2    0.04255 *
# Micro_Trtmt:Skin_PD            3.2137  1    0.07302 .
# Stressor:Skin_PD               2.9939  2    0.22381  
# Micro_Trtmt:Stressor:Skin_PD   4.4761  2    0.10666  