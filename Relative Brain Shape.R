#Relative Brain Shape

library(lme4)
library(car)
library(readr)
library(moments) 
library(psych)
library(pastecs)
library(ggplot2)
library(ggbiplot)
library(tidyverse)
library(modelbased)
library(dplyr)

file.choose()
df <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 3 - CORT & the Microbiota\\Aim 3 (Restarted)\\Relative Brain Shape\\Aim 3 - Relative Brain Shape.csv")
df$Water_Trtmt = factor(df$Water_Trtmt)
df$Stressor = factor(df$Stressor)
df$Trtmt_Combo = factor(df$Trtmt_Combo)
df$Replicate = factor(df$Replicate)

shapiro.test(df$TW)
shapiro.test(df$TL)
shapiro.test(df$OTW)
shapiro.test(df$OTL)
shapiro.test(df$DW)
shapiro.test(df$DL)
shapiro.test(df$MW)


shapiro.test(df$Log_TW)
#Only fail?
shapiro.test(df$Log_TL)
shapiro.test(df$Log_OTW)
shapiro.test(df$Log_OTL)
shapiro.test(df$Log_DW)
shapiro.test(df$Log_DL)
shapiro.test(df$Log_MW)

#Levene x TrtmtCombo
leveneTest(df$TW, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$TL, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$OTW, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$OTL, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$DW, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$DL, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$MW, df$Trtmt_Combo, center = mean, na.rm = TRUE)


leveneTest(df$Log_TW, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$Log_TL, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$Log_OTW, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$Log_OTL, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$Log_DW, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$Log_DL, df$Trtmt_Combo, center = mean, na.rm = TRUE)
leveneTest(df$Log_MW, df$Trtmt_Combo, center = mean, na.rm = TRUE)

#MANCOVA
man1 <- manova(cbind(Log_TW, Log_TL, Log_OTW, Log_OTL, Log_DW, Log_DL, Log_MW)~Water_Trtmt*Stressor*Log_Brain_Mass, data = df)
summary(man1)
summary.aov(man1)
#All pass!

#EMMs
modelTW <- lm(Log_TW ~ Water_Trtmt*Stressor*Log_Brain_Mass, data = df)
means_complexTW <- estimate_means(modelTW)

modelTL <- lm(Log_TL ~ Water_Trtmt*Stressor*Log_Brain_Mass, data = df)
means_complexTL <- estimate_means(modelTL)

modelOTW <- lm(Log_OTW ~ Water_Trtmt*Stressor*Log_Brain_Mass, data = df)
means_complexOTW <- estimate_means(modelOTW)

modelOTL <- lm(Log_OTL ~ Water_Trtmt*Stressor*Log_Brain_Mass, data = df)
means_complexOTL <- estimate_means(modelOTL)

modelDW <- lm(Log_DW ~ Water_Trtmt*Stressor*Log_Brain_Mass, data = df)
means_complexDW <- estimate_means(modelDW)

modelDL <- lm(Log_DL ~ Water_Trtmt*Stressor*Log_Brain_Mass, data = df)
means_complexDL <- estimate_means(modelDL)

modelMW <- lm(Log_MW ~ Water_Trtmt*Stressor*Log_Brain_Mass, data = df)
means_complexMW <- estimate_means(modelMW)

#Residuals
residuals <- man1$residuals
#write.table(residuals, file = "EM3dfwithresidualsshape.csv", sep = ",")

df.pca = df[,40:46]
KMO(df.pca)
cortest.bartlett(df.pca)
#Both pass

pca.p <- principal(df.pca, nfactors = 3, rotate = "varimax")

qplot(c(1:7), pca.p$values) +
  geom_line() +
  xlab("Principal Component") +
  ylab("Eigenvalue") +
  ylim(0,6)
#2 or 3 work

pca.p$values

pca.p$scores

Factor <- pca.p$scores
df <- cbind(df, Factor)
#write.table(df, file = "EM3dfwithfactorscores,2PCs.csv", sep = ",")

pca.p$loadings
#PC1: TW, TL, OTW
#PC2: OTL, DW, MW
#PC3: DL

pca1_glmm <- glmer(RC1~Water_Trtmt*Stressor + (1|Replicate), data = df, na.action = na.omit, family = "gaussian")
Anova(pca1_glmm)
# Response: RC1
#                       Chisq Df Pr(>Chisq)
# Water_Trtmt          0.2072  1     0.6490
# Stressor             3.6133  2     0.1642
# Water_Trtmt:Stressor 2.3303  2     0.3119

emmeans(pca1_glmm, list(pairwise ~ Stressor), adjust = "tukey")
# 1     estimate    SE   df t.ratio p.value
# A - B  -0.0672 0.369 18.2  -0.182  0.9819
# A - C  -0.6369 0.367 17.9  -1.736  0.2197
# B - C  -0.5697 0.367 17.8  -1.553  0.2911
# 
# Results are averaged over the levels of: Water_Trtmt 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates

pca2_glmm <- glmer(RC2~Water_Trtmt*Stressor + (1|Replicate), data = df, na.action = na.omit, family = "gaussian")
Anova(pca2_glmm)
# Response: RC2
#                       Chisq Df Pr(>Chisq)  
# Water_Trtmt          4.5229  1    0.03344 *
# Stressor             9.0789  2    0.01068 *
# Water_Trtmt:Stressor 3.5436  2    0.17003 

emmeans(pca2_glmm, list(pairwise ~ Stressor), adjust = "tukey")
# 1     estimate    SE   df t.ratio p.value
# A - B  -0.0938 0.226 18.1  -0.415  0.9098
# A - C   0.5108 0.223 17.4   2.296  0.0831
# B - C   0.6046 0.221 16.5   2.730  0.0367
# 
# Results are averaged over the levels of: Water_Trtmt 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 

pca3_glmm <- glmer(RC3~Water_Trtmt*Stressor + (1|Replicate), data = df, na.action = na.omit, family = "gaussian")
Anova(pca3_glmm)
# Response: RC3
#                       Chisq Df Pr(>Chisq)
# Water_Trtmt          0.0053  1     0.9422
# Stressor             0.9990  2     0.6068
# Water_Trtmt:Stressor 3.0005  2     0.2231

emmeans(pca3_glmm, list(pairwise ~ Stressor), adjust = "tukey")
# 1     estimate    SE   df t.ratio p.value
# A - B  -0.2508 0.341 18.2  -0.735  0.7461
# A - C   0.0619 0.339 17.8   0.183  0.9818
# B - C   0.3127 0.339 17.7   0.923  0.6334
# 
# Results are averaged over the levels of: Water_Trtmt 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates

#Plots
##PC-1

Stress <- c("Vehicle Control", "Predator Cues", "CORT")
#For labeling X axis

ggplot(df, aes(x= factor(Stressor, level=c('B', 'A','C')), y = RC1, fill = Water_Trtmt)) +
  geom_boxplot(width = 0.5, alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "PC-1") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14))

#PC2
ggplot(df, aes(x= factor(Stressor, level=c('B', 'A','C')), y = RC2, fill = Water_Trtmt)) +
  geom_boxplot(width = 0.5,  alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "PC-2") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14))

#PC2 with y axis shifted

ggplot(df, aes(x= factor(Stressor, level=c('B', 'A','C')), y = RC2, fill = Water_Trtmt)) +
  geom_boxplot(width = 0.5,  alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_y_continuous(limits=c(0, 2)) +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "PC-2") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16, colour = "black"))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14))

#PC3
ggplot(df, aes(x= factor(Stressor, level=c('B', 'A','C')), y = RC3, fill = Water_Trtmt)) +
  geom_boxplot(width = 0.5,  alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "PC-3") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14))

####Now will try analysis with only 2 PC's
##Re-did all analysis above except I changed parameters for 2 PCs
##Carrying the rest on down here

pca.p$loadings
#PC1: TW, TL, OTW
#PC2: OTL, DW, MW

pca1_glmm <- glmer(RC1~Water_Trtmt*Stressor + (1|Replicate), data = df, na.action = na.omit, family = "gaussian")
Anova(pca1_glmm)

pca2_glmm <- glmer(RC2~Water_Trtmt*Stressor + (1|Replicate), data = df, na.action = na.omit, family = "gaussian")
Anova(pca2_glmm)

Stress <- c("Vehicle Control", "Predator Cues", "CORT")
#For labeling X axis

ggplot(df, aes(x= factor(Stressor, level=c('B', 'A','C')), y = RC1, fill = Water_Trtmt)) +
  geom_boxplot(width = 0.5, alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "PC-1") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14))

#PC2
ggplot(df, aes(x= factor(Stressor, level=c('B', 'A','C')), y = RC2, fill = Water_Trtmt)) +
  geom_boxplot(width = 0.5,  alpha =0.35) +
  stat_boxplot(geom = "errorbar", width = .5) +
  theme_classic() +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Microbial Environment",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Stressor", y = "PC-2") +
  theme(aspect.ratio = 1) +
  scale_x_discrete(labels = Stress)+
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 20)) +
  theme(legend.position = "none") +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 14))
