

##################################################################
# [3.1] ANTIOXIDATIVE AFLATOXIN DEGRADATION/MODULATION IN FLAVI  #
##################################################################

# None-, Low-, and High-Aflatoxin Producer Flavi:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit <- read.csv("https://github.com/bkatati/afloxidate/Afla_modulation.csv")
# If error occurs, download csv file and create your own path on PC
head(fit)
dim(fit)

# Table 3 & 4: Pairwise comparison of aflatoxin degradation/modulation:
# ~~~~~~~~~~~~
Mytox <- subset(fit, Mytox == "none")

# NB: for 'low' and 'high' producers, use:
# Mytox <- subset(fit, Mytox == "low")
# Mytox <- subset(fit, Mytox == "high")

print(Mytox)
dim(Mytox)

x <- Mytox$B1
# NB: For Aflatoxin-G1, use "x <- Mytox$G1"
g <- factor(Mytox$Media)
pairwise.wilcox.test(x, g, p.adjust.method = "none", paired = FALSE)

# Figure 2. Boxplot showing change in quantity of aflatoxin-B1 and -G1:
# ~~~~~~~~
require(ggplot2)
y.expression <- expression("B1   ng/mL, × 10" ^	- 3) # symbol for minus is " - " (with "Windows 1252" encoding format)
boxplot(B1~Media, data = Mytox, xlab = "", las = 0, ylab = "", main = "Change in Aflatoxin B1 by Non-producer Flavi"); title(ylab = y.expression, line = 2.3); title(xlab = "Treatment", line = 2.5)

# Gene Expression due to Antioxidant:
# ***********************************

Genx <- read.csv("https://github.com/bkatati/afloxidate/GenOx.csv")
# If error occurs, download csv file and create your own path on PC
Genox <- subset(Genx, Strain == "M14F")
# NB: for isolate "E33C" use "Genox <- subset(Genx, Strain == "E33C")"
print(Genox)

# Dataframe:
Strain <- Genox$Strain
Gene <- Genox$Gene
Antioxidant <- Genox$AOx
RGE_aflD <- Genox$RGE_AflD
RGE_AflR <- Genox$RGE_AflR
AfB1 <- Genox$AfB1
AfG1 <- Genox$AfG1
gendf <- data.frame(Strain=Strain, AOx=Antioxidant, RGE_aflD = RGE_aflD, RGE_AflR = RGE_AflR, AFB1=AfB1,AFG1=AfG1)
print(gendf)

# Modified dataframe:
library(reshape2)
library(ggplot2)

modf <- melt(gendf, id.vars='AOx',
             measure.vars=c('RGE_aflD', 'RGE_AflR', 'AFB1','AFG1'))
head(modf)

p <- ggplot(modf) +
  geom_boxplot(aes(x=AOx, y=value, fill=variable)) + ggtitle("Isolate MLV14F structural (aflD) and regulatory (aflR) pathway
genes expression and aflatoxin production response to antioxidant") + theme_classic() +
  labs(x= "Antioxidant (Se) dose", y = "measure", fill="Legend", grouping("Gene")) + ylim(0, 17) + theme(plot.title = element_text(size = 13), axis.text = element_text(size = 13),text = element_text(size = 12)) # face = "bold")) #,face = "bold"))


# Plot
print(p)

# Spearman's Correlation For AflD and AflR Expression
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(Genx)
corr <- cor.test(x=Genx$RGE_AflD, y=Genx$RGE_AflR, method = "spearman")
corr

# Spearman's rank correlation rho:
# p-value < 0.001; Rho = 0.94, indicating strong correlation.

##############################################################
# [3.2] ROLE OF ANTIOXIDANT (Se) ON FUNGAL DARWINIAN FITNESS # 
##############################################################

# Data:
fit <- read.csv("https://github.com/bkatati/afloxidate/Fitness_antioxidant.csv")
# If error occurs, download csv file and create your own path on PC
head(fit)
dim(fit)

# Figure 4: Boxplot of change in Flavi fitness due to antioxidant (Se): 
# ~~~~~~~~
library(ggplot2)
y.expression <- expression("(Fitness)  Spores/mL, × 10" ^ - 6) # symbol for minus is " - " (with "Windows 1252" encoding format)
ggplot(fit, aes(x = Se, y = cfu, col = Variant)) + theme_classic() + geom_boxplot(size=0.6, alpha=1) +
  labs(title = "  Role of antioxidant (Se) in fungal Darwinian fitness", x = "Antioxidant Level", y = y.expression) + ylim(0, 3) + theme_classic() +
  theme(axis.title.x = element_text(), axis.title.y = element_text(),
        panel.grid.major = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.ticks = element_blank())

# Pairwise test, toxigenic vs atoxigenic Flavi due to antioxidant:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(fit)
dim(fit)
treat <- subset(fit, Se == "0.86 mg/kg")

# NB:
# for other two lower antioxidant concentrations (0.4 mg/kg and 0 mg/kg) use:
# treat <- subset(fit, Se == "0.40 mg/kg")
# treat <- subset(fit, Se == "0.0 mg/kg")

print(treat)
dim(treat)
atox <- treat$Variant=="atoxigenic"
tox <- treat$Variant=="toxigenic"
wilcox.test(treat$cfu[atox], treat$cfu[tox], paired=F)

# Spore count, sp/ml (geometric mean):

exp((mean(treat$lnCfu[atox])))

# ii) Toxigenic:
exp((mean(treat$lnCfu[tox])))

# Higher fitness of atoxigenic strains (geometric mean, 0.655 x 10???6 spores/ml)
# compared to toxigenic strains (geometric mean, 0.209 x 10^-6 spores/ml) due to antioxidant treatment.
# difference at 0 and 0.4 mg/kg antioxidant not significant (p = 0.74 and 0.91).

