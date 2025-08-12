# Astrangia poculata heat response
# Code for statistical analysis and reproducing plots for main text & supplement

# Load packages:
library(ggplot2)
library(csv)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(readr)
library(Rmisc)
library(lme4)
library(lmerTest)
library(rstatix)
library(gridExtra)
library(data.table)
library(car)
library(patchwork)
library(mgcv)
library(scales)
library(dplyr)
library(devtools)
library(emmeans)

####### Contents ####### 
## I. Tank temperatures and visual observations (Figures 1-2)
## II. Animal responses:
# 1. Plot of initial (ambient genet) symbiont density (for Figure 1)
# 2. Histology plots & stats (for Figure 3)
# 3. Temperature reaction norms (panels and stats) (for Figure 4)
# 4. Regressions: ambient genet symbiont density vs. ∆ physiology variables when heated (panels & stats) (for Figure 5)
# 5. Regression: ambient genet symbiont density vs. ∆ chlorophyll fluorescence (for Figure S2)
## III. pHi standard (for figure S1)

####### Temperatures and visual observations  ####### 

# Summary function
summarySE <- function(data = NULL,
                      measurevar,
                      groupvars = NULL,
                      na.rm=FALSE,
                      conf.interval=.95,
                      .drop=TRUE) {
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)
  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Temperature
# Load temperature data:
tempdata <- read.csv(file.choose(),
                    header = T)
# Stats:
tempdata$Treatment <- factor(tempdata$Treatment,
                             levels = c("Ambient",
                                        "High"))
tempdata$Tub <- factor(tempdata$Tub,
                       levels = c("1",
                                  "2",
                                  "3"))
tempdata$DateTime <- as.POSIXct(tempdata$DateTime)
templm <- lm(Temperature...C. ~ Treatment*DateTime,
             data = tempdata)
templmanova <- Anova(templm,
                     Type = "III")
templmanova

tempmeans <- emmeans(templm,
                     list(pairwise~Treatment),
                     adjust = "tukey")
tempmeans
# Create daily temp averages:
dailyavgdata <- tempdata %>%
  mutate(Date = floor_date(DateTime, "day")) %>%
  group_by(Treatment, Date) %>%
  dplyr::summarize(Daily.average.temp = mean(Temperature...C.),
                   sd = sd(Temperature...C.),
                   n = n(),
                   se = sd / sqrt(n))
dailyavgdata$Day <- rep(1:19, times = 2)
# Plot temps for the first 18 days:
dailyavgdata <- dailyavgdata %>% subset(Day<19)
dailyavgplot <- ggplot(data = dailyavgdata,
                       aes(x = Day,
                           y = Daily.average.temp,
                           color = Treatment,
                           fill = Treatment)) +
  geom_point() +
  geom_errorbar(aes(ymin = Daily.average.temp - se,
                    ymax = Daily.average.temp + se),
                width = 0.5) +
  geom_line() +
  xlab("Day") +
  ylab("Temperature (°C)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        legend.position = c(0.11, 0.77),
        legend.key = element_rect(fill = NA)) +
  scale_color_manual(values = c("Ambient" = "#0571B0",
                                "High" = "#F4A582")) +
  scale_fill_manual(values = c("Ambient" = "#0571B0",
                               "High" = "#F4A582")) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18)) +
  scale_y_continuous(breaks = c(21,23,25,27,29,31)) +
  coord_cartesian(ylim = c(21,31), xlim = c(0,18))
dailyavgplot

## Visual observations
# Load visual observation data:
dodata <- read.csv("ApocDailyObservations.csv")

## Stats
dodata$Treatment <- factor(dodata$Treatment,
                           levels = c("Ambient",
                                      "High"))
dodata$Tub <- factor(dodata$Tub,
                     levels = c("1",
                                "2",
                                "3"))
dodata$Colony <- factor(dodata$Colony)

# Extension statistics
extensionlm <- lmer(Extension ~ Treatment*Day + (1|Colony) + (1|Tub),
                    data = dodata)
extensionanova <- Anova(extensionlm,
                        Type = "III")
extensionanova
## treatment, day, and their interaction are all significant
# pairwise analysis:
dodata$Day <- as.factor(dodata$Day)
ext.lm <- lm(Extension ~ Treatment*Day, data = dodata)
ext.tukey <- emmeans(ext.lm, list(pairwise ~ Treatment*Day), simple = "Treatment", adjust = "tukey")
ext.tukey
# Significantly different on days 3, 4, 7, 10, 17

# Color statistics
colorlm <- lmer(Color ~ Treatment*Day + (1|Colony) + (1|Tub:Treatment),
                data = dodata)
coloranova <- Anova(colorlm,
                    Type = "III")
coloranova
## date and treatment x date are significant, but not treatment by itself
# pairwise analysis:
color.lm <- lm(Color ~ Treatment*Day, data = dodata)
color.tukey <- emmeans(color.lm, list(pairwise ~ Treatment*Day), simple = "Treatment", adjust = "tukey")
color.tukey
# Significantly different on days 13-17

## Plot polyp extension
# Create data summary
extensionsum <- summarySE(data = dodata,
                          measurevar = "Extension",
                          groupvars = c("Treatment",
                                        "Day"))
# Plot
extensionplot <- ggplot(data = extensionsum,
                        aes(x = Day,
                            y = Extension,
                            color = Treatment,
                            fill = Treatment)) +
  geom_point() +
  geom_errorbar(aes(ymin = Extension - se,
                    ymax = Extension + se),
                width = 0.5) +
  geom_line() +
  geom_vline(xintercept = 6.5,
             linetype = "dotted",
             color = "black",
             size = 1) +
  xlab("Day") +
  ylab("Polyp Extension (%)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  scale_color_manual(values = c("Ambient" = "#0571B0",
                                "High" = "#F4A582")) +
  scale_fill_manual(values = c("Ambient" = "#0571B0",
                               "High" = "#F4A582")) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18)) +
  coord_cartesian(ylim = c(0,100), xlim = c(0,18)) +
  annotate("text", x = 3, y = 52, size = 5, label = "*") +
  annotate("text", x = 4, y = 29, size = 5, label = "*") +
  annotate("text", x = 7, y = 77, size = 5, label = "*") +
  annotate("text", x = 10, y = 51, size = 5, label = "*") +
  annotate("text", x = 17, y = 20, size = 5, label = "*") 
extensionplot
## Plot color
# Create data summary
colorsum <- summarySE(data = dodata,
                      measurevar = "Color",
                      groupvars = c("Treatment",
                                    "Day"))

#Plot
colorplot <- ggplot(data = colorsum,
                    aes(x = Day,
                        y = Color,
                        color = Treatment,
                        fill = Treatment)) +
  geom_point() +
  geom_errorbar(aes(ymin = Color - se,
                    ymax = Color + se),
                width = 0.5) +
  geom_line() +
  xlab("Day") +
  ylab("Color Score (%)") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(color = "black"), 
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  scale_color_manual(values = c("Ambient" = "#0571B0",
                                "High" = "#F4A582")) +
  scale_fill_manual(values = c("Ambient" = "#0571B0",
                               "High" = "#F4A582")) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16,18)) +
  coord_cartesian(ylim = c(0,100), xlim = c(0,18)) +
  annotate("text", x = 13, y = 90, size = 5, label = "*") +
  annotate("text", x = 14, y = 91, size = 5, label = "*") +
  annotate("text", x = 15, y = 90, size = 5, label = "*") +
  annotate("text", x = 16, y = 90, size = 5, label = "*") +
  annotate("text", x = 17, y = 88, size = 5, label = "*") 
colorplot

### For Figure 2: all temperature, extension, color data
ggarrange(dailyavgplot, extensionplot, colorplot, nrow = 3, ncol = 1)

####### Physiology, histology, and pHi ####### 
# Load phys, histology, and pHi data 
apoc <- read.csv("ApocHeatStress.csv")
amb.apoc <- apoc %>% subset(Treatment == "Ambient")

### For Fig 1
# Symbiont density distribution for ambient-treated colonies (proxy for initial symb density)
Symbs.hist <- ggplot(amb.apoc, aes(x=Symbs.cm2/100000, fill = Treatment)) +
  geom_histogram(alpha=0.6, position = "identity", binwidth=0.5, color="black")+
  labs(y="Count", x=expression(atop("Symbiont Density", paste((10^5~cells~cm^-2)))))+
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 16),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_fill_manual(values = c("#0571B0")) +
  coord_cartesian(ylim = c(0,6))
Symbs.hist
# nonparametric test of symbdens by treatment
wilcox.test(apoc$Symbs.cm2 ~ apoc$Treatment)
# p=0.22

### Histology metrics (Fig 3)

# Symbiont integrity
SymbI <- ggplot(apoc, aes(x=Treatment, y=Symbiont.Health, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = "Symbiont Integrity") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  annotate("text", x = 1.5, y = 4.9, label = "* p=0.036", size = 5)+
  coord_cartesian(ylim = c(0,5))
SymbI
# linear mixed effects model
SymbI.lmer <- lmer(Symbiont.Health ~ Treatment + (1|Colony), data = apoc)
summary(SymbI.lmer) # significant effect of treatment (*p=0.0362)

# Epidermal clarity
Epidermis <- ggplot(apoc, aes(x=Treatment, y=Epidermis.Integrity, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = "Epidermal Clarity") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  annotate("text", x = 1.5, y = 5, label = "n.s.", size = 5)+
  coord_cartesian(ylim = c(0,5))
Epidermis
# linear mixed effects model
Epidermis.lmer <- lmer(Epidermis.Integrity ~ Treatment + (1|Colony), data = apoc)
summary(Epidermis.lmer) # no significant effect of treatment (p=0.398)

# Tissue integrity
TissueInteg <- ggplot(apoc, aes(x=Treatment, y=Tissue.Integrity, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = "Tissue Integrity") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  annotate("text", x = 1.5, y = 4.9, label = "*** p<0.001", size = 5)+
  coord_cartesian(ylim = c(0,5))
TissueInteg
# linear mixed effects model
TissueInteg.lmer <- lmer(Tissue.Integrity ~ Treatment + (1|Colony), data = apoc)
summary(TissueInteg.lmer) # highly sig effect of treatment (***p<0.001)

# Cellular integrity (lack of necrosis & granularity)
Necrosis <- ggplot(apoc, aes(x=Treatment, y=Necrosis.Granularity, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = "Cellular Integrity") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank()) +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  annotate("text", x = 1.5, y = 5, label = "n.s.", size = 5)+
  coord_cartesian(ylim = c(0,5))
Necrosis
# linear mixed effects model
Necrosis.lmer <- lmer(Necrosis.Granularity ~ Treatment + (1|Colony), data = apoc)
summary(Necrosis.lmer) # no significant effect of treatment (p=0.189)

# Epidermal thickness
EpidermThick <- ggplot(apoc, aes(x=Treatment, y=Epi.mm*1000, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = "Epiderm Thickness (µm)", x = "Temperature") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_x_discrete(labels=c("Ambient" = "22°C", "High" = "30°C")) +
  annotate("text", x = 1.5, y = 15.5, label = "p=0.087", size = 5)+
  coord_cartesian(ylim = c(0,16))
EpidermThick
# linear mixed effects model
EpidermThick.lmer <- lmer(Epi.mm ~ Treatment + (1|Colony), data = apoc)
summary(EpidermThick.lmer) # no significant effect of treatment (p=0.087)

# Eggs
eggs <- apoc %>% dplyr::select(Treatment,Colony,Gametes) %>%
  reshape(idvar = "Colony", timevar = "Treatment", direction = "wide") %>% drop_na() %>% dplyr::select(-c("Colony"))
eggs.bar <- barplot(colSums(gametes), col = c("#0571B0","#F4A582"))

# Eggs per square mm
Eggs.per.mm2 <- ggplot(apoc, aes(x=Treatment, y=Eggs.mm2, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = "Eggs mm-2", x = "Temperature") +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582")) +
  scale_x_discrete(labels=c("Ambient" = "22°C", "High" = "30°C")) +
  annotate("text", x = 1.5, y = 62, label = "* p=0.047", size = 5) +
  coord_cartesian(ylim = c(0,63))
Eggs.per.mm2

# All histology metrics patchworked together
# (using TissueInteg to hold space for eggs.bar, which is not a ggplot object and therefore can't be incorporated via patchwork)
(SymbI | TissueInteg) / (Necrosis | Epidermis) / (EpidermThick | (TissueInteg / Eggs.per.mm2))

#### For Figure 4: How did temperature affect phys variables?
# symb density:
# linear mixed model
symbs.lmer <- lmer(Symbs.cm2 ~ Treatment + (1|Colony), data = apoc)
summary(symbs.lmer)
# plot
SymbsByTreatment <- ggplot(apoc, aes(x=Treatment, y=Symbs.cm2/100000, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = expression(atop("Symbiont Density", paste((10^5~cells~cm^-2)))), x = "Temperature") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  annotate("text", x = 1.5, y = 13, label = "n.s.")
SymbsByTreatment

# chlorophyll:
# linear mixed model
chl.lmer <- lmer(rel.chl ~ Treatment + (1|Colony), data = apoc)
summary(chl.lmer)
# plot
ChlByTreatment <- ggplot(apoc, aes(x=Treatment, y=rel.chl, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = expression(atop("Chlorophyll", paste("(fluorescence"~cell^-1~")"))), x = "Temperature") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  annotate("text", x = 1.5, y = 8900, label = "p<0.001")+
  coord_cartesian(ylim=c(0,8900))
ChlByTreatment

# protein:
# linear mixed model
prot.lmer <- lmer(prot.ug.cm2 ~ Treatment + (1|Colony), data = apoc)
summary(prot.lmer)
# plot
ProtByTreatment <- ggplot(apoc, aes(x=Treatment, y=prot.ug.cm2, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = expression(atop("Protein", paste("(µg"~cm^-2~")"))), x = "Temperature") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  annotate("text", x = 1.5, y = 700, label = "p<0.001")+
  coord_cartesian(ylim=c(0,700))
ProtByTreatment

# calcification
# linear mixed model
calc.lmer <- lmer(Calcif.g.cm2.day ~ Treatment + (1|Colony), data = apoc)
summary(calc.lmer)
# plot
CalcByTreatment <- ggplot(apoc, aes(x=Treatment, y=Calcif.g.cm2.day, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = expression(atop("Calcification", paste((g~cm^-2~d^-1)))), x = "Temperature") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582")) +
  annotate("text", x = 1.5, y = 0.0165, label = "p=0.009") +
  coord_cartesian(ylim=c(-0.007,0.017)) +
  geom_hline(yintercept=0,linetype="dashed")
CalcByTreatment

# Symbiocyte pHi
# linear mixed model
symbpHi.lmer <- lmer(Symb.pHi ~ Treatment + (1|Colony), data = apoc)
summary(symbpHi.lmer)
# plot
SymbpHiByTreatment <- ggplot(apoc, aes(x=Treatment, y=Symb.pHi, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = expression(atop("Symbiocyte", paste(Intracellular~pH))), x = "Temperature") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582")) +
  annotate("text", x = 1.5, y = 7.6, label = "n.s.") +
  coord_cartesian(ylim=c(6,7.6))
SymbpHiByTreatment

# Nonsymbiocyte pHi
# linear mixed model
nonsymbpHi.lmer <- lmer(Nonsymb.pHi ~ Treatment + (1|Colony), data = apoc)
summary(nonsymbpHi.lmer)
# plot
NonsymbpHiByTreatment <- ggplot(apoc, aes(x=Treatment, y=Nonsymb.pHi, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = expression(atop("Nonsymbiocyte", paste(Intracellular~pH))), x = "Temperature") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582")) +
  annotate("text", x = 1.5, y = 7.6, label = "n.s.") +
  coord_cartesian(ylim=c(6,7.6))
NonsymbpHiByTreatment

# All together (aspect ratio for export: 7.9 x 5.1 in)
AllRxns <- ggarrange(SymbsByTreatment, ChlByTreatment, ProtByTreatment, 
                     CalcByTreatment, SymbpHiByTreatment, NonsymbpHiByTreatment,
                     nrow = 2, ncol = 3, align = "hv")
AllRxns

#### For Figure 5: Ambient symbiont density vs. change after heating for each variable
apoc.wide <- apoc %>% select(Treatment,Colony,Symbs.cm2,rel.chl,prot.ug.cm2,Symb.pHi,Nonsymb.pHi,Calcif.g.cm2.day,Epi.mm) %>%
  reshape(idvar = "Colony", timevar = "Treatment", direction = "wide")
apoc.wide$delta.Symbs.cm2 <- apoc.wide$Symbs.cm2.High - apoc.wide$Symbs.cm2.Ambient
apoc.wide$delta.relchl <- apoc.wide$rel.chl.High - apoc.wide$rel.chl.Ambient
apoc.wide$delta.Prot <- apoc.wide$prot.ug.cm2.High - apoc.wide$prot.ug.cm2.Ambient
apoc.wide$delta.Symb.pHi <- apoc.wide$Symb.pHi.High - apoc.wide$Symb.pHi.Ambient
apoc.wide$delta.Nonsymb.pHi <- apoc.wide$Nonsymb.pHi.High - apoc.wide$Nonsymb.pHi.Ambient
apoc.wide$delta.Calcif.g.cm2.day <- apoc.wide$Calcif.g.cm2.day.High - apoc.wide$Calcif.g.cm2.day.Ambient
apoc.wide$delta.EpidermalThickness <- apoc.wide$Epi.mm.High - apoc.wide$Epi.mm.Ambient
# Now plot ambient symbionts vs. delta for all variables
# symbiont density:
symbdens.vs.deltasymdens <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Symbs.cm2/100000))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Symbiont Density", paste((10^5~cells~cm^-2)))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 9, y = 10, label = "** p=0.006") +
  scale_x_continuous(breaks = c(2,4,6,8,10,12)) 
symbdens.vs.deltasymdens
symbdens.vs.deltasymdens.lm <- lm(delta.Symbs.cm2 ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.deltasymdens.lm)
# Corals with more symbionts tended to lose them (negative delta)
plot(symbdens.vs.deltasymdens.lm)

# epidermal thickness:
epithick.vs.deltasymdens <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.EpidermalThickness*1000))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Epidermal thickness", paste("(µm)"))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 9, y = 4.8, label = "* p=0.015") + 
  scale_x_continuous(breaks = c(2,4,6,8,10,12))
epithick.vs.deltasymdens
epithick.vs.deltasymdens.lm <- lm(delta.EpidermalThickness ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(epithick.vs.deltasymdens.lm)
plot(epithick.vs.deltasymdens.lm)

# protein:
symbdens.vs.prot <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Prot))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Host Protein", paste((µg~cm^-2)))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA))+
  annotate("text", x = 9, y = 60, label = "* p=0.014")+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))
symbdens.vs.prot
symbdens.vs.prot.lm <- lm(delta.Prot ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.prot.lm)
# Corals with more symbionts at the start lost more protein when heated

# total calcification:
symbdens.vs.calcif <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Calcif.g.cm2.day))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Symbiocyte", paste("Intracellular pH"))), 
       x = expression(atop("Symbiont Density", paste((10^5~cells~cm^-2))))) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.0155, label = "n.s.")+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))
symbdens.vs.calcif
symbdens.vs.calcif.lm <- lm(delta.Calcif.g.cm2.day ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.calcif.lm)
# no relationship

# symbiocyte pHi:
symbdens.vs.sympHi <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Symb.pHi))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Symbiocyte", paste("Intracellular pH"))), 
       x = expression(atop("Symbiont Density", paste((10^5~cells~cm^-2))))) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 1.55, label = "n.s.")+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))
symbdens.vs.sympHi
symbdens.vs.symphi.lm <- lm(delta.Symb.pHi ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.symphi.lm)
# no relationship

# nonsymbiocyte pHi:
symbdens.vs.nonsympHi <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Nonsymb.pHi))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Symbiocyte", paste("Intracellular pH"))), 
       x = expression(atop("Symbiont Density", paste((10^5~cells~cm^-2))))) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.65, label = "n.s.")+
  scale_x_continuous(breaks = c(2,4,6,8,10,12))
symbdens.vs.nonsympHi
symbdens.vs.nonsymphi.lm <- lm(delta.Nonsymb.pHi ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.nonsymphi.lm)
# no relationship
# All together: ambient genet symbiont density ~ ∆ of other variables
symbdens.vs.deltasymdens + epithick.vs.deltasymdens + symbdens.vs.prot + symbdens.vs.calcif + symbdens.vs.sympHi + symbdens.vs.nonsympHi

#### For Figure S2: ambient genet symbiont density ~ ∆ chlorophyll
symbdens.vs.deltachl <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.relchl))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Relative Chlorophyll", paste(Fluorescence))), 
       x = expression(paste(Symbiont~Density~(10^5~cells~cm^-2)))) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 9, y = -1000, label = "n.s.") +
  scale_x_continuous(breaks = c(2,4,6,8,10,12)) 
symbdens.vs.deltachl
symbdens.vs.deltachl.lm <- lm(delta.relchl ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.deltachl.lm)


###### pHi standard: plotting all Apoc cell measurements (for Figure S1) #####

# Read data
std <- read.csv("ApocStdAllMeasurements.csv")
std$curve <- as.factor(std$curve)

# Plot all raw datapoints colored by curve
# plus averages for each solution pH for each curve
RawRatios <- ggplot(data = std, aes(x = solution.pH, y = avg.yr.ratio)) +
  geom_point(alpha = 0.3, size = 2, aes(color=curve)) +
  stat_summary(geom = "point", fun = "mean", size = 3, color = "black", shape = 21, aes(fill = curve)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("#40B0A6", "#E1BE6A")) +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A")) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  coord_cartesian(xlim = c(6,8.55))
RawRatios

# Calculate log ratios
# Trim off ends and split by curve
logdf <- std %>% subset(solution.pH < 8.2 & solution.pH > 6.2)
# Calculate logs based on raw extremes
logdf$logged.ratio <- log10((logdf$avg.yr.ratio-0.52593)/(1.73127-logdf$avg.yr.ratio)*1.83614)
# Plot
LogRatios <- ggplot(data = logdf, aes(x = solution.pH, y = logged.ratio)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "solid")+
  geom_point(alpha = 0.3, size = 2, aes(color=curve)) +
  stat_summary(geom = "point", fun = "mean", size = 3, color = "black", shape = 21, aes(fill = curve)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_manual(values = c("#40B0A6", "#E1BE6A")) +
  scale_fill_manual(values = c("#40B0A6", "#E1BE6A")) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA)) +
  coord_cartesian(xlim = c(6,8.55))
LogRatios

# Both graphs
RawRatios + LogRatios
