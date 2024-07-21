# Astrangia poculata heat response
# Luella Allen-Waller
# 2024-05-07

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
library(MuMIn)
library(devtools)
library(ggbiplot)
library(vegan)
library(emmeans)

setwd("~/Library/CloudStorage/Box-Box/grp-sas-bio-barottlab/Data/2023 Astrangia heat stress")

######## Temperature and visual observations

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

######## Physiology and pHi
# Load phys and pHi data and filter to only analyze colonies with matched genotypes in both treatments
apoc <- read.csv("ApocHeatStress.csv") %>% drop_na(Colony)
amb.apoc <- apoc %>% subset(Treatment == "Ambient")

# For Fig 1: 
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

######## Histology metrics (Fig 3)

# Symbiont quality
SymbQ <- ggplot(apoc, aes(x=Treatment, y=Symbiont.Health, color=Treatment)) +
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
SymbQ
# linear mixed effects model
SymbQ.lmer <- lmer(Symbiont.Health ~ Treatment + (1|Colony), data = apoc)
summary(SymbQ.lmer) # significant effect of treatment (*p=0.0362)

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
(SymbQ | TissueInteg) / (Necrosis | Epidermis) / (EpidermThick | (TissueInteg / Eggs.per.mm2))


###### Ambient symbiont density vs. delta for each variable
apoc.wide <- apoc %>% select(Treatment, Colony, Symbs.cm2, rel.chl, prot.ug.cm2, Symb.pHi, Nonsymb.pHi, Calcif.percent.change.day,
                             Symbiont.Health,	Epidermis.Integrity, Tissue.Integrity, Necrosis.Granularity, Overall.Health.Score, 
                             Gametes, Egg.Size, Egg.Count, Eggs.mm2, Epi.mm) %>%
  reshape(idvar = "Colony", timevar = "Treatment", direction = "wide")
apoc.wide$delta.Symbs.cm2 <- apoc.wide$Symbs.cm2.High - apoc.wide$Symbs.cm2.Ambient
apoc.wide$delta.relchl <- apoc.wide$rel.chl.High - apoc.wide$rel.chl.Ambient
apoc.wide$delta.Prot <- apoc.wide$prot.ug.cm2.High - apoc.wide$prot.ug.cm2.Ambient
apoc.wide$delta.Symb.pHi <- apoc.wide$Symb.pHi.High - apoc.wide$Symb.pHi.Ambient
apoc.wide$delta.Nonsymb.pHi <- apoc.wide$Nonsymb.pHi.High - apoc.wide$Nonsymb.pHi.Ambient
apoc.wide$delta.Calcif.percent.change.day <- apoc.wide$Calcif.percent.change.day.High - apoc.wide$Calcif.percent.change.day.Ambient
apoc.wide$delta.Symbperprot <- (apoc.wide$Symbs.cm2.High/apoc.wide$prot.ug.cm2.High) - (apoc.wide$Symbs.cm2.Ambient/apoc.wide$prot.ug.cm2.Ambient)
apoc.wide$delta.Symbiont.Health <- apoc.wide$Symbiont.Health.High - apoc.wide$Symbiont.Health.Ambient
apoc.wide$delta.Epidermis.Integrity <- apoc.wide$Epidermis.Integrity.High - apoc.wide$Epidermis.Integrity.Ambient
apoc.wide$delta.Tissue.Integrity <- apoc.wide$Tissue.Integrity.High - apoc.wide$Tissue.Integrity.Ambient
apoc.wide$delta.Necrosis.Granularity <- apoc.wide$Necrosis.Granularity.High - apoc.wide$Necrosis.Granularity.Ambient
apoc.wide$delta.Overall.Health.Score <- apoc.wide$Overall.Health.Score.High - apoc.wide$Overall.Health.Score.Ambient
apoc.wide$delta.Gametes <- apoc.wide$Gametes.High - apoc.wide$Gametes.Ambient
apoc.wide$delta.Egg.Size <- apoc.wide$Egg.Size.High - apoc.wide$Egg.Size.Ambient
apoc.wide$delta.Egg.Count <- apoc.wide$Egg.Count.High - apoc.wide$Egg.Count.Ambient
apoc.wide$delta.Eggs.mm2 <- apoc.wide$Eggs.mm2.High - apoc.wide$Eggs.mm2.Ambient
apoc.wide$delta.Epi.mm <- apoc.wide$Epi.mm.High - apoc.wide$Epi.mm.Ambient

# Plot ambient symbionts vs. delta for all variables

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
  annotate("text", x = 10, y = 10, label = "** p=0.006") +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.deltasymdens
symbdens.vs.deltasymdens.lm <- lm(delta.Symbs.cm2 ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.deltasymdens.lm)
# Corals with more symbionts tended to lose them (negative delta)
plot(symbdens.vs.deltasymdens.lm)

# Symbionts per protein:
symbdens.vs.deltasymperprot <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Symbperprot))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Symbiont Density", paste((cells~µg~host~protein^-1)))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 10, y = 2000, label = "** p=0.0059") +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.deltasymperprot
symbdens.vs.deltasymperprot.lm <- lm(delta.Symbperprot ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.deltasymperprot.lm)
# not significant
# did delta sym correlate with delta sym per prot?
deltasymbdens.vs.deltasymperprot.lm <- lm(delta.Symbperprot ~ delta.Symbs.cm2, data = apoc.wide)
summary(deltasymbdens.vs.deltasymperprot.lm)
# yes, still significantly correlated

# chlorophyll:
symbdens.vs.relchl <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.relchl))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Relative Chlorophyll", paste("Fluorescence"))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 10, y = -2000, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.relchl
symbdens.vs.relchl.lm <- lm(delta.relchl ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.relchl.lm)
# No relationship between symbiont density at the start and ∆ chlorophyll per symbiont

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
  annotate("text", x = 10, y = 60, label = "* p=0.014")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.prot
symbdens.vs.prot.lm <- lm(delta.Prot ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.prot.lm)
# Corals with more symbionts at the start lost more protein when heated

# total calcification:
symbdens.vs.calcif <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Calcif.percent.change.day))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Calcification", paste(("% change"~day^-1)))), 
       x = expression(atop("Symbiont Density", paste((10^5~cells~cm^-2))))) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.45, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.calcif
symbdens.vs.calcif.lm <- lm(delta.Calcif.percent.change.day ~ Symbs.cm2.Ambient, data = apoc.wide)
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
  annotate("text", x = 11, y = 1.4, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.sympHi
symbdens.vs.symphi.lm <- lm(delta.Symb.pHi ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.symphi.lm)
# no relationship

# nonsymbiocyte pHi:
symbdens.vs.nonsympHi <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Nonsymb.pHi))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Nonsymbiocyte", paste("Intracellular pH"))), 
       x = expression(atop("Symbiont Density", paste((10^5~cells~cm^-2))))) +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.6, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.nonsympHi
symbdens.vs.nonsymphi.lm <- lm(delta.Nonsymb.pHi ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.nonsymphi.lm)
# no relationship

# Symbiont integrity
symbdens.vs.deltasymbhealth <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Symbiont.Health))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Symbiont", paste("Integrity"))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.2, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.deltasymbhealth
symbdens.vs.deltasymbhealth.lm <- lm(delta.Symbiont.Health ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.deltasymbhealth.lm)
# no relationship

# Epidermis integrity
symbdens.vs.delta.epi.integ <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Epidermis.Integrity))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Epidermis", paste("Integrity"))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.2, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.delta.epi.integ
symbdens.vs.delta.epi.integ.lm <- lm(delta.Epidermis.Integrity ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.delta.epi.integ.lm) # no relationship

# Tissue integrity
symbdens.vs.delta.tissue.integ <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Tissue.Integrity))+
geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Tissue", paste("Integrity"))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.2, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.delta.tissue.integ
symbdens.vs.delta.tissue.integ.lm <- lm(delta.Tissue.Integrity ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.delta.tissue.integ.lm) # no relationship

# Necrosis
symbdens.vs.delta.necrosis <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Necrosis.Granularity))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Host", paste("Necrosis"))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.2, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.delta.necrosis
symbdens.vs.delta.necrosis.lm <- lm(delta.Necrosis.Granularity ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.delta.necrosis.lm) # no relationship

# Overall Health Score
symbdens.vs.delta.health <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Overall.Health.Score))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Overall Host", paste("Health Score"))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.2, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.delta.health
symbdens.vs.delta.health.lm <- lm(delta.Overall.Health.Score ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.delta.health.lm) # no relationship

# Eggs per area
symbdens.vs.delta.eggcount <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Eggs.mm2))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Host Egg Count", paste((Eggs~mm^-2)))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 11, y = 0.2, label = "n.s.")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.delta.eggcount
symbdens.vs.delta.eggcount.lm <- lm(delta.Eggs.mm2 ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.delta.eggcount.lm) # no relationship

# Epidermal thickness
symbdens.vs.delta.epithick <- ggplot(apoc.wide, aes(x=Symbs.cm2.Ambient/100000, y=delta.Epi.mm*1000))+
  geom_smooth(method = "lm", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(color="black", alpha = 0.6) +
  labs(y = expression(atop("∆ Host Epidermal", paste(Thickness~(µm)))), 
       x = "") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none", legend.key = element_rect(fill = NA)) +
  annotate("text", x = 10, y = 5, label = "* p=0.015")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12))
symbdens.vs.delta.epithick
symbdens.vs.delta.epithick.lm <- lm(delta.Epi.mm ~ Symbs.cm2.Ambient, data = apoc.wide)
summary(symbdens.vs.delta.epithick.lm) # negative relationship (p=0.015)

# For Figure 4: ∆ symbiont, biomass, and pH-dependent data regressed by initial symbiont density:
symbdens.vs.deltasymdens + symbdens.vs.delta.epithick + symbdens.vs.prot + symbdens.vs.calcif + symbdens.vs.sympHi + symbdens.vs.nonsympHi


####### For Figure 5A (Symbiont density responses)
# How did temperature affect symb density?
SymbsByTreatment <- ggplot(apoc, aes(x=Treatment, y=Symbs.cm2/100000, color = Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group=Colony), color = "black", alpha = 0.2) +
  geom_point(alpha = 0.5, size = 2, position = position_dodge2(width = 0.04)) +
  labs(y = expression(atop("Symbiont density", paste((10^5~cells~cm^-2)))), x = "Temperature") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582")) +
  annotate("text", x = 1.5, y = 13, label = "n.s.")
SymbsByTreatment
# linear mixed model
symbs.lmer <- lmer(Symbs.cm2 ~ Treatment + (1|Colony), data = apoc)
summary(symbs.lmer) # no significant effects


######### For Figure 6: Box plots faceted by whether symbionts were lost or gained
# Generate Symbs.Response
apoc.wide$Symbs.Response <- ifelse(apoc.wide$delta.Symbs.cm2<0, "Lost",
                                   ifelse(apoc.wide$delta.Symbs.cm2>0, "Gained", 
                                          "neutral"))
symbsresponse <- apoc.wide %>%
  select(Colony, Symbs.Response)
apoc <- left_join(apoc, symbsresponse, by = "Colony")
apoc.gained <- apoc %>% subset(Symbs.Response == "Gained")
apoc.lost <- apoc %>% subset(Symbs.Response == "Lost")

# Chlorophyll:
ChlByTreatment <- ggplot(apoc, aes(x=Treatment, y=rel.chl, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Relative Chlorophyll", paste("Fluorescence"))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582")) +
  scale_shape_manual(values = c(19,1))
ChlByTreatment
# stats
chl.gain.lm <- lm(rel.chl ~ Treatment, data = apoc.gained)
summary(chl.gain.lm)
# significant heat detriment to chl (***p<0.001) for corals that gained symbionts
chl.loss.lm <- lm(rel.chl ~ Treatment, data = apoc.lost)
summary(chl.loss.lm)
# significant heat detriment to chl (***p<0.001) for corals that lost symbionts
# so, heat = lower chl per symb regardless of whether symbs were lost or gained

# Protein:
ProtByTreatment <- ggplot(apoc, aes(x=Treatment, y=prot.ug.cm2, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Host Protein", paste((µg~cm^-2)))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582")) +
  scale_shape_manual(values = c(19,1))
ProtByTreatment
# stats
prot.gain.lm <- lm(prot.ug.cm2 ~ Treatment, data = apoc.gained)
summary(prot.gain.lm)
# no effect of heat (p=0.2) if symbionts were gained
prot.loss.lm <- lm(prot.ug.cm2 ~ Treatment, data = apoc.lost)
summary(prot.loss.lm)
# if corals lost symbionts, they also lost protein (**p=0.0033)

# Total Calcification:
CalcifByTreatment <- ggplot(apoc, aes(x=Treatment, y=Calcif.percent.change.day, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Total Calcification", paste(("%"~increase~d^-1)))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  geom_hline(yintercept=0, color = "darkgrey", linetype="dashed")+
  scale_shape_manual(values = c(19,1))
CalcifByTreatment
# stats
calcif.gain.lm <- lm(Calcif.percent.change.day ~ Treatment, data = apoc.gained)
summary(calcif.gain.lm)
# significant heat benefit to calcification (**p=0.0065) for corals that gained symbionts
calcif.loss.lm <- lm(Calcif.percent.change.day ~ Treatment, data = apoc.lost)
summary(calcif.loss.lm)
# no significant heat benefit to calcification

# Symb pHi:
SympHiByTreatment <- ggplot(apoc, aes(x=Treatment, y=Symb.pHi, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Symbiocyte", paste("Intracellular pH"))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1))
SympHiByTreatment
# stats
sympHi.gain.lm <- lm(Symb.pHi ~ Treatment, data = apoc.gained)
summary(sympHi.gain.lm)
# no effect of heat (p=0.3)
sympHi.loss.lm <- lm(Symb.pHi ~ Treatment, data = apoc.lost)
summary(sympHi.loss.lm)
# no effect of heat (p=0.99)

# Nonsymb pHi:
NonsympHiByTreatment <- ggplot(apoc, aes(x=Treatment, y=Nonsymb.pHi, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Nonsymbiocyte", paste("Intracellular pH"))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1))
NonsympHiByTreatment
# stats
nonsympHi.gain.lm <- lm(Nonsymb.pHi ~ Treatment, data = apoc.gained)
summary(nonsympHi.gain.lm)
# no effect of heat (p=0.24)
nonsympHi.loss.lm <- lm(Nonsymb.pHi ~ Treatment, data = apoc.lost)
summary(nonsympHi.loss.lm)
# no effect of heat (p=0.98)

# Symbiont integrity:
SymbIntegByTreatment <- ggplot(apoc, aes(x=Treatment, y=Symbiont.Health, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Symbiont", paste("Integrity"))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1))
SymbIntegByTreatment
# stats
symb.integ.gain.lm <- lm(Symbiont.Health ~ Treatment, data = apoc.gained)
summary(symb.integ.gain.lm)
# no effect of heat (p=0.23)
symb.integ.loss.lm <- lm(Symbiont.Health ~ Treatment, data = apoc.lost)
summary(symb.integ.loss.lm) # no effect of heat (p=0.12)

# Host epidermal clarity:
EpidermalClarityByTreatment <- ggplot(apoc, aes(x=Treatment, y=Epidermis.Integrity, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Host Epidermal", paste("Clarity"))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1))
EpidermalClarityByTreatment
# stats
EpidermalClarity.gain.lm <- lm(Epidermis.Integrity ~ Treatment, data = apoc.gained)
summary(EpidermalClarity.gain.lm)
# no effect of heat (p=0.56)
EpidermalClarity.loss.lm <- lm(Epidermis.Integrity ~ Treatment, data = apoc.lost)
summary(EpidermalClarity.loss.lm) # no effect of heat (p=0.83)

# Host tissue integrity:
TissueIntegrityByTreatment <- ggplot(apoc, aes(x=Treatment, y=Tissue.Integrity, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Host Tissue", paste("Integrity"))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1))
TissueIntegrityByTreatment
# stats
TissueIntegrity.gain.lm <- lm(Tissue.Integrity ~ Treatment, data = apoc.gained)
summary(TissueIntegrity.gain.lm)
# no effect of heat (p=0.056)
Tissue.Integrity.loss.lm <- lm(Tissue.Integrity ~ Treatment, data = apoc.lost)
summary(Tissue.Integrity.loss.lm) # no effect of heat (p=0.27)

# Necrosis
NecrosisByTreatment <- ggplot(apoc, aes(x=Treatment, y=Necrosis.Granularity, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Host Necrosis", paste("Score"))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1))
NecrosisByTreatment
# stats
NecrosisByTreatment.gain.lm <- lm(Necrosis.Granularity ~ Treatment, data = apoc.gained)
summary(NecrosisByTreatment.gain.lm)
# no effect of heat (p=0.39)
NecrosisByTreatment.loss.lm <- lm(Necrosis.Granularity ~ Treatment, data = apoc.lost)
summary(NecrosisByTreatment.loss.lm) # no effect of heat (p=0.82)

# Gametes:
GametesByTreatment <- ggplot(apoc, aes(x=Treatment, y=Gametes, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Egg Count", paste((eggs~mm^-2)))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1))
GametesByTreatment
# stats
# gained symbionts:
Gametes.gain.lm <- lm(Gametes ~ Treatment, data = apoc.gained)
summary(Gametes.gain.lm)
# no effect of heat (p=0.68)
Gametes.loss.lm <- lm(Gametes ~ Treatment, data = apoc.lost)
summary(Gametes.loss.lm)
# not quite an effect of heat (p=0.067)

# Eggs per mm2:
EggsByTreatment <- ggplot(apoc, aes(x=Treatment, y=Eggs.mm2, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Egg Count", paste((eggs~mm^-2)))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1))
EggsByTreatment
# stats
Eggs.mm2.gain.lm <- lm(Eggs.mm2 ~ Treatment, data = apoc.gained)
summary(Eggs.mm2.gain.lm)
# no effect of heat (p=0.78)
Eggs.mm2.loss.lm <- lm(Eggs.mm2 ~ Treatment, data = apoc.lost)
summary(Eggs.mm2.loss.lm)
# not quite an effect of heat (p=0.097)

# Epidermal thickness:
EpiThickByTreatment <- ggplot(apoc, aes(x=Treatment, y=Epi.mm, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(facets="Symbs.Response")+
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05), aes(shape=Symbs.Response)) +
  labs(y = expression(atop("Epidermal Thickness", paste((mm)))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1))
EpiThickByTreatment
# stats
Epi.mm.gain.lm <- lm(Epi.mm ~ Treatment, data = apoc.gained)
summary(Epi.mm.gain.lm)
# not quite an effect of heat (p=0.069)
Epi.mm.loss.lm <- lm(Epi.mm ~ Treatment, data = apoc.lost)
summary(Epi.mm.loss.lm)
# no effect of heat (p=0.33)

# All together:
ggarrange(ChlByTreatment, EpiThickByTreatment, ProtByTreatment, 
          CalcifByTreatment, SympHiByTreatment, NonsympHiByTreatment,
          nrow = 2, ncol = 3)

# Epidermal thickness (simple)
EpiThickByTreatmentSimple <- ggplot(apoc, aes(x=Treatment, y=Epi.mm, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05)) +
  labs(y = expression(atop("Epidermal Thickness", paste((mm)))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582"))+
  scale_shape_manual(values = c(19,1)) +
  annotate("text", x = 1.5, y = 0.015, label = "n.s.")
EpiThickByTreatmentSimple
# stats
Epi.mm.lmer <- lmer(Epi.mm ~ Treatment + (1|Colony), data = apoc)
summary(Epi.mm.lmer) # not quite an effect of temperature (p=0.0866)

###### For supplemental symbiont distress figure
ChlByTreatmentSimple <- ggplot(apoc, aes(x=Treatment, y=rel.chl, color=Treatment)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_line(aes(group = Colony), color = "black", alpha=0.2) +
  geom_point(alpha = 0.5, position = position_dodge2(width = 0.05)) +
  labs(y = expression(atop("Chlorophyll (Relative", paste(Fluorescence~Cell^-1~")"))), x = "") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "none") +
  scale_color_manual(values = c("#0571B0","#F4A582")) +
  scale_shape_manual(values = c(19,1)) +
  annotate("text", x = 1.5, y = 8500, label = "p<0.001") +
  coord_cartesian(ylim = c(0,9000))
ChlByTreatmentSimple
# Stats
Relchl.lmer <- lmer(rel.chl ~ Treatment + (1|Colony), data = apoc)
summary(Relchl.lmer) # significant decrease (p<0.001)
# Figure S2
ChlByTreatmentSimple + symbdens.vs.relchl
