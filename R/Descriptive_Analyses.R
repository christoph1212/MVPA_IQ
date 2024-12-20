####################################
####### Descriptive Analysis ####### 
####################################

## Load relevant packages for further analysis
if (!require("pacman")) install.packages("pacman")
pacman::p_load(psych, apaTables, tidyverse, effsize, ggsignif, ggpubr, psychometric, )

options(scipen=999)

# ----------------------------------
# ---------- Demographics ----------
# ----------------------------------

## Load Participant IDs for Pre- and Post-Subs
subs_pre_EO <- read.csv("Data/subs_pre_EO.txt", header = FALSE)
subs_pre_EO <- subs_pre_EO %>%
  rename(ID = V1)

subs_pre_EC <- read.csv("Data/subs_pre_EC.txt", header = FALSE)
subs_pre_EC <- subs_pre_EC %>%
  rename(ID = V1)

subs_post_EO <- read.csv("Data/subs_post_EO.txt", header = FALSE)
subs_post_EO <- subs_post_EO %>%
  rename(ID = V1)

subs_post_EC <- read.csv("Data/subs_post_EC.txt", header = FALSE)
subs_post_EC <- subs_post_EC %>%
  rename(ID = V1)

## Descriptive Analysis

# Add Behavioral Data

behav_data = read.csv("Data/BehavioralData/Behavioral_Data.csv", header = T, na.strings=c("","NA"))
behav_data[, 2:12][behav_data[,2:12] == "NaN"] <- NA

behav_data$Gender <- as.factor(behav_data$Gender)

# Create data frames for each condition
behav_pre_EO = merge(behav_data, subs_pre_EO, by='ID')
behav_pre_EC = merge(behav_data, subs_pre_EC, by='ID')
behav_pre = merge(behav_pre_EO, behav_pre_EC)
behav_post_EO = merge(behav_data, subs_post_EO, by='ID')
behav_post_EC = merge(behav_data, subs_post_EC, by='ID')
behav_post = merge(behav_post_EO, behav_post_EC)
behav_fullSample = full_join(behav_pre, behav_post)

describe(behav_fullSample$Age)

behav_fullSample %>%
  group_by(Gender) %>%
  summarise(
    Ratio = n()/nrow(.)*100,
    N = n()
  )

describe(behav_fullSample)

describeBy(behav_fullSample, group = behav_fullSample$Gender)

# Compare Pre and Post Data
t.test(behav_fullSample$Exhausted_Post, behav_fullSample$Exhausted_Pre, paired = TRUE)
effsize::cohen.d(behav_fullSample$Exhausted_Post, behav_fullSample$Exhausted_Pre, paired = TRUE)

t.test(behav_fullSample$Tired_Post, behav_fullSample$Tired_Pre, paired = TRUE)
effsize::cohen.d(behav_fullSample$Tired_Post, behav_fullSample$Tired_Pre, paired = TRUE)

t.test(behav_fullSample$Post_Task_Sleepiness, behav_fullSample$Pre_Task_Sleepiness, paired = TRUE)
effsize::cohen.d(behav_fullSample$Post_Task_Sleepiness, behav_fullSample$Pre_Task_Sleepiness, paired = TRUE)

# ----------------------------------
# -------- Behavioral Data ---------
# ----------------------------------

# Bivariate Scatter Plots, Histograms, and Distribution
pairs.panels(behav_fullSample[,4:12])

behav_fullSample %>%
  summarise(
    SpearmanBrown_Pre = (2 * cor(Tired_Pre, Exhausted_Pre)) / (1 + cor(Tired_Pre, Exhausted_Pre))
  )

behav_fullSample %>%
  summarise(
    SpearmanBrown_Post = (2 * cor(Tired_Post, Exhausted_Post)) / (1 + cor(Tired_Post, Exhausted_Post))
  )

behav_fullSample %>%
  summarise(
    SpearmanBrown_IQ = (2 * cor(gf_score, gc_score_z, use = "complete.obs")) / (1 + cor(gf_score, gc_score_z, use = "complete.obs"))
  )

# Principal Component Analysis with Parallel Analysis for Sleepiness
fa.parallel(behav_fullSample[,4:7], fa = "pc", main = "Parallel Analysis Sleepiness")
fa_sleepiness_results <- principal(behav_fullSample[,4:7], nfactors = 1, rotate = "varimax")
fa.diagram(fa_sleepiness_results)

# Principal Component Analysis with Parallel Analysis for Intelligence
fa.parallel(behav_fullSample[,8:9], fa = "pc", main = "Parallel Analysis Intelligence")
fa_intelligence_results <- principal(behav_fullSample[,8:9], nfactors = 1, rotate = "varimax")
fa.diagram(fa_intelligence_results)

# Calculate McDonald's Omega and Cronbach's Alpha 
omega(behav_fullSample[,4:7])

# APA-Style Correlation Table
# apa.cor.table(behav_fullSample[,4:12], filename = "Results/Correlation Table.doc", show.conf.interval = F)



IST_full <- read.csv("Data/BehavioralData/IST_fulltable.csv", header = T)

IST_full$fluid <- rowSums(IST_full[,2:21])

IST_full$cryst <- rowSums(IST_full[,22:105])

IST_filtered <- IST_full %>%
  filter(ID %in% behav_fullSample$ID)

omega(IST_filtered[,2:21])

selectivity <- item.total(IST_filtered[,2:21])$Item.Total

selectivity_corrected <- rangeCorrection(selectivity, sdu = 3.43, sdr = sd(rowSums(IST_filtered[,2:21])))

omega(IST_filtered[,22:105])

# Plot the Data for Tiredness, Exhaustion, Sleepiness
# Spaltennamen für Task_Sleepiness umbenennen, damit Pre/Post am Ende stehen
behav_fullSample <- behav_fullSample %>%
  rename(Sleepiness_Pre = Pre_Task_Sleepiness,
         Sleepiness_Post = Post_Task_Sleepiness)

# Jetzt ins lange Format umwandeln
behav_long <- behav_fullSample %>%
  pivot_longer(cols = c(Tired_Pre, Tired_Post, Exhausted_Pre, Exhausted_Post, 
                        Sleepiness_Pre, Sleepiness_Post),
               names_to = c("Measure", "Time"),
               names_sep = "_",
               values_to = "Score")

# Um die Zeitpunkte und Maße zu ordnen
behav_long$Time <- factor(behav_long$Time, levels = c("Pre", "Post"))
behav_long$Measure <- factor(behav_long$Measure, levels = c("Tired", "Exhausted", "Sleepiness"))

# Tiredness
tired_plot <- ggplot(behav_long %>% filter(Measure == "Tired"), aes(x = Time, y = Score)) +
  geom_violin(aes(fill = Time), alpha = 0.5) + 
  geom_jitter(width = 0.2, height = 0.05, alpha = 0.2) +
  labs(title = "Pre vs Post Tiredness", 
       x = element_blank(), 
       y = "Tiredness",
       tag = "A") +
  theme_classic() +
  theme(legend.title=element_blank()) +
  geom_signif(comparisons = list(c("Pre", "Post")), map_signif_level = T)

tired_plot

# Exhaustion
exhaustion_plot <- ggplot(behav_long %>% filter(Measure == "Exhausted"), aes(x = Time, y = Score)) +
  geom_violin(aes(fill = Time), alpha = 0.5) + 
  geom_jitter(width = 0.2, height = 0.05, alpha = 0.2) +
  labs(title = "Pre vs Post Exhaustion", 
       x = element_blank(), 
       y = "Exhaustion",
       tag = "B") +
  theme_classic() +
  theme(legend.title=element_blank()) +
  geom_signif(comparisons = list(c("Pre", "Post")), map_signif_level = T)

exhaustion_plot

# Sleepiness
sleepiness_plot <- ggplot(behav_long %>% filter(Measure == "Sleepiness"), aes(x = Time, y = Score)) +
  geom_violin(aes(fill = Time), alpha = 0.5) + 
  geom_jitter(width = 0.2, height = 0.05, alpha = 0.2) +
  labs(title = "Pre vs Post Sleepiness", 
       x = element_blank(), 
       y = "Sleepiness",
       tag = "C") +
  theme_classic() +
  theme(legend.title=element_blank()) +
  geom_signif(comparisons = list(c("Pre", "Post")), map_signif_level = T)

sleepiness_plot

ggarrange(tired_plot, exhaustion_plot, sleepiness_plot, nrow = 1)
