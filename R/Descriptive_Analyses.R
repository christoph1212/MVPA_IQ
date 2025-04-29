####################################
####### Descriptive Analysis ####### 
####################################

## Load relevant packages for further analysis
if (!require("pacman")) install.packages("pacman")
pacman::p_load(psych, apaTables, tidyverse, effsize, ggsignif, ggpubr, psychometric, 
               stringr, ggstance, ggh4x)

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

# Add Sociodemographic Data

socio_table = read.csv("Data/SocioDemographics.txt", header = T, na.strings = c("", "NA"))

## Convert columns to factors and rename levels
socio_table$Gender <- as.factor(socio_table$Gender)
socio_table$HormonalContraceptives <- as.factor(socio_table$HormonalContraceptives)
socio_table$Ethnicity <- as.factor(socio_table$Ethnicity)
socio_table$MaritalStatus <- as.factor(socio_table$MaritalStatus)
socio_table$Occupancy <- as.factor(socio_table$Occupancy)
socio_table$HighestDegree <- as.factor(socio_table$HighestDegree)
socio_table$monthlyBruttoIncome <- as.factor(socio_table$monthlyBruttoIncome)

levels(socio_table$Gender)[levels(socio_table$Gender) == 1] <- "Female"
levels(socio_table$Gender)[levels(socio_table$Gender) == 2] <- "Male"
levels(socio_table$HormonalContraceptives)[levels(socio_table$HormonalContraceptives) == 1] <- "Yes"
levels(socio_table$HormonalContraceptives)[levels(socio_table$HormonalContraceptives) == 2] <- "No"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 1] <- "European"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 2] <- "Arabic"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 3] <- "African"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 4] <- "Asian"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 5] <- "Hisp"
levels(socio_table$Ethnicity)[levels(socio_table$Ethnicity) == 6] <- "other"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 1] <- "Single"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 2] <- "Relationship"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 3] <- "legal Partnership"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 4] <- "Divorced"
levels(socio_table$MaritalStatus)[levels(socio_table$MaritalStatus) == 5] <- "Widowed"
levels(socio_table$Occupancy)[levels(socio_table$Occupancy) == 1] <- "Working"
levels(socio_table$Occupancy)[levels(socio_table$Occupancy) == 2] <- "Student"
levels(socio_table$Occupancy)[levels(socio_table$Occupancy) == 3] <- "unemployed"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 1] <- "Hauptschule"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 2] <- "Realschule"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 3] <- "Abitur"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 4] <- "Hochschulabschluss"
levels(socio_table$HighestDegree)[levels(socio_table$HighestDegree) == 5] <- "Promotion"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 1] <- "<500"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 2] <- "500-1000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 3] <- "1000-2000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 4] <- "2000-3000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 5] <- ">3000"
levels(socio_table$monthlyBruttoIncome)[levels(socio_table$monthlyBruttoIncome) == 6] <- "no info"

# Filter included participants
socio_fullSample <- socio_table %>% 
  filter(ID %in% behav_fullSample$ID)

# Calculate descriptive statistics of the sample demographics
# absolute Values
table(socio_fullSample$HighestDegree)
table(socio_fullSample$Occupancy)
table(socio_fullSample$Ethnicity)

# relative Values
prop.table(table(socio_fullSample$HighestDegree))
prop.table(table(socio_fullSample$Occupancy))
prop.table(table(socio_fullSample$Ethnicity))


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

# Plot the Data for gf, gc, Tiredness, Exhaustion, Sleepiness

# Prep for distribution Plots
# gf
data_long_gf <- behav_fullSample %>%
  dplyr::select(gf_score) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Score")

data_long_gf$Variable <- factor(data_long_gf$Variable, levels = c("gf_score"))

data_hist_gf <- ggplot(data_long_gf, aes(x = Score, fill = Variable)) +
  geom_histogram(color = "white", binwidth = 1) +
  geom_boxplot(aes(y = 0), width = 16, color = "black", alpha = 0.6, position = position_nudge(y = -12)) +
  facet_wrap(~Variable, scales = "free_x",
             labeller = as_labeller(c(gf_score = "Fluid Intelligence"))) +
  scale_fill_manual(values = c("gf_score" = "#FFC000")) +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(0, 20))
    )
  ) +
  labs(x = element_blank(), 
       y = element_blank()) +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_blank(), 
        legend.position = "none")

# gc
data_long_gc <- behav_fullSample %>%
  dplyr::select(gc_score) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Score")

data_long_gc$Variable <- factor(data_long_gc$Variable, levels = c("gc_score"))

data_hist_gc <- ggplot(data_long_gc, aes(x = Score, fill = Variable)) +
  geom_histogram(color = "white") +
  geom_boxplot(aes(y = 0), width = 13, color = "black", alpha = 0.6, position = position_nudge(y = -9)) +
  facet_wrap(~Variable, scales = "free_x",
             labeller = as_labeller(c(gc_score = "Crystallized Intelligence"))) +
  scale_fill_manual(values = c("gc_score" = "#4F81BD")) +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(0, 80))
    )
  ) +
  labs(x = element_blank(), 
       y = element_blank()) +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_blank(), 
        legend.position = "none")

# Pre-Task Sleepiness
data_long_s_pre <- behav_fullSample %>%
  dplyr::select(Pre_Task_Sleepiness) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Score")

data_long_s_pre$Variable <- factor(data_long_s_pre$Variable, levels = c("Pre_Task_Sleepiness"))

data_hist_s_pre <- ggplot(data_long_s_pre, aes(x = Score, fill = Variable)) +
  geom_histogram(color = "white", binwidth = 0.5) +
  geom_boxplot(aes(y = 0), width = 30, color = "black", alpha = 0.6, position = position_nudge(y = -21)) +
  facet_wrap(~Variable, scales = "free_x",
             labeller = as_labeller(c(
               Pre_Task_Sleepiness = "Pre-Task Sleepiness"
             ))) +
  scale_fill_manual(values = c("Pre_Task_Sleepiness" = "#D74E2D")) +
  labs(x = element_blank(), 
       y = element_blank()) +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_blank(), 
        legend.position = "none")

# Post-Task Sleepiness
data_long_s_post <- behav_fullSample %>%
  dplyr::select(Post_Task_Sleepiness) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Score")

data_long_s_post$Variable <- factor(data_long_s_post$Variable, levels = c("Post_Task_Sleepiness"))

data_hist_s_post <- ggplot(data_long_s_post, aes(x = Score, fill = Variable)) +
  geom_histogram(color = "white", binwidth = 0.5) +
  geom_boxplot(aes(y = 0), width = 22, color = "black", alpha = 0.6, position = position_nudge(y = -16)) +
  facet_wrap(~Variable, scales = "free_x",
             labeller = as_labeller(c(
               Post_Task_Sleepiness = "Post-Task Sleepiness"
             ))) +
  scale_fill_manual(values = c("Post_Task_Sleepiness" = "#D74E2D")) +
  labs(x = element_blank(), 
       y = element_blank()) +
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        legend.title = element_blank(), 
        legend.position = "none")

# Pack them together
library(patchwork)

hist_boxplt <- (data_hist_gf | data_hist_gc) / (data_hist_s_pre | data_hist_s_post)
ggsave("Histogram_Boxplt.jpeg", hist_boxplt, 
       dpi=600, width = 16, height = 10, path = "Results/Plots_FullSample")


# Prep for remaining Plots
behav_fullSample <- behav_fullSample %>%
  rename(Sleepiness_Pre = Pre_Task_Sleepiness,
         Sleepiness_Post = Post_Task_Sleepiness)

behav_long <- behav_fullSample %>%
  pivot_longer(cols = c(Tired_Pre, Tired_Post, Exhausted_Pre, Exhausted_Post, 
                        Sleepiness_Pre, Sleepiness_Post),
               names_to = c("Measure", "Time"),
               names_sep = "_",
               values_to = "Score")


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



logs <- read.csv("Data/Preprocessed/LogFiles/LogFiles_all.csv")

included_logs <- logs %>%
  filter(!str_detect(ID, "excluded"))

mean_ICA_EO <- mean(included_logs$ICA_Bad_Components_EO, na.rm = T)
mean_ICA_EC <- mean(included_logs$ICA_Bad_Components_EC, na.rm = T)

included_logs$Bad_Channel_Count <- sapply(strsplit(ifelse(is.na(included_logs$Bad_Channels), "", included_logs$Bad_Channels), ","), length)
mean_BadC <- mean(included_logs$Bad_Channel_Count)

mean_Epochs_EO <- mean(included_logs$Epochs_EO, na.rm = T)
mean_Epochs_EC <- mean(included_logs$Epochs_EC, na.rm = T)


# -----------------------------------
# -------- Preproc Log Data ---------
# -----------------------------------

mean_ICA_EO <- mean(logs$ICA_Bad_Components_EO, na.rm = T)
mean_ICA_EC <- mean(logs$ICA_Bad_Components_EC, na.rm = T)


mean_Epochs_EO <- mean(logs$Epochs_EO, na.rm = T)
mean_Epochs_EC <- mean(logs$Epochs_EC, na.rm = T)


logs$Bad_Channel_Count <- sapply(strsplit(ifelse(is.na(logs$Bad_Channels), "", logs$Bad_Channels), ","), length)
mean_BadC <- mean(logs$Bad_Channel_Count)

