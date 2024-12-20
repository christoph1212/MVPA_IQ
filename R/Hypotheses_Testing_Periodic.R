
#clear the workspace
rm(list=ls())

# Remove scientific notation
options(scipen=999)

# Install and load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, R.matlab)
source("./R/spearman_ci.R")

############################################
############# Data Preparation #############
############################################


# ----------------- Pre EO -----------------

# Get Participants
Subs_Pre_EO <- read.csv("Data/subs_pre_EO.txt", header = F)

# Create empty data frame
Pre_EO_Data <- data.frame(matrix(ncol = 6, nrow = nrow(Subs_Pre_EO)))

names(Pre_EO_Data) <- c('ID', 'mean_theta', 'mean_alpha', 
                        'fluid_score', 'crystallized_score', 'sleepiness_score')

Pre_EO_Data$ID <- Subs_Pre_EO[,]

# Get EEG Data and add to data frame
EEG_Pre_EO <- readMat("Data/Ready_for_DDTBOX/EEG_sorted_cond/Full_Sample/Periodic/Pre/EyesOpen/eeg_sorted_cond.mat")

EEG_Pre_EO <- EEG_Pre_EO$eeg.sorted.cond[[1]]

for (i in 1:nrow(Pre_EO_Data)) {
  
  # Theta
  Pre_EO_Data[i,2] <- mean(EEG_Pre_EO[[1]][8:16,,i])
  
  # Alpha
  Pre_EO_Data[i,3] <- mean(EEG_Pre_EO[[1]][17:26,,i])
  
}

# Get Behavioral Data and add to data frame
Fluid_Pre_EO <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Fluid/Pre/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat")

Pre_EO_Data$fluid_score <- data.frame(Fluid_Pre_EO[[1]][[1]])[,]

Cryst_Pre_EO <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Crystallized/Pre/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat")

Pre_EO_Data$crystallized_score <- data.frame(Cryst_Pre_EO[[1]][[1]])[,]

Sleep_Pre_EO <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Pre_Sleepiness/Pre/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat")

Pre_EO_Data$sleepiness_score <- data.frame(Sleep_Pre_EO[[1]][[1]])[,]


# ----------------- Pre EC -----------------

# Get Participants
Subs_Pre_EC <- read.csv("Data/subs_pre_EC.txt", header = F)

# Create empty data frame
Pre_EC_Data <- data.frame(matrix(ncol = 6, nrow = nrow(Subs_Pre_EC)))

names(Pre_EC_Data) <- c('ID', 'mean_theta', 'mean_alpha', 
                        'fluid_score', 'crystallized_score', 'sleepiness_score')

Pre_EC_Data$ID <- Subs_Pre_EC[,]

# Get EEG Data and add to data frame
EEG_Pre_EC <- readMat("Data/Ready_for_DDTBOX/EEG_sorted_cond/Full_Sample/Periodic/Pre/EyesClosed/eeg_sorted_cond.mat")

EEG_Pre_EC <- EEG_Pre_EC$eeg.sorted.cond[[1]]

for (i in 1:nrow(Pre_EC_Data)) {
  
  # Theta
  Pre_EC_Data[i,2] <- mean(EEG_Pre_EC[[1]][8:16,,i])
  
  # Alpha
  Pre_EC_Data[i,3] <- mean(EEG_Pre_EC[[1]][17:26,,i])
  
}

# Get Behavioral Data and add to data frame
Fluid_Pre_EC <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Fluid/Pre/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat")

Pre_EC_Data$fluid_score <- data.frame(Fluid_Pre_EC[[1]][[1]])[,]

Cryst_Pre_EC <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Crystallized/Pre/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat")

Pre_EC_Data$crystallized_score <- data.frame(Cryst_Pre_EC[[1]][[1]])[,]

Sleep_Pre_EC <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Pre_Sleepiness/Pre/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat")

Pre_EC_Data$sleepiness_score <- data.frame(Sleep_Pre_EC[[1]][[1]])[,]


# ----------------- Post EO -----------------

# Get Participants
Subs_Post_EO <- read.csv("Data/subs_post_EO.txt", header = F)

# Create empty data frame
Post_EO_Data <- data.frame(matrix(ncol = 6, nrow = nrow(Subs_Post_EO)))

names(Post_EO_Data) <- c('ID', 'mean_theta', 'mean_alpha', 
                         'fluid_score', 'crystallized_score', 'sleepiness_score')

Post_EO_Data$ID <- Subs_Post_EO[,]

# Get EEG Data and add to data frame
EEG_Post_EO <- readMat("Data/Ready_for_DDTBOX/EEG_sorted_cond/Full_Sample/Periodic/Post/EyesOpen/eeg_sorted_cond.mat")

EEG_Post_EO <- EEG_Post_EO$eeg.sorted.cond[[1]]

for (i in 1:nrow(Post_EO_Data)) {
  
  # Theta
  Post_EO_Data[i,2] <- mean(EEG_Post_EO[[1]][8:16,,i])
  
  # Alpha
  Post_EO_Data[i,3] <- mean(EEG_Post_EO[[1]][17:26,,i])
  
}

# Get Behavioral Data and add to data frame
Fluid_Post_EO <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Fluid/Post/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat")

Post_EO_Data$fluid_score <- data.frame(Fluid_Post_EO[[1]][[1]])[,]

Cryst_Post_EO <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Crystallized/Post/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat")

Post_EO_Data$crystallized_score <- data.frame(Cryst_Post_EO[[1]][[1]])[,]

Sleep_Post_EO <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Post_Sleepiness/Post/EyesOpen/eeg_sorted_cond_regress_sorted_cond.mat")

Post_EO_Data$sleepiness_score <- data.frame(Sleep_Post_EO[[1]][[1]])[,]


# ----------------- Post EC -----------------

# Get Participants
Subs_Post_EC <- read.csv("Data/subs_post_EC.txt", header = F)

# Create empty data frame
Post_EC_Data <- data.frame(matrix(ncol = 6, nrow = nrow(Subs_Post_EC)))

names(Post_EC_Data) <- c('ID', 'mean_theta', 'mean_alpha', 
                         'fluid_score', 'crystallized_score', 'sleepiness_score')

Post_EC_Data$ID <- Subs_Post_EC[,]

# Get EEG Data and add to data frame
EEG_Post_EC <- readMat("Data/Ready_for_DDTBOX/EEG_sorted_cond/Full_Sample/Periodic/Post/EyesClosed/eeg_sorted_cond.mat")

EEG_Post_EC <- EEG_Post_EC$eeg.sorted.cond[[1]]

for (i in 1:nrow(Post_EC_Data)) {
  
  # Theta
  Post_EC_Data[i,2] <- mean(EEG_Post_EC[[1]][8:16,,i])
  
  # Alpha
  Post_EC_Data[i,3] <- mean(EEG_Post_EC[[1]][17:26,,i])
  
}

# Get Behavioral Data and add to data frame
Fluid_Post_EC <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Fluid/Post/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat")

Post_EC_Data$fluid_score <- data.frame(Fluid_Post_EC[[1]][[1]])[,]

Cryst_Post_EC <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Crystallized/Post/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat")

Post_EC_Data$crystallized_score <- data.frame(Cryst_Post_EC[[1]][[1]])[,]

Sleep_Post_EC <- readMat("Data/Ready_for_DDTBOX/Behavioral_scores/Full_Sample/Post_Sleepiness/Post/EyesClosed/eeg_sorted_cond_regress_sorted_cond.mat")

Post_EC_Data$sleepiness_score <- data.frame(Sleep_Post_EC[[1]][[1]])[,]


Pre_Data <- inner_join(Pre_EO_Data, Pre_EC_Data, by = c("ID", "fluid_score", 'crystallized_score', "sleepiness_score"), suffix = c("_EO", "_EC"))
Post_Data <- inner_join(Post_EO_Data, Post_EC_Data, by = c("ID", "fluid_score", 'crystallized_score', "sleepiness_score"), suffix = c("_EO", "_EC"))
Data <- inner_join(Pre_Data, Post_Data, by = c("ID"), suffix = c("_Pre", "_Post"))
Data$fluid_score_Post <- NULL
Data$crystallized_score_Post <- NULL
names(Data)[names(Data) == 'fluid_score_Pre'] <- 'fluid_score'
names(Data)[names(Data) == 'crystallized_score_Pre'] <- 'crystallized_score'
Data$mean_alpha <- rowMeans(Data[, c("mean_alpha_EO_Pre", "mean_alpha_EC_Pre", "mean_alpha_EO_Post", "mean_alpha_EC_Post")], na.rm = TRUE)
Data$mean_theta <- rowMeans(Data[, c("mean_theta_EO_Pre", "mean_theta_EC_Pre", "mean_theta_EO_Post", "mean_theta_EC_Post")], na.rm = TRUE)
Data <- Data[, c(1,4,5,6,11,14,15,2,3,7,8,9,10,12,13)]

############################################
############### Hypothesis 2 ###############
############################################

corr_theta_fluid <- cor.test(Data$mean_theta,Data$fluid_score, method = "spearman", exact = F, alternative = "less")
spearman_ci(Data$mean_theta,Data$fluid_score,corr_theta_fluid)

corr_alpha_fluid <- cor.test(Data$mean_alpha,Data$fluid_score, method = "spearman", exact = F, alternative = "greater")
spearman_ci(Data$mean_alpha,Data$fluid_score,corr_alpha_fluid)

corr_theta_sleep_Pre <- cor.test(rowMeans(Pre_Data[, c("mean_theta_EO", "mean_theta_EC")]), Pre_Data$sleepiness_score, method = "spearman", exact = F, alternative = "greater")
spearman_ci(rowMeans(Pre_Data[, c("mean_theta_EO", "mean_theta_EC")]), Pre_Data$sleepiness_score,corr_theta_sleep_Pre)

corr_theta_sleep_Post <- cor.test(rowMeans(Post_Data[, c("mean_theta_EO", "mean_theta_EC")]), Post_Data$sleepiness_score, method = "spearman", exact = F, alternative = "greater")
spearman_ci(rowMeans(Post_Data[, c("mean_theta_EO", "mean_theta_EC")]), Post_Data$sleepiness_score, corr_theta_sleep_Post)

corr_alpha_sleep_Pre <- cor.test(rowMeans(Pre_Data[, c("mean_alpha_EO", "mean_alpha_EC")]), Pre_Data$sleepiness_score, method = "spearman", exact = F, alternative = "greater")
spearman_ci(rowMeans(Pre_Data[, c("mean_alpha_EO", "mean_alpha_EC")]), Pre_Data$sleepiness_score, corr_alpha_sleep_Pre)

corr_alpha_sleep_Post <- cor.test(rowMeans(Post_Data[, c("mean_alpha_EO", "mean_alpha_EC")]), Post_Data$sleepiness_score, method = "spearman", exact = F, alternative = "greater")
spearman_ci(rowMeans(Post_Data[, c("mean_alpha_EO", "mean_alpha_EC")]), Post_Data$sleepiness_score, corr_alpha_sleep_Post)

p.adjust(c(corr_theta_fluid$p.value, corr_alpha_fluid$p.value, corr_theta_sleep_Pre$p.value, corr_theta_sleep_Post$p.value, corr_alpha_sleep_Pre$p.value, corr_alpha_sleep_Post$p.value), method = "holm")
p.adjust(c(corr_theta_fluid$p.value, corr_alpha_fluid$p.value, corr_theta_sleep_Post$p.value, corr_alpha_sleep_Post$p.value), method = "holm")
