# This script was modified from the original script from Jach et al. (2020).
# Plots for the periodic and aperiodic signal components and further plots for
# the supplements were added.
#
# ----------------------------------------------------------------------
# ----------------        MVPA Output Graphs      ----------------------
# ----------------------------------------------------------------------

# Script written by Hayley Jach

# Modified by Christoph Fruehlinger, 2024

# ---------- Housekeeping --------------
#clear the workspace
rm(list=ls())

# Remove scientific notation
options(scipen=999)

# Install necessary packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, reshape2, cowplot)

# Load necessary libraries
library(tidyverse) # contains tidyr, ggplot, dplyr etc.
library(reshape2) # for melting data etc.
library(cowplot) #for combining ggplots

# --------------------------------------------------
# ---------- Get individual 4-plots ---------------   
# --------------------------------------------------

myFolders <- list.files(path = "Results", pattern = "MVPA*")
foldersplits <- strsplit(myFolders, '_')

for (i_folders in 1:length(myFolders)) {
  
  print(paste("Preparing Plots for the", foldersplits[[i_folders]][3], "Sample."))
  
  # Create Output Folder
  if (!dir.exists(paste0("Results/Plots_", foldersplits[[i_folders]][3]))) {
    dir.create(paste0("Results/Plots_", foldersplits[[i_folders]][3]), recursive = TRUE)
  }
  
  # list files in the data 10 folds graphing directory
  myFiles <- list.files(path = paste0("Results/",myFolders[i_folders]), pattern=".*csv")
  
  # Exclude Parameter files
  myFiles <- myFiles[-grep("_Parameters_", myFiles)]
  
  # Make empty list to add data frames
  data_list <- list()
  
  # make empty list to add plots
  myplots <- list() 
  
  # Remove suffix
  myFilesNoCSV <- str_remove(myFiles, ".mat_averageCV.csv") 
  
  filesplits <- strsplit(myFilesNoCSV, '_')
  
  # Make list of file names
  FileNames <- character(length(myFilesNoCSV))
  PlotTitles <- character(length(myFilesNoCSV))
  Behav <- character(length(myFilesNoCSV))
  Condition <- character(length(myFilesNoCSV))
  Sample <- character(length(myFilesNoCSV))
  
  for (i in seq_along(myFilesNoCSV)) {
    
    if (i_folders == 2) {
      
      if (filesplits[[i]][6] == "Sleepiness") {
        
        filesplits[[i]][8] <- ifelse(filesplits[[i]][8] == "EyesOpen", "EO", filesplits[[i]][8])
        filesplits[[i]][8] <- ifelse(filesplits[[i]][8] == "EyesClosed", "EC", filesplits[[i]][8])
        
        FileNames[i] <- paste(filesplits[[i]][4], filesplits[[i]][6], filesplits[[i]][7], filesplits[[i]][8], filesplits[[i]][2], collapse = ' ')
        PlotTitles[i] <- paste(filesplits[[i]][6], filesplits[[i]][7], filesplits[[i]][8], collapse = ' ')
        Behav[i] <- paste(filesplits[[i]][6])
        Condition[i] <- paste(filesplits[[i]][7], filesplits[[i]][8], collapse = ' ')
        Sample[i] <- paste(filesplits[[i]][2])
        
      } else {
        
        filesplits[[i]][7] <- ifelse(filesplits[[i]][7] == "EyesOpen", "EO", filesplits[[i]][7])
        filesplits[[i]][7] <- ifelse(filesplits[[i]][7] == "EyesClosed", "EC", filesplits[[i]][7])
        
        FileNames[i] <- paste(filesplits[[i]][4], filesplits[[i]][5], filesplits[[i]][6], filesplits[[i]][7], filesplits[[i]][2], collapse = ' ')
        PlotTitles[i] <- paste(filesplits[[i]][5], filesplits[[i]][6], filesplits[[i]][7], collapse = ' ')
        Behav[i] <- paste(filesplits[[i]][5])
        Condition[i] <- paste(filesplits[[i]][6], filesplits[[i]][7], collapse = ' ')
        Sample[i] <- paste(filesplits[[i]][2])
        
      }
      
    } else {
      
      filesplits[[i]][6] <- ifelse(filesplits[[i]][6] == "EyesOpen", "EO", filesplits[[i]][6])
      filesplits[[i]][6] <- ifelse(filesplits[[i]][6] == "EyesClosed", "EC", filesplits[[i]][6])
      
      FileNames[i] <- paste(filesplits[[i]][3], filesplits[[i]][4], filesplits[[i]][5], filesplits[[i]][6], filesplits[[i]][2], collapse = ' ')
      PlotTitles[i] <- paste(filesplits[[i]][3], filesplits[[i]][4], filesplits[[i]][5], filesplits[[i]][6], collapse = ' ')
      Behav[i] <- paste(filesplits[[i]][4])
      Condition[i] <- paste(filesplits[[i]][5], filesplits[[i]][6], collapse = ' ')
      Sample[i] <- paste(filesplits[[i]][2])
      
    }
    
  }
  
  
  # Make list of color values for plots
  if (i_folders == 2) {
    color_vals <- c("#4F81BD","#4F81BD","#4F81BD","#4F81BD",
                    "#FFC000","#FFC000","#FFC000","#FFC000",
                    "#D74E2D","#D74E2D","#D74E2D","#D74E2D",
                    "#4F81BD","#4F81BD","#4F81BD","#4F81BD",
                    "#FFC000","#FFC000","#FFC000","#FFC000",
                    "#D74E2D","#D74E2D","#D74E2D","#D74E2D",
                    "#4F81BD","#4F81BD","#4F81BD","#4F81BD",
                    "#FFC000","#FFC000","#FFC000","#FFC000",
                    "#D74E2D","#D74E2D","#D74E2D","#D74E2D")
  } else {
    color_vals <- c("#4F81BD","#4F81BD","#4F81BD","#4F81BD",
                    "#FFC000","#FFC000","#FFC000","#FFC000",
                    "#4F81BD","#4F81BD","#4F81BD","#4F81BD",
                    "#FFC000","#FFC000","#FFC000","#FFC000",
                    "#4F81BD","#4F81BD","#4F81BD","#4F81BD",
                    "#FFC000","#FFC000","#FFC000","#FFC000")
  }
  
  # -------- Make graphs for all Behavioral variables, all conditions --------
  
  # Loop across all 
  # empty data frame for correlation heatmap (supplementary material)
  hm <- data.frame(matrix(ncol = length(myFiles), nrow = 30))
  
  for (i in 1:length(myFiles)) { # for i in 1 to the number of files
    
    # read in data
    data <- read.csv( 
      paste0("Results/", myFolders[i_folders],"/", myFiles[i]), header = FALSE)
    
    # Add column names
    names(data) <- c('1','2','3','4','5','6','7','8','9','10','P1','P2','P3','P4',
                     'P5','P6','P7','P8','P9','P10') # P = "permuted"
    
    # Create mean scores for variables
    data$Mean<- rowMeans(data[,c(1:10)]) 
    # Create mean permuted scores
    data$PermMean <- rowMeans(data[,c(11:20)]) 
    # Standard deviation 
    data$SD <- apply(data[,c(1:10)], 1, sd) 
    # Standard error
    data$SE <- data$SD/sqrt(10)
    # Standard deviation, permuted labels
    data$Perm_SD <- apply(data[,c(11:20)], 1, sd) 
    # Standard error, permuted labels
    data$Perm_SE <- data$Perm_SD/sqrt(10) 
    
    # Make dataframes with the means/SEs
    mean_data <- data[,c("Mean", "PermMean")]
    se_data <- data[,c("SE", "Perm_SE")]
    
    #Add ID column
    mean_data <- tibble::rowid_to_column(mean_data, "ID") 
    se_data <- tibble::rowid_to_column(se_data, "ID") 
    
    # pull everything into longform
    mean_melt <- melt(mean_data, id="ID") 
    se_melt <- melt(se_data, id="ID") 
    melt <- mean_melt
    melt$SE <- se_melt$value # "melt" is the main DF from hereon in
    melt$upper.se = melt$value + melt$SE # error bars, upper
    melt$lower.se = melt$value - melt$SE # error bars, lower
    melt$upper.percentile = melt$value + 1.96 * melt$SE
    melt$lower.percentile = melt$value - 1.96 * melt$SE
    melt <- melt[1:30,]
    
    data_list[[i]] <- melt
    
    hm[, i] <- melt$value
    colnames(hm)[i] <- FileNames[[i]]
    
    # Make GG Plot
    tempplot <- ggplot() +
      geom_line(
        data = melt, aes(x = ID, y = value, color = variable), linewidth = 0.75) +
      geom_ribbon(data = melt, aes(x=ID, ymin = lower.percentile, ymax = upper.percentile, 
                                   color = variable, fill = variable), 
                  linetype = "blank", alpha=.40) + 
      theme(legend.title = element_blank(),
            legend.position="none") +
      scale_y_continuous(breaks=seq(-.20,25,.20), limits = c(-.20,.25)) +
      labs(title = toString(PlotTitles[i]), y = "Correlation", 
           x = "Frequency (Hz)") +
      geom_hline(yintercept = .20, linetype = "dashed", color = "#8c8c8c") + 
      scale_color_manual(name="Decoding Accuracy", 
                         labels = c("M", "P M"),
                         values = c(toString(color_vals[i]), "grey50")) + 
      scale_fill_manual(name = "Decoding Accuracy",
                        labels = c("M", "P M"),
                        values = c(toString(color_vals[i]), "grey50"))
    
    # Add into list
    myplots[[i]] <- tempplot
    # Add names to list
    names(myplots[i]) <- toString(myFiles[i])
    
    
  } # end for i
  
  # For Supplement
  hm <- cbind(Frequency = 1:nrow(hm), hm)
  
  # Create new labels for axis description
  new_labels <- sapply(names(hm), function(name) {
    paste(unlist(strsplit(name, " "))[2:4], collapse = " ")
  })
  
  # Safe hm data frame for each sample
  if (i_folders == 1) {
    
    decoding_performance_female <- hm
    
  } else if (i_folders == 2) {
    
    decoding_performance_full <- hm
    
  } else {
  
    decoding_performance_male <- hm
      
  }
  
  
  # Combine all data frames into a single data frame
  combined_supp <- bind_rows(data_list, .id = "FileID")
  
  combined_supp <- combined_supp %>%
    mutate(FileNames = FileNames[as.numeric(FileID)], 
           PlotTitles = PlotTitles[as.numeric(FileID)],
           Behav = Behav[as.numeric(FileID)],
           Condition = Condition[as.numeric(FileID)])
  
  # Full Sample (i_folders = 2) also contains Sleepiness data
  if (i_folders == 2) {
    
    range_aper = 1:13
    range_per = 14:25
    range_total = 26:37
    
    behav_order <- c("Fluid", "Crystallized", "Sleepiness")
    combined_supp_aper <- combined_supp[1:360,]
    combined_supp_per <- combined_supp[361:720,]
    combined_supp_total <- combined_supp[721:1080,]
    
  } else {
    
    range_aper = 1:9
    range_per = 10:17
    range_total = 18:25
    
    behav_order <- c("Fluid", "Crystallized")
    combined_supp_aper <- combined_supp[1:240,]
    combined_supp_per <- combined_supp[241:480,]
    combined_supp_total <- combined_supp[481:720,]
    
  }
  
  condition_order <- c("Pre EO", "Pre EC", "Post EO", "Post EC")
  
  fileid_colors <- setNames(color_vals, unique(combined_supp$FileID))
  
  # Create Heatmaps for Supplement
  heatmap_aper <- ggplot(melt(hm[,range_aper], id.vars = "Frequency"), aes(x = variable, y = Frequency, fill = value)) + 
    geom_tile(colour="white", linewidth=0.25) + 
    scale_fill_viridis_c(name = "Correlation", guide = guide_colorbar(title.theme = element_text(size = 14)))  +
    scale_y_continuous(breaks = c(1, 10, 20, 30), expand = c(0, 0), limits = c(0.5, 30.5)) +
    scale_x_discrete(labels = new_labels) +
    theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    labs(y = "Frequency (Hz)")
  
  ggsave("Supp_Aperiodic_Heatmap.jpeg", heatmap_aper, 
         dpi=300, width = 10, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3]))
  
  heatmap_per <- ggplot(melt(hm[,c(1,range_per)], id.vars = "Frequency"), aes(x = variable, y = Frequency, fill = value)) + 
    geom_tile(colour="white", linewidth=0.25) + 
    scale_fill_viridis_c(name = "Correlation", guide = guide_colorbar(title.theme = element_text(size = 14))) +
    scale_y_continuous(breaks = c(1, 10, 20, 30), expand = c(0, 0), limits = c(0.5, 30.5)) +
    scale_x_discrete(labels = new_labels) +
    theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    labs(y = "Frequency (Hz)")
  
  ggsave("Supp_Periodic_Heatmap.jpeg", heatmap_per, 
         dpi=300, width = 10, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3]))
  
  heatmap_total <- ggplot(melt(hm[,c(1,range_total)], id.vars = "Frequency"), aes(x = variable, y = Frequency, fill = value)) + 
    geom_tile(colour="white", linewidth=0.25) + 
    scale_fill_viridis_c(name = "Correlation", guide = guide_colorbar(title.theme = element_text(size = 14))) +
    scale_y_continuous(breaks = c(1, 10, 20, 30), expand = c(0, 0), limits = c(0.5, 30.5)) +
    scale_x_discrete(labels = new_labels) +
    theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    labs(y = "Frequency (Hz)")
  
  ggsave("Supp_Total_Heatmap.jpeg", heatmap_total, 
         dpi=300, width = 10, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3]))
  
  
  # Plots for each Behavioral variable and condition
  grid_aper <- ggplot(combined_supp_aper, aes(x = ID, y = value, color = factor(FileID))) +
    geom_line(linewidth = 0.75) +
    geom_ribbon(aes(ymin = lower.percentile, ymax = upper.percentile, fill = factor(FileID)),
                linetype = "blank", alpha = 0.40) +
    scale_y_continuous(breaks = seq(-0.20, 0.25, 0.20), limits = c(-0.20, 0.25)) +
    geom_hline(yintercept = 0.20, linetype = "dashed", color = "#8c8c8c") +
    scale_color_manual(name = "FileID", values = fileid_colors) +
    scale_fill_manual(name = "FileID", values = fileid_colors) +
    facet_grid(factor(Condition, condition_order) ~ factor(Behav, behav_order)) +
    labs(x = "Frequency (Hz)", y = "Correlation") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 20),
          strip.text = element_text(size = 20),
          legend.title = element_blank(), 
          legend.position = "none")
  
  ggsave("Supp_Grid_aper.jpeg", grid_aper, 
         dpi=300, width = 16, height = 10, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
  
  
  grid_per <- ggplot(combined_supp_per, aes(x = ID, y = value, color = factor(FileID))) +
    geom_line(linewidth = 0.75) +
    geom_ribbon(aes(ymin = lower.percentile, ymax = upper.percentile, fill = factor(FileID)),
                linetype = "blank", alpha = 0.40) +
    scale_y_continuous(breaks = seq(-0.20, 0.25, 0.20), limits = c(-0.20, 0.25)) +
    geom_hline(yintercept = 0.20, linetype = "dashed", color = "#8c8c8c") +
    scale_color_manual(name = "FileID", values = fileid_colors) +
    scale_fill_manual(name = "FileID", values = fileid_colors) +
    facet_grid(factor(Condition, condition_order) ~ factor(Behav, behav_order)) +
    labs(x = "Frequency (Hz)", y = "Correlation") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 20),
          strip.text = element_text(size = 20),
          legend.title = element_blank(), 
          legend.position = "none")
  
  ggsave("Supp_Grid_per.jpeg", grid_per, 
         dpi=300, width = 16, height = 10, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
  
  
  grid_total <- ggplot(combined_supp_total, aes(x = ID, y = value, color = factor(FileID))) +
    geom_line(linewidth = 0.75) +
    geom_ribbon(aes(ymin = lower.percentile, ymax = upper.percentile, fill = factor(FileID)),
                linetype = "blank", alpha = 0.40) +
    scale_y_continuous(breaks = seq(-0.20, 0.25, 0.20), limits = c(-0.20, 0.25)) +
    geom_hline(yintercept = 0.20, linetype = "dashed", color = "#8c8c8c") +
    scale_color_manual(name = "FileID", values = fileid_colors) +
    scale_fill_manual(name = "FileID", values = fileid_colors) +
    facet_grid(factor(Condition, condition_order) ~ factor(Behav, behav_order)) +
    labs(x = "Frequency (Hz)", y = "Correlation") +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title = element_text(size = 20),
          strip.text = element_text(size = 20),
          legend.title = element_blank(), 
          legend.position = "none")
  
  ggsave("Supp_Grid_total.jpeg", grid_total, 
         dpi=300, width = 16, height = 10, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
  
  
  # -------- Combine plots and save ------------
  
  if (i_folders == 2) {
    # Crystallized
    crystallized_aper <- plot_grid(
      myplots[[4]],myplots[[3]], myplots[[2]], myplots[[1]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("cryst4plots_aper.jpeg", crystallized_aper, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    crystallized_per <- plot_grid(
      myplots[[16]],myplots[[15]], myplots[[14]], myplots[[13]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("cryst4plots_per.jpeg", crystallized_per, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    crystallized_total <- plot_grid(
      myplots[[28]],myplots[[27]], myplots[[26]], myplots[[25]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("cryst4plots_total.jpeg", crystallized_total, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3]))
    
    # Fluid
    fluid_aper <- plot_grid(
      myplots[[8]],myplots[[7]], myplots[[6]], myplots[[5]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("fluid4plots_aper.jpeg", fluid_aper, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    fluid_per <- plot_grid(
      myplots[[20]],myplots[[19]], myplots[[18]], myplots[[17]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("fluid4plots_per.jpeg", fluid_per, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    fluid_total <- plot_grid(
      myplots[[32]],myplots[[31]], myplots[[30]], myplots[[29]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("fluid4plots_total.jpeg", fluid_total, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3]))
    
    # Sleepiness
    sleep_aper <- plot_grid(
      myplots[[12]],myplots[[11]], myplots[[10]], myplots[[9]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("sleep4plots_aper.jpeg", sleep_aper, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    sleep_per <- plot_grid(
      myplots[[24]],myplots[[23]], myplots[[22]], myplots[[21]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("sleep4plots_per.jpeg", sleep_per, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    sleep_total <- plot_grid(
      myplots[[36]],myplots[[35]], myplots[[34]], myplots[[33]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("sleep4plots_total.jpeg", sleep_total, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3]))
    
  } else {
    # Crystallized
    crystallized_aper <- plot_grid(
      myplots[[4]],myplots[[3]], myplots[[2]], myplots[[1]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("cryst4plots_aper.jpeg", crystallized_aper, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    crystallized_per <- plot_grid(
      myplots[[12]],myplots[[11]], myplots[[10]], myplots[[9]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("cryst4plots_per.jpeg", crystallized_per, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    crystallized_total <- plot_grid(
      myplots[[20]],myplots[[19]], myplots[[18]], myplots[[17]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("cryst4plots_total.jpeg", crystallized_total, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    # Fluid
    fluid_aper <- plot_grid(
      myplots[[8]],myplots[[7]], myplots[[6]], myplots[[5]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("fluid4plots_aper.jpeg", fluid_aper, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    fluid_per <- plot_grid(
      myplots[[16]],myplots[[15]], myplots[[14]], myplots[[13]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("fluid4plots_per.jpeg", fluid_per, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
    fluid_total <- plot_grid(
      myplots[[24]],myplots[[23]], myplots[[22]], myplots[[21]], 
      nrow = 2, labels = c('A', 'B', 'C', 'D'), label_size = 12)
    ggsave("fluid4plots_total.jpeg", fluid_total, 
           dpi=300, width = 8, height = 5, path = paste0("Results/Plots_", foldersplits[[i_folders]][3])) 
    
  }
  
  
  # --------------------------------------------------
  # -------------  Get average plots  ----------------   
  # --------------------------------------------------
  
  # make empty data frames for values & permuted values
  mean_vals <- data.frame(matrix(NA, nrow = 30, ncol = length(FileNames)))
  mean_perm_vals <- data.frame(matrix(NA, nrow = 30, ncol = length(FileNames)))
  
  # add column names
  colnames(mean_vals) <- FileNames
  colnames(mean_perm_vals) <- paste0(FileNames, "_perm")
  
  
  for (i in 1:length(myFiles)) { # for i in the number of files
    # Read in file
    data <- read.csv(
      paste0("Results/", myFolders[i_folders], "/", myFiles[i]), header = FALSE)
    
    # Create mean scores for variables
    Mean <- rowMeans(data[,c(1:10)]) 
    mean_vals[i] <- Mean
  }  # end for i in number of files
  
  if (i_folders == 2) {
    
    # ---------------------------------------
    # ---------- Crystallized ---------------
    # ---------------------------------------
    
    # get means/sds/SEs
    mean_cryst_aper <- apply(mean_vals[,c(1:4)], 1, mean)
    sd_cryst_aper <- apply(mean_vals[,c(1:4)], 1, sd)
    se_cryst_aper <- sd_cryst_aper/sqrt(4) #standard error
    
    # make mean values
    cryst_means_aper <- as.data.frame(mean_cryst_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SE values
    cryst_SEs_aper <- as.data.frame(se_cryst_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    cryst_means_aper$se <- cryst_SEs_aper$value
    # make upper and lower band of SD
    cryst_means_aper$upper.percentile = cryst_means_aper$value + 1.96 * cryst_means_aper$se
    cryst_means_aper$lower.percentile = cryst_means_aper$value - 1.96 * cryst_means_aper$se
    cryst_means_aper$behav <- "Crystallized"
    cryst_means_aper$signal <- "Aperiodic"
    cryst_means_aper$sample <- "Full"
    
    # get means/sds/SEs
    mean_cryst_per <- apply(mean_vals[,c(13:16)], 1, mean)
    sd_cryst_per <- apply(mean_vals[,c(13:16)], 1, sd)
    se_cryst_per <- sd_cryst_per/sqrt(4) #standard error
    
    # make mean values
    cryst_means_per <- as.data.frame(mean_cryst_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    cryst_SEs_per <- as.data.frame(se_cryst_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    cryst_means_per$se <- cryst_SEs_per$value
    # make upper and lower band of SD
    cryst_means_per$upper.percentile = cryst_means_per$value + 1.96 * cryst_means_per$se
    cryst_means_per$lower.percentile = cryst_means_per$value - 1.96 * cryst_means_per$se
    cryst_means_per$behav <- "Crystallized"
    cryst_means_per$signal <- "Periodic"
    cryst_means_per$sample <- "Full"
    
    
    # get means/sds/SEs
    mean_cryst_total <- apply(mean_vals[,c(25:28)], 1, mean)
    sd_cryst_total <- apply(mean_vals[,c(25:28)], 1, sd)
    se_cryst_total <- sd_cryst_total/sqrt(4) #standard error
    
    # make mean values
    cryst_means_total <- as.data.frame(mean_cryst_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    cryst_SEs_total <- as.data.frame(se_cryst_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    cryst_means_total$se <- cryst_SEs_total$value
    # make upper and lower band of SD
    cryst_means_total$upper.percentile = cryst_means_total$value + 1.96 * cryst_means_total$se
    cryst_means_total$lower.percentile = cryst_means_total$value - 1.96 * cryst_means_total$se
    cryst_means_total$behav <- "Crystallized"
    cryst_means_total$signal <- "Total"
    cryst_means_total$sample <- "Full"
    
    # ------------------------------------
    # -------------- Fluid ---------------
    # ------------------------------------
    
    # get means/sds/SEs
    mean_fluid_aper <- apply(mean_vals[,c(5:8)], 1, mean)
    sd_fluid_aper <- apply(mean_vals[,c(5:8)], 1, sd)
    se_fluid_aper <- sd_fluid_aper/sqrt(4) #standard error
    
    # make mean values
    fluid_means_aper <- as.data.frame(mean_fluid_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SE values
    fluid_SEs_aper <- as.data.frame(se_fluid_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    fluid_means_aper$se <- fluid_SEs_aper$value
    # make upper and lower band of SD
    fluid_means_aper$upper.percentile = fluid_means_aper$value + 1.96 * fluid_means_aper$se
    fluid_means_aper$lower.percentile = fluid_means_aper$value - 1.96 * fluid_means_aper$se
    fluid_means_aper$behav <- "Fluid"
    fluid_means_aper$signal <- "Aperiodic"
    fluid_means_aper$sample <- "Full"
    
    # get means/sds/SEs
    mean_fluid_per <- apply(mean_vals[,c(17:20)], 1, mean)
    sd_fluid_per <- apply(mean_vals[,c(17:20)], 1, sd)
    se_fluid_per <- sd_fluid_per/sqrt(4) #standard error
    
    # make mean values
    fluid_means_per <- as.data.frame(mean_fluid_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    fluid_SEs_per <- as.data.frame(se_fluid_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    fluid_means_per$se <- fluid_SEs_per$value
    # make upper and lower band of SD
    fluid_means_per$upper.percentile = fluid_means_per$value + 1.96 * fluid_means_per$se
    fluid_means_per$lower.percentile = fluid_means_per$value - 1.96 * fluid_means_per$se
    fluid_means_per$behav <- "Fluid"
    fluid_means_per$signal <- "Periodic"
    fluid_means_per$sample <- "Full"
    
    
    # get means/sds/SEs
    mean_fluid_total <- apply(mean_vals[,c(29:32)], 1, mean)
    sd_fluid_total <- apply(mean_vals[,c(29:32)], 1, sd)
    se_fluid_total <- sd_fluid_total/sqrt(4) #standard error
    
    # make mean values
    fluid_means_total <- as.data.frame(mean_fluid_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    fluid_SEs_total <- as.data.frame(se_fluid_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    fluid_means_total$se <- fluid_SEs_total$value
    # make upper and lower band of SD
    fluid_means_total$upper.percentile = fluid_means_total$value + 1.96 * fluid_means_total$se
    fluid_means_total$lower.percentile = fluid_means_total$value - 1.96 * fluid_means_total$se
    fluid_means_total$behav <- "Fluid"
    fluid_means_total$signal <- "Total"
    fluid_means_total$sample <- "Full"
    
    
    # -----------------------------------
    # ----------- Sleepiness ------------
    # -----------------------------------
    
    # get means/sds/SEs
    mean_sleep_aper <- apply(mean_vals[,c(9:12)], 1, mean)
    sd_sleep_aper <- apply(mean_vals[,c(9:12)], 1, sd)
    se_sleep_aper <- sd_sleep_aper/sqrt(4) #standard error
    
    # make mean values
    sleep_means_aper <- as.data.frame(mean_sleep_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SE values
    sleep_SEs_aper <- as.data.frame(se_sleep_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    sleep_means_aper$se <- sleep_SEs_aper$value
    # make upper and lower band of SD
    sleep_means_aper$upper.percentile = sleep_means_aper$value + 1.96 * sleep_means_aper$se
    sleep_means_aper$lower.percentile = sleep_means_aper$value - 1.96 * sleep_means_aper$se
    sleep_means_aper$behav <- "Sleepiness"
    sleep_means_aper$signal <- "Aperiodic"
    sleep_means_aper$sample <- "Full"
    
    # get means/sds/SEs
    mean_sleep_per <- apply(mean_vals[,c(21:24)], 1, mean)
    sd_sleep_per <- apply(mean_vals[,c(21:24)], 1, sd)
    se_sleep_per <- sd_sleep_per/sqrt(4) #standard error
    
    # make mean values
    sleep_means_per <- as.data.frame(mean_sleep_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    sleep_SEs_per <- as.data.frame(se_sleep_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    sleep_means_per$se <- sleep_SEs_per$value
    # make upper and lower band of SD
    sleep_means_per$upper.percentile = sleep_means_per$value + 1.96 * sleep_means_per$se
    sleep_means_per$lower.percentile = sleep_means_per$value - 1.96 * sleep_means_per$se
    sleep_means_per$behav <- "Sleepiness"
    sleep_means_per$signal <- "Periodic"
    sleep_means_per$sample <- "Full"
    
    
    # get means/sds/SEs
    mean_sleep_total <- apply(mean_vals[,c(33:36)], 1, mean)
    sd_sleep_total <- apply(mean_vals[,c(33:36)], 1, sd)
    se_sleep_total <- sd_sleep_total/sqrt(4) #standard error
    
    # make mean values
    sleep_means_total <- as.data.frame(mean_sleep_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    sleep_SEs_total <- as.data.frame(se_sleep_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    sleep_means_total$se <- sleep_SEs_total$value
    # make upper and lower band of SD
    sleep_means_total$upper.percentile = sleep_means_total$value + 1.96 * sleep_means_total$se
    sleep_means_total$lower.percentile = sleep_means_total$value - 1.96 * sleep_means_total$se
    sleep_means_total$behav <- "Sleepiness"
    sleep_means_total$signal <- "Total"
    sleep_means_total$sample <- "Full"
    
    
    combined <- rbind(fluid_means_total, fluid_means_per, fluid_means_aper,
                      cryst_means_total, cryst_means_per, cryst_means_aper,
                      sleep_means_total, sleep_means_per, sleep_means_aper)
    
    # Make list of color values for plots
    color_vals <- c("#FFC000","#4F81BD","#D74E2D")
    
    
    fileid_colors <- setNames(color_vals, unique(combined$behav))
    
    mean_plot <- ggplot(combined, aes(x = ID, y = value, color = behav)) +
      geom_line(linewidth = 0.75) +
      geom_ribbon(aes(ymin = lower.percentile, ymax = upper.percentile, fill = behav),
                  linetype = "blank", alpha = 0.40) +
      scale_y_continuous(breaks = seq(-0.20, 0.25, 0.20), limits = c(-0.20, 0.25)) +
      geom_hline(yintercept = 0.20, linetype = "dashed", color = "#8c8c8c") +
      scale_color_manual(name = "Behav", values = fileid_colors) +
      scale_fill_manual(name = "Behav", values = fileid_colors) +
      facet_grid(factor(signal, c("Total", "Periodic", "Aperiodic")) ~ behav) +
      labs(x = "Frequency (Hz)", y = "Correlation") +
      theme(axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            axis.title = element_text(size = 20),
            strip.text = element_text(size = 20),
            legend.title = element_blank(), 
            legend.position = "none")
    
    ggsave("Mean_Grid.jpeg", mean_plot, 
           dpi=300, width = 16, height = 10, path = paste0("Results/Plots_", foldersplits[[i_folders]][3]))
    
  } else {
    
    # ---------------------------------------
    # ---------- Crystallized ---------------
    # ---------------------------------------
    
    # get means/sds/SEs
    mean_cryst_aper <- apply(mean_vals[,c(1:4)], 1, mean)
    sd_cryst_aper <- apply(mean_vals[,c(1:4)], 1, sd)
    se_cryst_aper <- sd_cryst_aper/sqrt(4) #standard error
    
    # make mean values
    cryst_means_aper <- as.data.frame(mean_cryst_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SE values
    cryst_SEs_aper <- as.data.frame(se_cryst_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    cryst_means_aper$se <- cryst_SEs_aper$value
    # make upper and lower band of SD
    cryst_means_aper$upper.percentile = cryst_means_aper$value + 1.96 * cryst_means_aper$se
    cryst_means_aper$lower.percentile = cryst_means_aper$value - 1.96 * cryst_means_aper$se
    cryst_means_aper$behav <- "Crystallized"
    cryst_means_aper$signal <- "Aperiodic"
    cryst_means_aper$sample <- Sample[1]
    
    # get means/sds/SEs
    mean_cryst_per <- apply(mean_vals[,c(9:12)], 1, mean)
    sd_cryst_per <- apply(mean_vals[,c(9:12)], 1, sd)
    se_cryst_per <- sd_cryst_per/sqrt(4) #standard error
    
    # make mean values
    cryst_means_per <- as.data.frame(mean_cryst_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    cryst_SEs_per <- as.data.frame(se_cryst_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    cryst_means_per$se <- cryst_SEs_per$value
    # make upper and lower band of SD
    cryst_means_per$upper.percentile = cryst_means_per$value + 1.96 * cryst_means_per$se
    cryst_means_per$lower.percentile = cryst_means_per$value - 1.96 * cryst_means_per$se
    cryst_means_per$behav <- "Crystallized"
    cryst_means_per$signal <- "Periodic"
    cryst_means_per$sample <- Sample[1]
    
    
    # get means/sds/SEs
    mean_cryst_total <- apply(mean_vals[,c(17:20)], 1, mean)
    sd_cryst_total <- apply(mean_vals[,c(17:20)], 1, sd)
    se_cryst_total <- sd_cryst_total/sqrt(4) #standard error
    
    # make mean values
    cryst_means_total <- as.data.frame(mean_cryst_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    cryst_SEs_total <- as.data.frame(se_cryst_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    cryst_means_total$se <- cryst_SEs_total$value
    # make upper and lower band of SD
    cryst_means_total$upper.percentile = cryst_means_total$value + 1.96 * cryst_means_total$se
    cryst_means_total$lower.percentile = cryst_means_total$value - 1.96 * cryst_means_total$se
    cryst_means_total$behav <- "Crystallized"
    cryst_means_total$signal <- "Total"
    cryst_means_total$sample <- Sample[1]
    
    # ------------------------------------
    # -------------- Fluid ---------------
    # ------------------------------------
    
    # get means/sds/SEs
    mean_fluid_aper <- apply(mean_vals[,c(5:8)], 1, mean)
    sd_fluid_aper <- apply(mean_vals[,c(5:8)], 1, sd)
    se_fluid_aper <- sd_fluid_aper/sqrt(4) #standard error
    
    # make mean values
    fluid_means_aper <- as.data.frame(mean_fluid_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SE values
    fluid_SEs_aper <- as.data.frame(se_fluid_aper) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    fluid_means_aper$se <- fluid_SEs_aper$value
    # make upper and lower band of SD
    fluid_means_aper$upper.percentile = fluid_means_aper$value + 1.96 * fluid_means_aper$se
    fluid_means_aper$lower.percentile = fluid_means_aper$value - 1.96 * fluid_means_aper$se
    fluid_means_aper$behav <- "Fluid"
    fluid_means_aper$signal <- "Aperiodic"
    fluid_means_aper$sample <- Sample[1]
    
    # get means/sds/SEs
    mean_fluid_per <- apply(mean_vals[,c(13:16)], 1, mean)
    sd_fluid_per <- apply(mean_vals[,c(13:16)], 1, sd)
    se_fluid_per <- sd_fluid_per/sqrt(4) #standard error
    
    # make mean values
    fluid_means_per <- as.data.frame(mean_fluid_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    fluid_SEs_per <- as.data.frame(se_fluid_per) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    fluid_means_per$se <- fluid_SEs_per$value
    # make upper and lower band of SD
    fluid_means_per$upper.percentile = fluid_means_per$value + 1.96 * fluid_means_per$se
    fluid_means_per$lower.percentile = fluid_means_per$value - 1.96 * fluid_means_per$se
    fluid_means_per$behav <- "Fluid"
    fluid_means_per$signal <- "Periodic"
    fluid_means_per$sample <- Sample[1]
    
    
    # get means/sds/SEs
    mean_fluid_total <- apply(mean_vals[,c(21:24)], 1, mean)
    sd_fluid_total <- apply(mean_vals[,c(21:24)], 1, sd)
    se_fluid_total <- sd_fluid_total/sqrt(4) #standard error
    
    # make mean values
    fluid_means_total <- as.data.frame(mean_fluid_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    # make SD values
    fluid_SEs_total <- as.data.frame(se_fluid_total) %>%
      tibble::rowid_to_column("ID") %>%  melt(id = "ID")
    
    # add SDs & perm SDs to the main dataset
    fluid_means_total$se <- fluid_SEs_total$value
    # make upper and lower band of SD
    fluid_means_total$upper.percentile = fluid_means_total$value + 1.96 * fluid_means_total$se
    fluid_means_total$lower.percentile = fluid_means_total$value - 1.96 * fluid_means_total$se
    fluid_means_total$behav <- "Fluid"
    fluid_means_total$signal <- "Total"
    fluid_means_total$sample <- Sample[1]
    
    
    
    combined <- rbind(fluid_means_total, fluid_means_per, fluid_means_aper,
                      cryst_means_total, cryst_means_per, cryst_means_aper)
    
    # Make list of color values for plots
    color_vals <- c("#FFC000","#4F81BD")
    
    
    fileid_colors <- setNames(color_vals, unique(combined$behav))
    
    mean_plot <- ggplot(combined, aes(x = ID, y = value, color = behav)) +
      geom_line(linewidth = 0.75) +
      geom_ribbon(aes(ymin = lower.percentile, ymax = upper.percentile, fill = behav),
                  linetype = "blank", alpha = 0.40) +
      scale_y_continuous(breaks = seq(-0.20, 0.25, 0.20), limits = c(-0.20, 0.25)) +
      geom_hline(yintercept = 0.20, linetype = "dashed", color = "#8c8c8c") +
      scale_color_manual(name = "Behav", values = fileid_colors) +
      scale_fill_manual(name = "Behav", values = fileid_colors) +
      facet_grid(factor(signal, c("Total", "Periodic", "Aperiodic")) ~ behav) +
      labs(x = "Frequency (Hz)", y = "Correlation") +
      theme(axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            axis.title = element_text(size = 20),
            strip.text = element_text(size = 20),
            legend.title = element_blank(), 
            legend.position = "none")
    
    ggsave("Mean_Grid.jpeg", mean_plot, 
           dpi=300, width = 16, height = 10, path = paste0("Results/Plots_", foldersplits[[i_folders]][3]))
    
  }
}

print("Done.")