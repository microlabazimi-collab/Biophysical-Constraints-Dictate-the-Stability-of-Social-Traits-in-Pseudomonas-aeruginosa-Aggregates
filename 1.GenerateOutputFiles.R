library(stringr)
library(ggplot2)
library(spatstat)
library(tidyverse)
library(rstatix)
library(ggsignif)
library(ggpubr)
library(data.table)

rm(list = ls())

############################################################################
##################Define Main Data Output Function##########################
############################################################################

generate_output <- function(R, G) {
  
  # Process Red (R) and Green (G) data frames
  dfR1 <- R %>%
    mutate(col = LOCATION %>% str_sub(4, -2)) %>%
    separate(col, c("X", "Y"), sep = "\\s", convert = TRUE) %>%
    mutate(
      Y = Y %>% str_sub(3, -1),
      color = "red"
    )
  
  dfG1 <- G %>%
    mutate(col = LOCATION %>% str_sub(4, -2)) %>%
    separate(col, c("X", "Y"), sep = "\\s", convert = TRUE) %>%
    mutate(
      Y = Y %>% str_sub(3, -1),
      color = "green"
    )
  
  # Define box size based on image size
  bb <- box3(
    xrange = c(0, 134.95),
    yrange = c(0, 134.95),
    zrange = c(1, 11)
  )
  
  # Combine data and split by experiment and color
  dfA <- bind_rows(dfR1, dfG1)
  dfA_1 <- split(dfA, dfA$EXPERIMENT)
  dfA_2 <- lapply(dfA_1, function(x) split(x, x$color))
  
  # Initialize the output data frame
  output <- data.frame('Names' = character(), 'GreenNNDist' = numeric(), 
                       'RedNNDist' = numeric(), 'CrossNNDist' = numeric(),
                       'RedCount' = integer(), 'GreenCount' = integer(),
                       'Experiment' = character(), 'Treatment' = character())
  
  # Loop through the data
  for (i in seq_along(dfA_2)) {
    loopEnv <- new.env()
    
    for (color in names(dfA_2[[i]])) {
      dist.pp3 <- pp3(as.numeric(dfA_2[[i]][[color]]$X), 
                      as.numeric(dfA_2[[i]][[color]]$Y), 
                      as.numeric(dfA_2[[i]][[color]]$POSITION), bb)
      
      nndist <- nndist.pp3(dist.pp3, k = 1)
      nnmean <- mean(nndist)
      
      assign(paste(color, "dist.pp3", sep = "_"), dist.pp3, envir = loopEnv)
      assign(paste(color, "nndist", sep = "_"), nndist, envir = loopEnv)
      assign(paste(color, "nnmean", sep = "_"), nnmean, envir = loopEnv)
    }
    
    # Cross-distances between red and green cells
    x <- nncross.pp3(loopEnv$red_dist.pp3, loopEnv$green_dist.pp3, k = 1)
    red_count <- nrow(dfA_2[[i]]$red)
    green_count <- nrow(dfA_2[[i]]$green)
    
    max_length <- max(
      length(x$dist),
      length(loopEnv$green_nndist),
      length(loopEnv$red_nndist),
      red_count,
      green_count
    )
    
    # Create newdf, extending each column to the maximum length
    newdf <- tibble(
      "Names" = rep(names(dfA_2)[i], length.out = max_length),
      "GreenNNDist" = rep(loopEnv$green_nndist, length.out = max_length),
      "RedNNDist" = rep(loopEnv$red_nndist, length.out = max_length),
      "CrossNNDist" = c(x$dist, rep(NA, max_length - length(x$dist))),
      "RedCount" = rep(red_count, length.out = max_length),
      "GreenCount" = rep(green_count, length.out = max_length),
      "Experiment" = rep(names(dfA_2)[i], length.out = max_length),
      "Treatment" = rep(str_sub(names(dfA_2)[i], 1, str_locate(names(dfA_2)[i], "_")[1] - 1), 
                        length.out = max_length)
    )
    
    output <- bind_rows(output, newdf)
    gc()
  }
  
  return(output)
}

#############################################################################
####################### Define Plotting Function ############################
#############################################################################

# Define the function
plotting <- function(output) {
  
  for (treatment in unique(output$Treatment)) {
    treatment_data <- output[output$Treatment == treatment, ]
    
    # Calculate the mean and median values
    dist_mean <- mean(treatment_data$CrossNNDist, na.rm = TRUE)
    green_mean <- mean(treatment_data$GreenNNDist, na.rm = TRUE)
    red_mean <- mean(treatment_data$RedNNDist, na.rm = TRUE)
    dist_median <- median(treatment_data$CrossNNDist, na.rm = TRUE)
    
    # Plot the histogram for the 'CrossNNDist' variable
    plot <- ggplot(treatment_data, aes(x = CrossNNDist)) +
      geom_histogram(binwidth = 0.10, boundary = 0, fill = "blue", color = "black") +
      labs(title = paste("Nearest Neighbor XDist:", treatment), x = "Dist", y = "Frequency") +
      theme_minimal() +
      geom_point(aes(x = dist_mean, y = 5000), color = "black", size = 3) +
      geom_point(aes(x = green_mean, y = 5500), color = "green", size = 3) +
      geom_point(aes(x = red_mean, y = 4500), color = "red", size = 3) +
      geom_vline(aes(xintercept = dist_median), color = "black", linetype = "dashed", size = 1)
    
    # Save the plot as a PNG file
    ggsave(filename = paste0("Nearest_Neighbor_XDist", "_", treatment, ".png"), plot = plot, width = 8, height = 6, dpi = 300, bg = "white")
  }
}

#############################################################################
####################### Define T-Test Function ##############################
#############################################################################

# Function to perform t-tests with Bonferroni adjustment
t_tests <- function(output) {
  results <- data.frame(Treatment = character(),
                        Comparison = character(),
                        t_statistic = numeric(),
                        p_value = numeric(),
                        stringsAsFactors = FALSE)
  
  # Loop through unique treatments
  for (treatment in unique(output$Treatment)) {
    treatment_data <- output[output$Treatment == treatment, ]
    
    # Perform t-tests
    test_green_vs_cross <- t.test(treatment_data$GreenNNDist, treatment_data$CrossNNDist, na.rm = TRUE)
    test_red_vs_cross <- t.test(treatment_data$RedNNDist, treatment_data$CrossNNDist, na.rm = TRUE)
    
    # Collect results
    results <- rbind(results, data.frame(
      Treatment = treatment,
      Comparison = "GreenNNDist vs CrossNNDist",
      t_statistic = test_green_vs_cross$statistic,
      p_value = test_green_vs_cross$p.value
    ))
    
    results <- rbind(results, data.frame(
      Treatment = treatment,
      Comparison = "RedNNDist vs CrossNNDist",
      t_statistic = test_red_vs_cross$statistic,
      p_value = test_red_vs_cross$p.value
    ))
  }
  
  # Bonferroni adjustment for p-values
  results$p_value_adj <- p.adjust(results$p_value, method = "bonferroni")
  
  return(results)
}


############################################################################
################Read in Files & Call Functions##############################
############################################################################

################## Controls First ##########################################

############################################################################
###########################pao1 vs ssg OS-CAA###############################
############################################################################
#missing headers, not sure why, need to add in
R <- read.csv("csv_files_useme/controls/ssgvspao1_oscaa_mCherry.csv", header = FALSE)
G <- read.csv("csv_files_useme/controls/ssgvspao1_oscaa_gfp.csv", header = FALSE)

# Define new column names
new_colnames <- c('NAME', 'ASSOCIATION', 'EXPERIMENT', 'IMAGE', 'LOCATION', 'POSITION', 'ZSCORE')

# Assign new column names to both data frames
colnames(R) <- new_colnames
colnames(G) <- new_colnames

#call on function to generate main output data file
output_oscaa_controls <- generate_output(R, G)

#rename treatment for graphing and tables
output_oscaa_controls$Treatment <- gsub("ssg", "Asocial: Δssg GFP-Δssg mCherry", output_oscaa_controls$Treatment)
output_oscaa_controls$Treatment <- gsub("PAO1", "Asocial: PAO1 GFP-PAO1 mCherry", output_oscaa_controls$Treatment)
output_oscaa_controls$Treatment <- gsub("5050mix", "Asocial: PAO1 GFP-Δssg mCherry", output_oscaa_controls$Treatment)

#call plotting function
plotting(output_oscaa_controls)

#call ttest function
oscaa_controls_ttest_results <- t_tests(output_oscaa_controls)

#save output & ttest files as csv
write.csv(output_oscaa_controls,"oscaa_controls_output.csv")
write.csv(oscaa_controls_ttest_results,"oscaa_controls_ttest_results.csv")

#remove R & G data files from R environment
rm(R, G)

############################################################################
#############################pao1 vs ssg QSM2##############################
############################################################################
R <- read.csv("csv_files_useme/controls/pao1vsssg_qsm2_june2023_mC.csv", sep=",", header = TRUE)
G <- read.csv("csv_files_useme/controls/pao1vsssg_qsm2_june2023_gfp.csv", header = TRUE)

#call on function to generate main output data file
output_qsm2_controls <- generate_output(R, G)

#rename treatment for graphing and tables
output_qsm2_controls$Treatment <- gsub("ssggfpmC", "Social: Δssg GFP-Δssg mCherry", output_qsm2_controls$Treatment)
output_qsm2_controls$Treatment <- gsub("pao1gfpmC", "Social: PAO1 GFP-PAO1 mCherry", output_qsm2_controls$Treatment)
output_qsm2_controls$Treatment <- gsub("5050pao1gfp-ssgmC", "Social: PAO1 GFP-Δssg mCherry", output_qsm2_controls$Treatment)

#call plotting function
plotting(output_qsm2_controls)

#call ttest function
qsm2_controls_ttest_results <- t_tests(output_qsm2_controls)

#save output & ttest files as csv
write.csv(output_qsm2_controls,"qsm2_controls_output.csv")
write.csv(qsm2_controls_ttest_results,"qsm2_controls_ttest_results.csv")

#remove R & G data files from R environment
rm(R, G)

############################################################################
###########################pao1 vs wbpL QSM2################################
############################################################################
R <- read.csv("csv_files_useme/controls/12022022_wbpLvspao1_mCherry.csv", header = TRUE)
G <- read.csv("csv_files_useme/controls/12022022_wbpLvspao1_gfp.csv", header = TRUE)

#call on function to generate main output data file
output_qsm2_wbpLcontrols <- generate_output(R, G)

#rename treatment for graphing and tables
output_qsm2_wbpLcontrols$Treatment <- gsub("m9m-bsa-wbpL", "Social: ΔwbpL GFP-ΔwbpL mCherry", output_qsm2_wbpLcontrols$Treatment)
output_qsm2_wbpLcontrols$Treatment <- gsub("m9m-bsa-pao1", "Social: PAO1 GFP-PAO1 mCherry (wbpL experiment version)", output_qsm2_wbpLcontrols$Treatment)
output_qsm2_wbpLcontrols$Treatment <- gsub("m9m-bsa-pao1gfpwbpLmC", "Social: PAO1 GFP-ΔwbpL mCherry", output_qsm2_wbpLcontrols$Treatment)

#call plotting function
plotting(output_qsm2_wbpLcontrols)

#call ttest function
qsm2_wbpLcontrols_ttest_results <- t_tests(output_qsm2_wbpLcontrols)

#save output & ttest files as csv
write.csv(output_qsm2_wbpLcontrols,"qsm2_wbpLcontrols_output.csv")
write.csv(qsm2_wbpLcontrols_ttest_results,"qsm2_wbpLcontrols_ttest_results.csv")

#remove R & G data files from R environment
rm(R, G)


######################### Read in Experimental Files #######################


############################################################################
############################OS Cheats All###################################
############################################################################
R <- read.csv("csv_files_useme/oscaa_cheats/06082022_oscaacheats_mCherry.csv")
G <- read.csv("csv_files_useme/oscaa_cheats/06082022_oscaacheats_gfp.csv")

#call on function to generate main output data file
output_oscaa_cheats <- generate_output(R, G)

#unique <- unique(output_oscaa_cheats$Treatment)

#rename treatment for graphing and tables
output_oscaa_cheats$Treatment <- gsub("pao150-50lasRrhlR", "Asocial: PAO1 GFP-PAO1ΔlasRΔrhlR mCherry", output_oscaa_cheats$Treatment)
output_oscaa_cheats$Treatment <- gsub("pao150-50lasRrhlRssg", "Asocial: PAO1 GFP-ΔssgΔlasRΔrhlR mCherry", output_oscaa_cheats$Treatment)
output_oscaa_cheats$Treatment <- gsub("ssg50-50lasRrhlRssg", "Asocial: Δssg GFP-ΔssgΔlasRΔrhlR mCherry", output_oscaa_cheats$Treatment)
output_oscaa_cheats$Treatment <- gsub("ssg50-50lasRrhlR", "Asocial: Δssg GFP-PAO1ΔlasRΔrhlR mCherry", output_oscaa_cheats$Treatment)

#call plotting function
plotting(output_oscaa_cheats)

#call ttest function
oscaa_cheats_ttest_results <- t_tests(output_oscaa_cheats)

#save output & ttest files as csv
write.csv(output_oscaa_cheats,"os_cheats_output.csv")
write.csv(oscaa_cheats_ttest_results,"os_cheats_ttest_results.csv")

#remove R & G data files from R environment
rm(R, G)

############################################################################
########################## QSM2 Cheats #####################################
############################################################################
#exclude tryptone & pk results
#read in pk next

#R1/G1 = pao1 vs lasRrhlR

#mixed imaging day contains contains non pk pao1 vs lasRrhlR
R1a <- read.csv("csv_files_useme/cheating&pk/mixedimaging02092023_mC.csv")
G1a <- read.csv("csv_files_useme/cheating&pk/mixedimaging02092023_GFP.csv")

#keep only pao1 gfp vs lasRrhlRmc from mixed imaging
R1 <- R1a[grepl("pao1gfp50-50lasRrhlRmC-bsa1caa0.05", R1a$EXPERIMENT), ]
G1 <- G1a[grepl("pao1gfp50-50lasRrhlRmC-bsa1caa0.05", G1a$EXPERIMENT), ]

#R2/G2 = ssg vs lasRrhlRssg
R2a <- read.csv("csv_files_useme/cheating&pk/ssggfpvslasRrhlRssg03022023-mC.csv")
G2a <- read.csv("csv_files_useme/cheating&pk/ssggfpvslasRrhlRssg03022023-gfp.csv")

#keep only 50:50 data
R2 <- R2a[grepl("ssggfp50-50lasRrhlRssgmC-bsa1caa0.05", R2a$EXPERIMENT), ]
G2 <- G2a[grepl("ssggfp50-50lasRrhlRssgmC-bsa1caa0.05", G2a$EXPERIMENT), ]

#R3/G3 = ssg vs lasRrhlR
R3a <- read.csv("csv_files_useme/cheating&pk/ssggfplasRrhlRmC_bsa1caa0.05_mCherry.csv")
G3a <- read.csv("csv_files_useme/cheating&pk/ssggfplasRrhlRmC_bsa1caa0.05_gfp.csv")

#keep only 50:50 data
R3 <- R3a[grepl("ssggfp100-1lasRrhlRmC-bsa1caa0.05", R3a$EXPERIMENT), ]
G3 <- G3a[grepl("ssggfp100-1lasRrhlRmC-bsa1caa0.05", G3a$EXPERIMENT), ]

#R4/G4 = pao1 vs lasRrhlRssg
R4a <- read.csv("csv_files_useme/cheating&pk/pao1gfplasRrhlRssgmC_01262023_mC.csv")
G4a <- read.csv("csv_files_useme/cheating&pk/pao1gfplasRrhlRssgmC_01262023_gfp.csv")

#keep only 50:50 data
R4 <- R4a[grepl("pao1gfp50-50lasRrhlRssg-bsa1caa0.05", R4a$EXPERIMENT), ]
G4 <- G4a[grepl("pao1gfp50-50lasRrhlRssg-bsa1caa0.05", G4a$EXPERIMENT), ]

R <- rbind(R1,R2,R3,R4)
G <- rbind(G1,G2,G3,G4)

#a <- unique(R$EXPERIMENT) #there should be 24
#b <- unique(G$EXPERIMENT) #there should be 24

#remove intermediate R and G files
rm(R1,R1a,R2,R2a,R3,R3a,R4,R4a,G1,G1a,G2,G2a,G3,G3a,G4,G4a)

#call on function to generate main output data file
output_qsm2_cheats <- generate_output(R, G)

#rename treatment for graphing and tables
output_qsm2_cheats$Treatment <- gsub("pao1gfp50-50lasRrhlRmC-bsa1caa0.05", 
                                     "Social: PAO1 GFP-PAO1ΔlasRΔrhlR mCherry", output_qsm2_cheats$Treatment)
output_qsm2_cheats$Treatment <- gsub("pao1gfp50-50lasRrhlRssg-bsa1caa0.05", 
                                     "Social: PAO1 GFP-ΔssgΔlasRΔrhlR mCherry", output_qsm2_cheats$Treatment)
output_qsm2_cheats$Treatment <- gsub("ssggfp50-50lasRrhlRssgmC-bsa1caa0.05", 
                                     "Social: Δssg GFP-ΔssgΔlasRΔrhlR mCherry", output_qsm2_cheats$Treatment)
output_qsm2_cheats$Treatment <- gsub("ssggfp100-1lasRrhlRmC-bsa1caa0.05", 
                                     "Social: Δssg GFP-PAO1ΔlasRΔrhlR mCherry", output_qsm2_cheats$Treatment)

#call plotting function
plotting(output_qsm2_cheats)

#call ttest function
qsm2_cheats_ttest_results <- t_tests(output_qsm2_cheats)

#save output & ttest files as csv
write.csv(output_qsm2_cheats,"qsm2_cheats_output.csv")
write.csv(qsm2_cheats_ttest_results,"qsm2_cheats_ttest_results.csv")

#remove R & G data files from R environment
rm(R, G)

############################################################################
############################ QSM2 pk 0.5 Cheats ############################
############################################################################

#code in progress

#R1/G1 = pao1 vs lasRrhlR

R1a <- read.csv("csv_files_useme/cheating&pk/pao1gfplasRrhlRmC_02012023_mC.csv")[,-8]
G1a <- read.csv("csv_files_useme/cheating&pk/pao1gfplasRrhlRmC_02012023_gfp.csv")[,-8]

#keep only pk test from pao1 gfp vs lasRrhlRmc imaging
R1 <- R1a[grepl("pao1gfp50-50lasRrhlRmC-pk0.5", R1a$EXPERIMENT), ]
G1 <- G1a[grepl("pao1gfp50-50lasRrhlRmC-pk0.5", G1a$EXPERIMENT), ]

#R2/G2 = ssg vs lasRrhlRssg
R2a <- read.csv("csv_files_useme/cheating&pk/ssggfpvslasRrhlRssg03022023-mC.csv")[,-8]
G2a <- read.csv("csv_files_useme/cheating&pk/ssggfpvslasRrhlRssg03022023-gfp.csv")[,-8]

#keep only pk data
R2 <- R2a[grepl("ssggfp50-50lasRrhlRssgmC-pk0.5", R2a$EXPERIMENT), ]
G2 <- G2a[grepl("ssggfp50-50lasRrhlRssgmC-pk0.5", G2a$EXPERIMENT), ]

#R3/G3 = ssg vs lasRrhlR
R3a <- read.csv("csv_files_useme/cheating&pk/ssggfplasRrhlRmC_pk_mCherry.csv")[,-8]
G3a <- read.csv("csv_files_useme/cheating&pk/ssggfplasRrhlRmC_pk_gfp.csv")[,-8]

#keep only 0.5 ug/mL pk data
R3 <- R3a[grepl("ssggfp50-50lasRrhlRmC-pk0.5", R3a$EXPERIMENT), ]
G3 <- G3a[grepl("ssggfp50-50lasRrhlRmC-pk0.5", G3a$EXPERIMENT), ]

#R4/G4 = pao1 vs lasRrhlRssg
R4a <- read.csv("csv_files_useme/cheating&pk/pao1gfp-lasRrhlRssgmC-pk_mC.csv", header = FALSE)
G4a <- read.csv("csv_files_useme/cheating&pk/pao1gfp-lasRrhlRssgmC-pk_gfp.csv", header = FALSE)

#dumb formatting but necessary!

process_dataframe <- function(df,df1) {
  # Set the first row as column names and remove it from the data
  colnames(df) <- df[1, ]
  df <- df[-1, ]
  # Concatenate EXPERIMENT and IMAGE columns with a separator
  df$EXPERIMENT <- paste(df$EXPERIMENT, df$IMAGE, sep = "")
  # Concatenate LOCATION and POSITION columns with a separator
  df$IMAGE <- paste(df$LOCATION, df$POSITION, sep = "")
  # Copy ZSCORE column values into LOCATION column
  df$LOCATION <- paste(df$ZSCORE)
  # Rename column 8 to "Replace_Position" and update POSITION
  colnames(df)[8] <- "Replace_Position"
  df$POSITION <- paste(df$Replace_Position)
  # Rename column 9 to "Replace_Zscore" and update ZSCORE
  colnames(df)[9] <- "Replace_Zscore"
  df$ZSCORE <- paste(df$Replace_Zscore)
  # Select and return the required columns
  DF <- df %>% select(NAME, ASSOCIATION, EXPERIMENT, IMAGE, LOCATION, POSITION, ZSCORE)
  
  assign(df1, DF, envir = .GlobalEnv)
}

#Run all the function
process_dataframe(R4a, "R4b")
process_dataframe(G4a, "G4b")

#keep only 50:50 data
R4 <- R4b[grepl("pao1gfp50-50lasRrhlRssg-pk0.5", R4b$EXPERIMENT), ]
G4 <- G4b[grepl("pao1gfp50-50lasRrhlRssg-pk0.5", G4b$EXPERIMENT), ]


R <- rbind(R1,R2,R3,R4)
G <- rbind(G1,G2,G3,G4)

#a <- unique(R$EXPERIMENT) #there should be 24
#b <- unique(G$EXPERIMENT) #there should be 24

rm(R1,R1a,R2,R2a,R3,R3a,R4,R4a,G1,G1a,G2,G2a,G3,G3a,G4,G4a)

#call on function to generate main output data file
output_qsm2_pk_cheats <- generate_output(R, G)

unique <- unique(output_qsm2_pk_cheats$Treatment)

#rename treatment for graphing and tables
output_qsm2_pk_cheats$Treatment <- gsub("pao1gfp50-50lasRrhlRmC-pk0.5", 
                                     "PK: PAO1 GFP-PAO1ΔlasRΔrhlR mCherry", output_qsm2_pk_cheats$Treatment)
output_qsm2_pk_cheats$Treatment <- gsub("pao1gfp50-50lasRrhlRssg-pk0.5", 
                                     "PK: PAO1 GFP-ΔssgΔlasRΔrhlR mCherry", output_qsm2_pk_cheats$Treatment)
output_qsm2_pk_cheats$Treatment <- gsub("ssggfp50-50lasRrhlRssgmC-pk0.5", 
                                     "PK: Δssg GFP-ΔssgΔlasRΔrhlR mCherry", output_qsm2_pk_cheats$Treatment)
output_qsm2_pk_cheats$Treatment <- gsub("ssggfp50-50lasRrhlRmC-pk0.5", 
                                     "PK: Δssg GFP-PAO1ΔlasRΔrhlR mCherry", output_qsm2_pk_cheats$Treatment)

#call plotting function
plotting(output_qsm2_pk_cheats)

#call ttest function
qsm2_pk_cheats_ttest_results <- t_tests(output_qsm2_pk_cheats)

#save output & ttest files as csv
write.csv(output_qsm2_pk_cheats,"qsm2_pk_cheats_output.csv")
write.csv(qsm2_pk_cheats_ttest_results,"qsm2_pk_cheats_ttest_results.csv")

#remove R & G data files from R environment
rm(R, G)







