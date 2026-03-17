library(tidyverse)

rm(list = ls())

#read in output datafiles

oscaa_controls <- read.csv("oscaa_controls_output.csv")

qsm2_controls <- read.csv("qsm2_controls_output.csv")

cheats_os <- read.csv("os_cheats_output.csv")

cheats_qsm2 <- read.csv("qsm2_cheats_output.csv")

cheats_pk <- read.csv("qsm2_pk_cheats_output.csv")

#rbind all
output <- rbind(oscaa_controls, qsm2_controls, cheats_os, cheats_qsm2, cheats_pk)

################################################################
#remove NA values
#log transform data to normalize
#also remove outliers

# Remove NA values
output <- output %>%
  filter(!is.na(CrossNNDist) & !is.na(RedNNDist) & !is.na(GreenNNDist))

#log transform
output <- output %>%
  mutate(
    Log_CrossNNDist = log10(CrossNNDist),
    Log_RedNNDist = log10(RedNNDist),
    Log_GreenNNDist = log10(GreenNNDist)
  )

#Define a function to remove outliers
#Do by image
remove_outliers <- function(df) {
  df %>%
    group_by(Names) %>%
    filter(
      Log_CrossNNDist >= quantile(Log_CrossNNDist, 0.25, na.rm = TRUE) - 1.5 * IQR(Log_CrossNNDist, na.rm = TRUE),
      Log_CrossNNDist <= quantile(Log_CrossNNDist, 0.75, na.rm = TRUE) + 1.5 * IQR(Log_CrossNNDist, na.rm = TRUE),
      Log_GreenNNDist >= quantile(Log_GreenNNDist, 0.25, na.rm = TRUE) - 1.5 * IQR(Log_GreenNNDist, na.rm = TRUE),
      Log_GreenNNDist <= quantile(Log_GreenNNDist, 0.75, na.rm = TRUE) + 1.5 * IQR(Log_GreenNNDist, na.rm = TRUE),
      Log_RedNNDist >= quantile(Log_RedNNDist, 0.25, na.rm = TRUE) - 1.5 * IQR(Log_RedNNDist, na.rm = TRUE),
      Log_RedNNDist <= quantile(Log_RedNNDist, 0.75, na.rm = TRUE) + 1.5 * IQR(Log_RedNNDist, na.rm = TRUE)
    ) %>%
    ungroup()
}

# Remove outliers
filtered_output <- remove_outliers(output)


################################################################
#Calculate mean values by image, not raw data

mean_output <- filtered_output %>%
  group_by(Names, Treatment) %>%
  summarize(
    cross_mean = mean(Log_CrossNNDist, na.rm = TRUE),
    green_mean = mean(Log_GreenNNDist, na.rm = TRUE),
    red_mean = mean(Log_RedNNDist, na.rm = TRUE)
  )

meanofMEAN <- mean_output %>%
  group_by(Treatment) %>%
  summarize(
    cross_meanMEAN = mean(cross_mean, na.rm = TRUE),
    green_meanMEAN = mean(green_mean, na.rm = TRUE),
    red_meanMEAN = mean(red_mean, na.rm = TRUE)
  )

write_csv(meanofMEAN,"logtransformed_meansbyimage.csv")

#single tailed t-test to determine if
#crossdist is higher than red or green by treatment
results <- do.call(rbind, lapply(split(mean_output, mean_output$Treatment), function(df) {
  # Perform t-tests
  cross_vs_red <- t.test(df$cross_mean, df$red_mean, alternative = "greater", paired = FALSE)
  cross_vs_green <- t.test(df$cross_mean, df$green_mean, alternative = "greater", paired = FALSE)
  
  # Collect results
  data.frame(
    Treatment = unique(df$Treatment),
    Cross_vs_Red_p = cross_vs_red$p.value,
    Cross_vs_Green_p = cross_vs_green$p.value
  )
}))




#divide up groups so can perform holm adjustment
control_asocial_results <- results %>% filter(Treatment %in% c("Asocial: PAO1 GFP-PAO1 mCherry",
                                                               "Asocial: PAO1 GFP-Δssg mCherry",
                                                               "Asocial: Δssg GFP-Δssg mCherry"))

control_social_results <- results %>% filter(Treatment %in% c("Social: PAO1 GFP-PAO1 mCherry",
                                                           "Social: PAO1 GFP-Δssg mCherry",
                                                           "Social: Δssg GFP-Δssg mCherry"))


cheats_asocial_results <- results %>% filter(Treatment %in% c("Asocial: PAO1 GFP-ΔssgΔlasRΔrhlR mCherry",
                                                              "Asocial: Δssg GFP-ΔssgΔlasRΔrhlR mCherry",
                                                              "Asocial: PAO1 GFP-PAO1ΔlasRΔrhlR mCherry",
                                                              "Asocial: Δssg GFP-PAO1ΔlasRΔrhlR mCherry"))

cheats_social_results <- results %>% filter(Treatment %in% c("Social: PAO1 GFP-ΔssgΔlasRΔrhlR mCherry",
                                                              "Social: Δssg GFP-ΔssgΔlasRΔrhlR mCherry",
                                                              "Social: PAO1 GFP-PAO1ΔlasRΔrhlR mCherry",
                                                              "Social: Δssg GFP-PAO1ΔlasRΔrhlR mCherry"))

cheats_pk_results <- results %>% filter(Treatment %in% c("PK: PAO1 GFP-ΔssgΔlasRΔrhlR mCherry",
                                                             "PK: Δssg GFP-ΔssgΔlasRΔrhlR mCherry",
                                                             "PK: PAO1 GFP-PAO1ΔlasRΔrhlR mCherry",
                                                             "PK: Δssg GFP-PAO1ΔlasRΔrhlR mCherry"))

                                          

# List of data frames
data_frames <- list(control_asocial_results, control_social_results,
                    cheats_asocial_results, cheats_social_results,
                    cheats_pk_results)

# Loop through each data frame and apply the Holm adjustment
for (i in seq_along(data_frames)) {
  df <- data_frames[[i]]
  
  df$Adjusted_P_value_CrossvsRed <- p.adjust(df$Cross_vs_Red_p, method = "holm")
  df$Adjusted_P_value_CrossvsGreen <- p.adjust(df$Cross_vs_Green_p, method = "holm")
  
  # Save the adjusted data frame back to the list
  data_frames[[i]] <- df
}

# Combine the adjusted data frames back into the results data frame
adjusted_results <- bind_rows(data_frames)
# Now `adjusted_results` contains the Holm-adjusted p-values only for the relevant groups


write_csv(adjusted_results,"singletailttest_logtrasnformed_holmadjusted.csv")
