library(tidyverse)
library(densityClust)


#raw data!!

################################################################
#read in data

#controls

oscaa_controls <- read.csv("oscaa_controls_output.csv")

qsm2_controls <- read.csv("qsm2_controls_output.csv")

#cheats

cheats_os <- read.csv("os_cheats_output.csv")

cheats_qsm2 <- read.csv("qsm2_cheats_output.csv")

cheats_pk <- read.csv("qsm2_pk_cheats_output.csv")

################################################################

# Remove NA values
oscaa_controls <- oscaa_controls %>%
  filter(!is.na(CrossNNDist) & !is.na(RedNNDist) & !is.na(GreenNNDist))

qsm2_controls <- qsm2_controls %>%
  filter(!is.na(CrossNNDist) & !is.na(RedNNDist) & !is.na(GreenNNDist))

cheats_os <- cheats_os %>%
  filter(!is.na(CrossNNDist) & !is.na(RedNNDist) & !is.na(GreenNNDist))

cheats_pk <- cheats_pk %>%
  filter(!is.na(CrossNNDist) & !is.na(RedNNDist) & !is.na(GreenNNDist))

cheats_qsm2 <- cheats_qsm2 %>%
  filter(!is.na(CrossNNDist) & !is.na(RedNNDist) & !is.na(GreenNNDist))


################################################################

#write function

highlight_density <- function(data, treatment) {
  # Filter data for the current treatment
  data_subset <- data %>% filter(Treatment == treatment)
  
  # Calculate densities
  density_cross <- density(data_subset$CrossNNDist, na.rm = TRUE)
  density_red <- density(data_subset$RedNNDist, na.rm = TRUE)
  density_green <- density(data_subset$GreenNNDist, na.rm = TRUE)
  
  # Create a data frame for plotting
  density_data <- data.frame(
    x = density_cross$x,
    y_cross = density_cross$y,
    y_red = approx(density_red$x, density_red$y, xout = density_cross$x)$y,
    y_green = approx(density_green$x, density_green$y, xout = density_cross$x)$y
  )
  
  # Ensure NA values are treated as 0
  density_data$y_red[is.na(density_data$y_red)] <- 0
  density_data$y_green[is.na(density_data$y_green)] <- 0
  
  # Calculate the area not covered by RedNNDist or GreenNNDist to the right
  density_data$y_highlight <- ifelse(density_data$x > max(density_data$x[density_data$y_red > density_data$y_cross & 
                                                                           density_data$y_green > density_data$y_cross]), 
                                     pmax(0, density_data$y_cross - pmax(density_data$y_red, density_data$y_green)), 
                                     0)
  
  # Calculate the total area under CrossNNDist
  total_area_cross <- sum(density_data$y_cross) * diff(density_data$x[1:2])
  
  # Calculate the highlighted area
  highlighted_area <- sum(density_data$y_highlight) * diff(density_data$x[1:2])
  
  # Calculate the percentage of the highlighted area
  percentage_highlighted <- (highlighted_area / total_area_cross) * 100
  
  # Plot the densities and highlight the area not covered by RedNNDist or GreenNNDist to the right
  p <- ggplot(density_data, aes(x = x)) +
    geom_line(aes(y = y_cross, color = 'CrossNNDist')) +
    geom_line(aes(y = y_red, color = 'RedNNDist')) +
    geom_line(aes(y = y_green, color = 'GreenNNDist')) +
    geom_ribbon(aes(ymin = 0, ymax = y_highlight), fill = 'black', alpha = 0.3) +
    labs(title = paste('Raw Nearest Neighbor Density Plots:\n',treatment, 
                       '\nPercentage of CrossNNDist Above\nRedNNDist and GreenNNDist:', round(percentage_highlighted, 2), '%'),
         x = 'Value', y = 'Density') +
    theme_minimal() +
    scale_color_manual(values = c('CrossNNDist' = 'black', 'RedNNDist' = 'red', 'GreenNNDist' = 'green')) +
    scale_fill_manual(values = c('CrossNNDist' = 'black', 'RedNNDist' = 'red', 'GreenNNDist' = 'green'))

  
  # Save the plot
  ggsave(filename = paste0('Density_Plot_Treatment_', treatment, '.png'), plot = p, bg = 'white')
  
  return(p)
}

################################################################

#run function

# Loop through each Treatment group and plot
for (treatment in unique(oscaa_controls$Treatment)) {
  print(highlight_density(oscaa_controls, treatment))
}

# Loop through each Treatment group and plot
for (treatment in unique(qsm2_controls$Treatment)) {
  print(highlight_density(qsm2_controls, treatment))
}

# Loop through each Treatment group and plot
for (treatment in unique(cheats_os$Treatment)) {
  print(highlight_density(cheats_os, treatment))
}

# Loop through each Treatment group and plot
for (treatment in unique(cheats_qsm2$Treatment)) {
  print(highlight_density(cheats_qsm2, treatment))
}

# Loop through each Treatment group and plot
for (treatment in unique(cheats_pk$Treatment)) {
  print(highlight_density(cheats_pk, treatment))
}



