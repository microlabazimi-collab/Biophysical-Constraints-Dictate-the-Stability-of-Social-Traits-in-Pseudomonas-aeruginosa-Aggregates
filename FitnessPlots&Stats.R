library(tidyverse)
library(emmeans)

data <- read_csv("cheatingfitnesssummary.csv")

list <- c("WT vs ∆QS","WT vs ∆OSA∆QS","∆OSA vs ∆OSA∆QS","∆OSA vs ∆QS")

data1 <- data %>%
  filter(Treatment %in% list)

#rename to keep consistent with paper
data1$Treatment <- gsub("WT", "PAO1", data1$Treatment)
data1$Treatment <- gsub("OSA", "ssg", data1$Treatment)
data1$Treatment <- gsub("QS", "lasR∆rhlR", data1$Treatment)

#make graph
ggplot(data1, aes(x = Starting, y = Fitness, color = Treatment)) +
  geom_point(size = 3) +             # Add points
  geom_smooth(method = "lm", se = FALSE) +  # Add a regression line
  scale_y_log10() +                  # Transform y-axis to log10 scale
  theme_minimal() +                  # Apply a minimal theme
  labs(
    title = "Cheater Fitness vs Starting",
    x = "Starting",
    y = "Fitness (log10 scale)",
    color = "Treatment"
  )


#stats
#TEST 1: are Treatments signif different at Starting 0.1 (aka 1/10)
# Fit the linear model with interaction
model <- lm(Fitness ~ Treatment * Starting, data = data1)

# Create a new data frame for predictions
new_data <- data.frame(
  Treatment = unique(data1$Treatment),
  Starting = 0.1
)

# Predict Fitness and calculate confidence intervals
predictions <- predict(model, newdata = new_data, interval = "confidence")
results <- cbind(new_data, predictions)
print(results)

# Compare treatments at Starting = 0.1
emms <- emmeans(model, ~ Treatment, at = list(Starting = 0.1))
pairwise_results <- pairs(emms)
print(pairwise_results)


#TEST 2: are fitness and starting positively or negatively correlated?
#use pearson
#is the correlation significant?

correlation_results <- data1 %>%
  group_by(Treatment) %>%
  summarize(
    cor = cor(Starting, Fitness, method = "pearson"),
    p_value = cor.test(Starting, Fitness, method = "pearson")$p.value
  )

print(correlation_results)

#TEST 3: are fitness and starting correlated? use model slope

model_results <- data1 %>%
  group_by(Treatment) %>%
  summarize(
    slope = coef(lm(Fitness ~ Starting))[2], # Extract slope
    p_value = summary(lm(Fitness ~ Starting))$coefficients[2, 4] # Extract p-value
  )

print(model_results)
