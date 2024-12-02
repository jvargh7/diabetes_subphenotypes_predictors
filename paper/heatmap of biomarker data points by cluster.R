rm(list = ls());gc();source(".Rprofile")

library(dplyr)
library(ggplot2)
library(tidyr)

analytic_df <- readRDS(paste0(path_diabetes_subphenotypes_predictors_folder,"/working/cleaned/dsppre01_analytic df.RDS")) %>% 
  dplyr::filter(!is.na(cluster_study_id)) %>% 
  mutate(t = age - dmagediag) %>% 
  dplyr::filter(t >= -15, t < 0)

biomarkers <- c("hba1c", "insulinf", "glucosef", "glucose2h", "tgl", 
                "hdlc", "ldlc", "serumcreatinine", "urinecreatinine", 
                "egfr", "apo_a", "apo_b", "uric_acid", "vldlc")

# Pivot the data longer so each row is a single observation of a single biomarker
long_df <- analytic_df %>%
  select(cluster, t, all_of(biomarkers)) %>%
  pivot_longer(cols = -c(cluster, t), names_to = "biomarker", values_to = "value") %>%
  dplyr::filter(!is.na(value))  # Ensure to include only rows where biomarker values are not NA

# Counting the data points
count_data <- long_df %>%
  group_by(cluster, t, biomarker) %>%
  summarize(count = n(), .groups = "drop")

# Plotting
heatmap_plot <- ggplot(count_data, aes(y = biomarker, x = t, fill = count)) +
  geom_tile() +  # Use geom_tile for heatmap squares
  facet_wrap(~cluster, scales = "free_y") +  # Separate panel for each cluster
  theme_minimal() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 0, hjust = 1),  # Rotate x-axis text for readability
        axis.title = element_text(size = 12)) +
  labs(y = "Biomarker", x = "Years before diagnosis (t)", fill = "Count")

# Display the plot
print(heatmap_plot)
