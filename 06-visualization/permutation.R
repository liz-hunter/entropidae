library(tidyverse)

# use mean
read_entropy_summary <- df %>%
  group_by(SequenceID, taxid, name, cv_group) %>%
  summarise(mean_entropy = mean(entropy), .groups = "drop")

# use standard deviation
read_entropy_summary <- df %>%
  group_by(SequenceID, taxid, name, cv_group) %>%
  summarise(sd_entropy = sd(entropy), .groups = "drop")

# plot 
ggplot(read_entropy_summary, aes(x = cv_group, y = mean_entropy)) +  # or sd_entropy
  geom_violin(trim = FALSE, fill = "gray80") +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  theme_minimal() +
  labs(title = "Entropy scores by CV group",
       y = "Mean entropy per read",
       x = "Taxa CV group")

# test for significance (mean or sd)
read_entropy_summary %>%
  wilcox_test(mean_entropy ~ cv_group) %>%
  add_significance()

