# Load libraries
library(ggplot2)
library(dplyr)

# ---- Step 1: Visualize the data ----
plot_data <- df_joined %>% sample_n(500000)   # adjust size as needed

ggplot(plot_data, aes(x = stat_group, y = mean_dist, fill = stat_group)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribution of mean_dist by group")

# ---- Step 2: Observed difference in means ----
obs_diff <- plot_data %>%
  group_by(stat_group) %>%
  summarise(mean_val = mean(mean_dist, na.rm = TRUE)) %>%
  summarise(diff = diff(mean_val)) %>%  # high - low
  pull(diff)

obs_diff

# ---- Step 3: Permutation test ----
set.seed(123)   # for reproducibility
n_perm <- 10000
perm_diffs <- replicate(n_perm, {
  shuffled <- sample(plot_data$stat_group)  # shuffle labels
  plot_data %>%
    mutate(shuffled_group = shuffled) %>%
    group_by(shuffled_group) %>%
    summarise(mean_val = mean(mean_dist, na.rm = TRUE)) %>%
    summarise(diff = diff(mean_val)) %>%
    pull(diff)
})

# ---- Step 4: Two-sided p-value ----
p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
p_value <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)
p_value

# ---- Step 5: Visualize null distribution ----
data.frame(perm_diffs = perm_diffs) %>%
  ggplot(aes(x = perm_diffs)) +
  geom_histogram(bins = 50, fill = "grey70", color = "black") +
  geom_vline(xintercept = obs_diff, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = -obs_diff, color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = "Permutation test: Null distribution of mean differences",
       x = "Difference in mean(mean_dist) (high - low)",
       y = "Frequency")


group_means <- df_joined %>%
  group_by(stat_group) %>%
  summarise(mean_val = mean(mean_dist, na.rm = TRUE))

group_means
