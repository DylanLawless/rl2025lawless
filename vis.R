# reload data ----
results <- readRDS("./data/results.rds")
all_epoch_metrics <- readRDS("./data/all_epoch_metrics.rds")
all_train_history <- readRDS("./data/all_train_history.rds")
auc_results <- readRDS("./data/auc_results.rds")

# Plot study data ----
# Generate combined data for all noise levels
df_list <- lapply(noise_levels, function(nl) {
  df_tmp <- generate_data(noise_level = nl)
  df_tmp$noise_level <- factor(nl)
  df_tmp
})
df_overlay <- do.call(rbind, df_list)

# Overlay Distribution Plots ----
# ACMG Score Distribution Overlay
p_acmg_overlay <- ggplot(df_overlay, aes(x = GuRuScore, fill = noise_level)) +
  geom_histogram(position = "identity", binwidth = 1, alpha = 0.4, color = "black") +
  labs(title = "ACMG Score Distribution", x = "ACMG Score", y = "Count", fill = "Noise Level") 

# Population Frequency Distribution Overlay
p_freq_overlay <- ggplot(df_overlay, aes(x = PopulationFrequency, fill = noise_level)) +
  geom_histogram(position = "identity", binwidth = 0.05, alpha = 0.4, color = "black") +
  labs(title = "Population Frequency Distribution", x = "Population Frequency", y = "Count", fill = "Noise Level") 

# Gene Number Distribution (discrete) - use position dodge
p_gene_overlay <- ggplot(df_overlay, aes(x = factor(GeneNumber), fill = noise_level)) +
  geom_bar(position = "dodge", color = "black") +
  labs(title = "Gene Number Distribution", x = "Gene Number", y = "Count", fill = "Noise Level") 

# ClinVar Pathogenicity Distribution for known entries
df_known_overlay <- df_overlay %>% filter(!is.na(ClinVarPathogenicity))

p_path_overlay <- ggplot(df_known_overlay, aes(x = factor(ClinVarPathogenicity), fill = noise_level)) +
  geom_bar(position = "dodge", color = "black") +
  labs(title = "ClinVar Pathogenicity (Known)", x = "Pathogenicity", y = "Count", fill = "Noise Level") 

# Combine overlay distribution plots in a grid
overlay_distributions <- (p_acmg_overlay | p_freq_overlay) / (p_gene_overlay | p_path_overlay)

# Cov corr ----
# Define mapping for variable names to new headings
var_labels <- c("GuRuScore" = "Guru Score", 
                "PopulationFrequency" = "Pop Freq", 
                "GeneNumber" = "Gene", 
                "ClinVarPathogenicity" = "Pathogenicity")

# Helper function for correlation plot per noise level with text labels
make_corr_plot <- function(nl) {
  df_nl <- df_overlay %>% filter(!is.na(ClinVarPathogenicity), noise_level == nl)
  num_vars <- c("GuRuScore", "PopulationFrequency", "GeneNumber", "ClinVarPathogenicity")
  cor_mat <- cor(df_nl[, num_vars])
  cor_df <- melt(cor_mat)
  ggplot(cor_df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_x_discrete(labels = var_labels) +
    scale_y_discrete(labels = var_labels) +
    labs(title = paste("Correlation (Noise =", nl, ")"), x = "", y = "", fill = "Corr") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Helper function for covariance plot per noise level with text labels
make_cov_plot <- function(nl) {
  df_nl <- df_overlay %>% filter(!is.na(ClinVarPathogenicity), noise_level == nl)
  num_vars <- c("GuRuScore", "PopulationFrequency", "GeneNumber", "ClinVarPathogenicity")
  cov_mat <- cov(df_nl[, num_vars])
  cov_df <- melt(cov_mat)
  ggplot(cov_df, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
    scale_fill_gradient(low = "white", high = "darkred") +
    scale_x_discrete(labels = var_labels) +
    scale_y_discrete(labels = var_labels) +
    labs(title = paste("Covariance (Noise =", nl, ")"), x = "", y = "", fill = "Cov") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Generate correlation and covariance plots for each noise level
corr_plots <- lapply(levels(df_overlay$noise_level), make_corr_plot)
cov_plots <- lapply(levels(df_overlay$noise_level), make_cov_plot)

# Stack correlation plots vertically and covariance plots vertically
corr_stack <- wrap_plots(corr_plots, ncol = 1)
cov_stack <- wrap_plots(cov_plots, ncol = 1)

# Combine correlation and covariance stacks side-by-side
matrix_plots <- corr_stack | cov_stack


# Combine All Plots and Save ---- 

overlay_distributions <- overlay_distributions + plot_annotation(title = "Data Distributions Across Noise Levels")

matrix_plots <- matrix_plots + plot_annotation(title = "Corr / Cov Matrices Across Noise Levels")

ggsave(filename = "./figures/data_dist.png", 
       plot = overlay_distributions, width = 8, height = 5)

ggsave(filename = "./figures/data_matrices.png", 
       plot = matrix_plots, width = 9, height = 8)

# Plotting aggregated results ----
# Epoch-level metrics: Average TD Error per Epoch, faceted by alpha, beta, colored by noise level
p_epoch_td <- ggplot(all_epoch_metrics, aes(x = epoch, y = avg_td_error, color = factor(noise_level))) +
  geom_line() +
  facet_grid(alpha ~ beta, labeller = label_both) +
  labs(title = "Avg TD Error per Epoch", x = "Epoch", y = "Avg TD Error", color = "Noise Level") 

# Average Reward per Epoch
p_epoch_reward <- ggplot(all_epoch_metrics, aes(x = epoch, y = avg_reward, color = factor(noise_level))) +
  geom_line() +
  facet_grid(alpha ~ beta, labeller = label_both) +
  labs(title = "Avg Reward per Epoch", x = "Epoch", y = "Avg Reward", color = "Noise Level") 

# AUC heatmap (faceted by noise level)
# p_auc <- ggplot(auc_results, aes(x = factor(alpha), y = factor(beta), fill = AUC)) +
#   geom_tile() +
#   facet_wrap(~ noise_level, labeller = label_both) +
#   scale_fill_gradient(low = "lightblue", high = "darkblue") +
#   labs(title = "AUC by Alpha and Beta\n(faceted by Noise Level)", x = "Alpha", y = "Beta")

p_auc <- ggplot(auc_results, aes(x = factor(alpha), y = factor(beta), fill = avg_auc)) +
  geom_tile() +
  facet_wrap(~ noise_level, labeller = label_both) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "AUC by Alpha and Beta\n(faceted by Noise Level)", x = "Alpha", y = "Beta")


# Learning curve: cumulative average reward over samples
all_train_history <- all_train_history %>%
  group_by(noise_level, alpha, beta) %>%
  mutate(cum_avg_reward = cumsum(reward)/row_number())

p_learning <- ggplot(all_train_history, aes(x = sample, y = cum_avg_reward, color = factor(noise_level))) +
  geom_line(size = 0.5, alpha = 0.6) +
  facet_grid(alpha ~ beta, labeller = label_both) +
  labs(title = "Cumulative Average Reward Over Training Samples", x = "Training Sample Index",
       y = "Cumulative Avg Reward", color = "Noise Level")

p_learning

# Calibration plot for test set predictions aggregated over combinations
calib_df <- data.frame()
# For each combination, regenerate test predictions and store them
for (i in 1:nrow(param_grid)) {
  curr_noise <- param_grid$noise_level[i]
  curr_alpha <- param_grid$alpha[i]
  curr_beta <- param_grid$beta[i]
  
  variant_data <- generate_data(curr_noise)
  known_data <- subset(variant_data, !is.na(ClinVarPathogenicity))
  train_indices <- sample(1:nrow(known_data), 0.7 * nrow(known_data))
  train_set <- known_data[train_indices, ]
  test_set <- known_data[-train_indices, ]
  
  n_epochs <- 20
  n_states <- 4 * 3 * 2
  w <- rep(0, n_states)
  v <- rep(0, n_states)
  for (epoch in 1:n_epochs) {
    train_set <- train_set[sample(nrow(train_set)), ]
    for (j in 1:nrow(train_set)) {
      row <- train_set[j, ]
      state_idx <- get_state(row$GuRuScore, row$PopulationFrequency, row$GeneNumber) + 1
      p_action <- sigmoid(w[state_idx])
      action <- ifelse(runif(1) < p_action, 1, 0)
      reward <- ifelse(action == row$ClinVarPathogenicity, 1, -1)
      td_error <- reward - v[state_idx]
      w[state_idx] <- w[state_idx] + curr_alpha * td_error * (action - p_action)
      v[state_idx] <- v[state_idx] + curr_beta * td_error
    }
  }
  test_set$Predicted <- NA
  action_probs <- numeric(nrow(test_set))
  for (k in 1:nrow(test_set)) {
    row <- test_set[k, ]
    state_idx <- get_state(row$GuRuScore, row$PopulationFrequency, row$GeneNumber) + 1
    p_action <- sigmoid(w[state_idx])
    action_probs[k] <- p_action
    test_set$Predicted[k] <- ifelse(p_action > 0.5, 1, 0)
  }
  temp_df <- data.frame(p_action = action_probs,
                        true_label = test_set$ClinVarPathogenicity,
                        noise_level = curr_noise,
                        alpha = curr_alpha,
                        beta = curr_beta)
  calib_df <- rbind(calib_df, temp_df)
}

calib_df <- calib_df %>%
  mutate(bin = cut(p_action, breaks = seq(0,1,by=0.1), include.lowest = TRUE))

calib_summary <- calib_df %>%
  group_by(noise_level, alpha, beta, bin) %>%
  summarize(mean_pred = mean(p_action),
            observed = mean(true_label),
            count = n(), .groups = "drop")

p_calib <- ggplot(calib_summary, aes(x = mean_pred, y = observed)) +
  geom_point(size = 2, color = "darkgreen") +
  geom_line(color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  facet_grid(alpha ~ beta + noise_level, labeller = label_both) +
  labs(title = "Calibration Plot by Parameter Combination",
       x = "Mean Predicted Probability", y = "Observed Proportion")

roc_data <- calib_df %>%
  group_by(noise_level, alpha, beta) %>%
  do({
    roc_obj <- roc(.$true_label, .$p_action, levels = c(0,1))
    data.frame(specificity = rev(roc_obj$specificities),
               sensitivity = rev(roc_obj$sensitivities))
  })

p_roc <- ggplot(roc_data, aes(x = specificity, y = sensitivity, 
                              group = interaction(noise_level, alpha, beta),
                              color = factor(noise_level))) +
  geom_line(size = 1) +
  facet_grid(alpha ~ beta, labeller = label_both) +
  geom_abline(intercept = 1, slope = -1, linetype = "dashed", color = "gray") +
  labs(title = "ROC Curves by Parameter Combination",
       x = "Specificity", y = "Sensitivity", color = "Noise Level") 


# Get All Plots ----
# Save each plot as a PNG file in the "./figures/" directory
if (!dir.exists("./figures/")) {
  dir.create("./figures/")
}

ggsave(filename = "./figures/epoch_td.png", plot = p_epoch_td, width = 8, height = 6)
ggsave(filename = "./figures/epoch_reward.png", plot = p_epoch_reward, width = 8, height = 6)
ggsave(filename = "./figures/auc.png", plot = p_auc, width = 8, height = 4)
ggsave(filename = "./figures/learning.png", plot = p_learning, width = 8, height = 6)
ggsave(filename = "./figures/calibration.png", plot = p_calib, width = 16, height = 6)
ggsave(filename = "./figures/roc_curve.png", plot = p_roc, width = 8, height = 6)
# ggsave(filename = "./figures/final_layout.png", plot = final_layout, width = 10, height = 10)


