library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(caret)
library(pROC)
library(patchwork)
library(doParallel)
library(foreach)
library(reshape2)
library(knitr)

# Define parameter grids ----
noise_levels <- c(0.1, 0.2, 0.3)
alphas <- c(0.1, 0.2, 0.4)
betas <- c(0.05, 0.10, 0.25)

# Create a grid of all parameter combinations
param_grid <- expand.grid(noise_level = noise_levels,
                          alpha = alphas,
                          beta = betas,
                          stringsAsFactors = FALSE)

# Overall data generation parameters
set.seed(666)  # For reproducibility of data generation
total_entries <- 2000 # start time 22:22
known_ratio <- 0.5          # 50% known, 50% unknown
known_entries <- total_entries * known_ratio
unknown_entries <- total_entries - known_entries

# data ----
# Define gene pools and groups
non_pathogenic_genes <- 1:6   # More likely to have non-pathogenic variants
pathogenic_genes <- 4:10      # More likely to have pathogenic variants

groups <- list(
  group1 = list(pathogenicity = 0, acmg_range = 0:5,  freq_range = c(0.5, 1.0), gene_pool = non_pathogenic_genes),
  group2 = list(pathogenicity = 0, acmg_range = 6:10, freq_range = c(0.3, 0.7), gene_pool = non_pathogenic_genes),
  group3 = list(pathogenicity = 1, acmg_range = 11:15, freq_range = c(0.1, 0.5), gene_pool = pathogenic_genes),
  group4 = list(pathogenicity = 1, acmg_range = 16:20, freq_range = c(0.0, 0.3), gene_pool = pathogenic_genes)
)

n_per_group <- as.integer(known_entries / length(groups))  # force integer

# Function to generate synthetic data given a noise level
generate_data <- function(noise_level) {
  variant_data <- data.frame(
    VariantNumber = integer(),
    GeneNumber = integer(),
    ClinVarPathogenicity = integer(),
    GuRuScore = integer(),
    PopulationFrequency = numeric(),
    stringsAsFactors = FALSE
  )
  # Generate known data
  for (i in seq_along(groups)) {
    group <- groups[[i]]
    known_part <- data.frame(
      VariantNumber = seq(from = nrow(variant_data) + 1, length.out = n_per_group),
      GeneNumber = sample(group$gene_pool, n_per_group, replace = TRUE),
      ClinVarPathogenicity = rep(group$pathogenicity, n_per_group),
      GuRuScore = sample(group$acmg_range, n_per_group, replace = TRUE),
      PopulationFrequency = runif(n_per_group, group$freq_range[1], group$freq_range[2])
    )
    # Introduce noise: flip a fraction of the labels
    flip_indices <- sample(1:nrow(known_part), size = round(nrow(known_part) * noise_level))
    known_part$ClinVarPathogenicity[flip_indices] <- 1 - known_part$ClinVarPathogenicity[flip_indices]
    variant_data <- rbind(variant_data, known_part)
  }
  # Generate unknown data (set ClinVarPathogenicity as NA of integer type)
  unknown_data <- data.frame(
    VariantNumber = seq(from = nrow(variant_data) + 1, to = nrow(variant_data) + unknown_entries),
    GeneNumber = sample(1:10, unknown_entries, replace = TRUE),
    ClinVarPathogenicity = rep(NA_integer_, unknown_entries),
    GuRuScore = sample(0:20, unknown_entries, replace = TRUE),
    PopulationFrequency = runif(unknown_entries, 0, 1)
  )
  variant_data <- rbind(variant_data, unknown_data)
  rownames(variant_data) <- NULL
  return(variant_data)
}

# Define helper functions ----
sigmoid <- function(x) 1/(1+exp(-x))
get_state <- function(acmg, freq, gene) {
  acmg_bin <- if (acmg <= 5) 0 else if (acmg <= 10) 1 else if (acmg <= 15) 2 else 3
  freq_bin <- if (freq >= 0.5) 0 else if (freq >= 0.3) 1 else 2
  pathogenic_genes <- 4:10
  gene_risk <- ifelse(gene %in% pathogenic_genes, 1, 0)
  state <- acmg_bin * 6 + freq_bin * 2 + gene_risk
  return(state)
}

# Setup parallel backend ----
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Parallel loop over parameter combinations
results <- foreach(i = 1:nrow(param_grid), .combine = 'rbind', .packages = c("dplyr", "caret", "pROC")) %dopar% {
  curr_noise <- param_grid$noise_level[i]
  curr_alpha <- param_grid$alpha[i]
  curr_beta <- param_grid$beta[i]
  
  # Generate data with current noise level
  variant_data <- generate_data(curr_noise)
  
  # Split known data into training and test sets
  known_data <- subset(variant_data, !is.na(ClinVarPathogenicity))
  train_indices <- sample(1:nrow(known_data), 0.7 * nrow(known_data))
  train_set <- known_data[train_indices, ]
  test_set <- known_data[-train_indices, ]
  
  # # Initialise RL parameters and weights
  n_epochs <- 20
  n_states <- 4 * 3 * 2
  w <- rep(0, n_states)
  v <- rep(0, n_states)
  
  epoch_metrics <- data.frame()
  train_history <- data.frame()
  sample_counter <- 0
  
  # log_file <- sprintf("./figures/experiment_log_worker_%d.txt", Sys.getpid())
  # Create a log file unique to the current worker based on its process ID.
  log_file <- sprintf("./log/experiment_log_worker_%d.txt", Sys.getpid())
  
  for (epoch in 1:n_epochs) {
    print(paste("epoch:", epoch))
    td_errors <- numeric(nrow(train_set))
    rewards <- numeric(nrow(train_set))
    train_set <- train_set[sample(nrow(train_set)), ]
    for (j in 1:nrow(train_set)) {
      sample_counter <- sample_counter + 1
      row <- train_set[j, ]
      state_idx <- get_state(row$GuRuScore, row$PopulationFrequency, row$GeneNumber) + 1
      p_action <- sigmoid(w[state_idx])
      action <- ifelse(runif(1) < p_action, 1, 0)
      reward <- ifelse(action == row$ClinVarPathogenicity, 1, -1)
      td_error <- reward - v[state_idx]
      w[state_idx] <- w[state_idx] + curr_alpha * td_error * (action - p_action)
      v[state_idx] <- v[state_idx] + curr_beta * td_error
      td_errors[j] <- td_error
      rewards[j] <- reward
      train_history <- rbind(train_history, data.frame(epoch = epoch,
                                                       sample = j,
                                                       state = state_idx,
                                                       w_val = w[state_idx],
                                                       v_val = v[state_idx],
                                                       p_action = p_action,
                                                       reward = reward,
                                                       td_error = td_error,
                                                       noise_level = curr_noise,
                                                       alpha = curr_alpha,
                                                       beta = curr_beta))
    }
    epoch_metrics <- rbind(epoch_metrics, data.frame(epoch = epoch,
                                                     avg_td_error = mean(td_errors),
                                                     avg_reward = mean(rewards),
                                                     noise_level = curr_noise,
                                                     alpha = curr_alpha,
                                                     beta = curr_beta))
    
    # At the end of each epoch, log the epoch details:
    # cat(sprintf("Running: noise=%.2f, alpha=%.2f, beta=%.2f, Epoch %d \n",
                # Sys.getpid(), curr_noise, curr_alpha, curr_beta, epoch),
        # file = log_file, append = TRUE)
    
 
    # Inside the loop over epochs, log with the epoch counter:
    # log_message <- sprintf("Worker %d: Running combination: noise_level=%.2f, alpha=%.2f, beta=%.2f, Epoch %d\n", 
                           # Sys.getpid(), curr_noise, curr_alpha, curr_beta, epoch)
    # cat(log_message, file = log_file, append = TRUE)
    
    log_message <- sprintf("Worker %d: Running combination: noise_level=%.2f, alpha=%.2f, beta=%.2f, Epoch %d\n", 
                           Sys.getpid(), curr_noise, curr_alpha, curr_beta, epoch)
    cat(log_message, file = log_file, append = TRUE)
    
  }
  
  # Evaluate on test set using trained weights
  action_probs <- numeric(nrow(test_set))
  for (k in 1:nrow(test_set)) {
    row <- test_set[k, ]
    state_idx <- get_state(row$GuRuScore, row$PopulationFrequency, row$GeneNumber) + 1
    p_action <- sigmoid(w[state_idx])
    action_probs[k] <- p_action
    test_set$Predicted[k] <- ifelse(p_action > 0.5, 1, 0)
  }
  
  roc_obj <- roc(response = test_set$ClinVarPathogenicity,
                 predictor = action_probs,
                 levels = c(0,1))
  auc_val <- auc(roc_obj)
  
  # Combine results into one data frame (adding a column to distinguish parameter combo)
  epoch_metrics$combination <- paste(curr_noise, curr_alpha, curr_beta, sep = "_")
  train_history$combination <- paste(curr_noise, curr_alpha, curr_beta, sep = "_")
  data.frame(noise_level = curr_noise,
             alpha = curr_alpha,
             beta = curr_beta,
             avg_auc = auc_val,
             epoch_metrics = I(list(epoch_metrics)),
             train_history = I(list(train_history)))
  
}




# Stop the cluster when done ----
stopCluster(cl)

print("if problem do: kill -9 $(ps aux | grep  R.framework | awk '{print $2}')")

# Combine aggregated results from each parameter combination
# Here we extract the nested data frames and then bind rows together
all_epoch_metrics <- do.call(rbind, results$epoch_metrics)
all_train_history <- do.call(rbind, results$train_history)
auc_results <- results %>% dplyr::select(noise_level, alpha, beta, avg_auc)

# (The rest of your plotting code remains unchanged.)



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

# # Bin the training history by grouping sample indices into bins (e.g. every 100 samples)
# library(dplyr)
# binned_history <- all_train_history %>%
#   mutate(bin = sample %/% 100) %>%  # Create bins (change 100 to desired bin width)
#   group_by(noise_level, alpha, beta, bin) %>%
#   summarise(mean_cum_avg = mean(cum_avg_reward), .groups = "drop")
# 
# # Plot the binned cumulative average reward
# p_learning <- ggplot(binned_history, aes(x = bin * 100, y = mean_cum_avg, color = factor(noise_level))) +
#   geom_line(size = 1) +
#   facet_grid(alpha ~ beta, labeller = label_both) +
#   labs(title = "Binned Cumulative Average Reward Over Training Samples", 
#        x = "Training Sample Index", 
#        y = "Mean Cumulative Average Reward", 
#        color = "Noise Level")
# p_learning


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

# Ensure the figures directory exists
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



# Plot data ----

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

