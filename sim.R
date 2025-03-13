
set.seed(666)  # For reproducibility

# Parameter to control total sample size and noise level
total_entries <- 1000
known_ratio <- 0.5          # 50% known, 50% unknown
known_entries <- total_entries * known_ratio
unknown_entries <- total_entries - known_entries
noise_level <- 0.2          # 20% of known labels will be flipped

# Initialize an empty data frame
variant_data <- data.frame(
  VariantNumber = integer(),
  GeneNumber = integer(),
  ClinVarPathogenicity = integer(),
  ACMGScore = integer(),
  PopulationFrequency = numeric()
)

# Define gene pools: higher number genes are more likely to be pathogenic
non_pathogenic_genes <- 1:6   # More likely to have non-pathogenic variants
pathogenic_genes <- 4:10      # More likely to have pathogenic variants

# Define groups with conditions (4 groups)
groups <- list(
  group1 = list(pathogenicity = 0, acmg_range = 0:5,  freq_range = c(0.5, 1.0), gene_pool = non_pathogenic_genes),
  group2 = list(pathogenicity = 0, acmg_range = 6:10, freq_range = c(0.3, 0.7), gene_pool = non_pathogenic_genes),
  group3 = list(pathogenicity = 1, acmg_range = 11:15, freq_range = c(0.1, 0.5), gene_pool = pathogenic_genes),
  group4 = list(pathogenicity = 1, acmg_range = 16:20, freq_range = c(0.0, 0.3), gene_pool = pathogenic_genes)
)

# Known entries: evenly split among groups
n_per_group <- known_entries / length(groups)

# Generate known data
for (i in seq_along(groups)) {
  group <- groups[[i]]
  known_part <- data.frame(
    VariantNumber = seq(from = nrow(variant_data) + 1, length.out = n_per_group),
    GeneNumber = sample(group$gene_pool, n_per_group, replace = TRUE),
    ClinVarPathogenicity = rep(group$pathogenicity, n_per_group),
    ACMGScore = sample(group$acmg_range, n_per_group, replace = TRUE),
    PopulationFrequency = runif(n_per_group, group$freq_range[1], group$freq_range[2])
  )
  # Introduce noise: flip a fraction of the labels
  flip_indices <- sample(1:nrow(known_part), size = round(nrow(known_part) * noise_level))
  known_part$ClinVarPathogenicity[flip_indices] <- 1 - known_part$ClinVarPathogenicity[flip_indices]
  variant_data <- rbind(variant_data, known_part)
}

# Generate unknown data
unknown_data <- data.frame(
  VariantNumber = seq(from = nrow(variant_data) + 1, to = nrow(variant_data) + unknown_entries),
  GeneNumber = sample(1:10, unknown_entries, replace = TRUE),
  ClinVarPathogenicity = rep(NA, unknown_entries),
  ACMGScore = sample(0:20, unknown_entries, replace = TRUE),
  PopulationFrequency = runif(unknown_entries, 0, 1)
)

# Combine known and unknown datasets
variant_data <- rbind(variant_data, unknown_data)
rownames(variant_data) <- NULL

# Display the first few rows
head(variant_data)




set.seed(123)
library(ggplot2)
library(dplyr)
library(caret)
library(patchwork)

# Split known data into training and test sets
known_data <- subset(variant_data, !is.na(ClinVarPathogenicity))
train_indices <- sample(1:nrow(known_data), 0.7 * nrow(known_data))
train_set <- known_data[train_indices, ]
test_set <- known_data[-train_indices, ]

# Initialise key parameters and reinforcement learning (RL) weights

alpha <- 0.2       # Actor learning rate: controls how quickly the policy (action weights) is updated
beta <- 0.10       # Critic learning rate: controls how quickly the value function (critic weights) is updated
n_epochs <- 20     # Number of training epochs: each epoch represents one complete pass through the training data
n_states <- 4 * 3 * 2  # Total number of discrete states: calculated from 4 bins for ACMGScore, 3 bins for PopulationFrequency, and 2 levels for Gene risk, resulting in 24 unique states
w <- rep(0, n_states)  # Initialize the actor (policy) weights as a vector of zeros for each state
v <- rep(0, n_states)  # Initialise the critic (value function) weights as a vector of zeros for each state

# Data frames to log metrics
epoch_metrics <- data.frame(epoch = integer(), avg_td_error = numeric(), avg_reward = numeric())
train_history <- data.frame(epoch = integer(), sample = integer(), state = integer(),
                            w_val = numeric(), v_val = numeric(), p_action = numeric(),
                            reward = numeric(), td_error = numeric())

# Training loop with detailed logging
sample_counter <- 0
for (epoch in 1:n_epochs) {
  td_errors <- numeric(nrow(train_set))
  rewards <- numeric(nrow(train_set))
  # Shuffle training set each epoch
  train_set <- train_set[sample(nrow(train_set)), ]
  for (i in 1:nrow(train_set)) {
    sample_counter <- sample_counter + 1
    row <- train_set[i, ]
    state_idx <- get_state(row$ACMGScore, row$PopulationFrequency, row$GeneNumber) + 1
    p_action <- sigmoid(w[state_idx])
    action <- ifelse(runif(1) < p_action, 1, 0)
    true_label <- row$ClinVarPathogenicity
    reward <- if (action == true_label) 1 else -1
    td_error <- reward - v[state_idx]
    # Update weights
    w[state_idx] <- w[state_idx] + alpha * td_error * (action - p_action)
    v[state_idx] <- v[state_idx] + beta * td_error
    td_errors[i] <- td_error
    rewards[i] <- reward
    # Log sample metrics
    train_history <- rbind(train_history, data.frame(epoch = epoch,
                                                     sample = i,
                                                     state = state_idx,
                                                     w_val = w[state_idx],
                                                     v_val = v[state_idx],
                                                     p_action = p_action,
                                                     reward = reward,
                                                     td_error = td_error))
    # Print key iteration information every 1000 samples
    if (sample_counter %% 1000 == 0) {
      cat(sprintf("Epoch %d, Sample %d, State %d, p_action: %.3f, Reward: %d, TD error: %.3f\n",
                  epoch, sample_counter, state_idx, p_action, reward, td_error))
    }
  }
  epoch_metrics <- rbind(epoch_metrics, data.frame(epoch = epoch,
                                                   avg_td_error = mean(td_errors),
                                                   avg_reward = mean(rewards)))
  cat(sprintf("Epoch %d complete: Avg TD error: %.3f, Avg reward: %.3f\n",
              epoch, mean(td_errors), mean(rewards)))
}


# Plot 1: Average TD Error per Epoch
p1 <- ggplot(epoch_metrics, aes(x = epoch, y = avg_td_error)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Average TD Error per Epoch", x = "Epoch", y = "Avg TD Error") +
  theme_minimal()
p1

# Plot 2: Average Reward per Epoch
p2 <- ggplot(epoch_metrics, aes(x = epoch, y = avg_reward)) +
  geom_line(color = "red", size = 1) +
  labs(title = "Average Reward per Epoch", x = "Epoch", y = "Avg Reward") +
  theme_minimal()

# Plot 3: Evolution of Actor Weight for a Representative State (state index 1)
# rep_state_history <- train_history %>% filter(state == 1)
# Use dplyr::filter explicitly to avoid conflict with ensembldb::filter
rep_state_history <- dplyr::filter(train_history, state == 1)

p3 <- ggplot(rep_state_history, aes(x = 1:nrow(rep_state_history), y = w_val)) +
  geom_line(color = "darkgreen", size = 0.8) +
  labs(title = "Evolution of Actor Weight (w) for State 1", x = "Sample Count", y = "w") +
  theme_minimal()

# Plot 4: Distribution of Action Probabilities over all training samples
p4 <- ggplot(train_history, aes(x = p_action)) +
  geom_histogram(bins = 20, fill = "purple", color = "black") +
  labs(title = "Distribution of Action Probabilities", x = "P(action=1)", y = "Count") +
  theme_minimal()

final_metric_plot <- (p1 + p2) / (p3 + p4)
print(final_metric_plot)

# Evaluate the trained model on the test set
test_set$Predicted <- NA
action_probs <- numeric(nrow(test_set))
for (i in 1:nrow(test_set)) {
  row <- test_set[i, ]
  state_idx <- get_state(row$ACMGScore, row$PopulationFrequency, row$GeneNumber) + 1
  p_action <- sigmoid(w[state_idx])
  action_probs[i] <- p_action
  test_set$Predicted[i] <- ifelse(p_action > 0.5, 1, 0)
}

conf_mat <- confusionMatrix(as.factor(test_set$Predicted), as.factor(test_set$ClinVarPathogenicity))
print(conf_mat)

conf_df <- as.data.frame(as.table(conf_mat$table))
p_conf <- ggplot(conf_df, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq), color = "white", size = 6) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Confusion Matrix", x = "Predicted", y = "Reference") +
  theme_minimal()

# Plot distribution of action probabilities for test examples
p_test <- ggplot(data.frame(p_action = action_probs), aes(x = p_action)) +
  geom_histogram(bins = 20, fill = "orange", color = "black") +
  labs(title = "Test Set:\nDistribution of Action Probabilities", x = "P(action=1)", y = "Count") +
  theme_minimal()

print(p_test)




final_metric_plot <- (p1 + p2) / (p3 + p4) + (p_conf + p_test)
print(final_metric_plot)













# measure ----

library(pROC)
library(caret)
library(dplyr)
library(ggplot2)

#---------------------------
# ROC Curve and AUC
#---------------------------
roc_obj <- roc(response = test_set$ClinVarPathogenicity,
               predictor = action_probs,
               levels = c(0,1))
auc_val <- auc(roc_obj)
cat(sprintf("AUC: %.3f\n", auc_val))

p_roc <- ggplot(data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)),
  aes(x = specificity, y = sensitivity)) +
  geom_line(color = "darkred", size = 1) +
  geom_abline(intercept = 1, slope = -1, linetype = "dashed") +
  labs(title = sprintf("ROC Curve (AUC = %.3f)", auc_val),
       x = "Specificity", y = "Sensitivity") +
  theme_minimal()

print(p_roc)

#---------------------------
# Precision, Recall, F1 Score
#---------------------------
cm <- confusionMatrix(as.factor(test_set$Predicted), 
                      as.factor(test_set$ClinVarPathogenicity), 
                      positive = "1")
print(cm)

precision <- cm$byClass["Pos Pred Value"]
recall <- cm$byClass["Sensitivity"]
f1 <- 2 * precision * recall / (precision + recall)
cat(sprintf("Precision: %.3f, Recall: %.3f, F1 Score: %.3f\n", precision, recall, f1))

#---------------------------
# Learning Curves Over Time
# (Epoch-level metrics already plotted: avg TD error and avg reward)
# Additionally, plot cumulative average reward over samples during training
train_history$cum_avg_reward <- cumsum(train_history$reward) / 1:nrow(train_history)
p_learning <- ggplot(train_history, aes(x = 1:nrow(train_history), y = cum_avg_reward)) +
  geom_line(color = "steelblue", size = 1) +
  labs(title = "Cumulative Average Reward Over Training Samples",
       x = "Training Sample Index", y = "Cumulative Average Reward") +
  theme_minimal()
print(p_learning)

#---------------------------
# Calibration Plot
# Bin predicted probabilities and compare observed frequencies
calib_df <- data.frame(p_action = action_probs,
                       true_label = test_set$ClinVarPathogenicity)
# Create bins
calib_df <- calib_df %>%
  mutate(bin = cut(p_action, breaks = seq(0,1,by=0.1), include.lowest = TRUE))
calib_summary <- calib_df %>%
  group_by(bin) %>%
  summarize(mean_pred = mean(p_action),
            observed = mean(true_label),
            count = n())

p_calib <- ggplot(calib_summary, aes(x = mean_pred, y = observed)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_line(color = "darkgreen") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(title = "Calibration Plot",
       x = "Mean Predicted Probability", y = "Observed Proportion") +
  theme_minimal()
print(p_calib)

#---------------------------
# Combine All Plots into One Final Layout
#---------------------------
final_metric_plot2   <- (p_roc / p_learning / p_calib)
print(final_metric_plot2)
final_metric_plot

