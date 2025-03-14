library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(caret)
library(pROC)
library(patchwork)
library(doParallel)
library(foreach)
library(reshape2)
library(knitr)
library(gganimate)
library(dplyr)

# In this version we update the train_history to save features for interpretation in annimations 

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
# total_entries <- 2000 # start time 13:50 end 14:03 (13 minutes)
total_entries <- 3000 # start time 14:10 end 14:35 (25 minutes)
# total_entries <- 100
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
      
      train_history <- rbind(train_history, data.frame(
        epoch = epoch,
        sample = j,
        state = state_idx,
        GeneNumber = row$GeneNumber,   # Include GeneNumber here
        w_val = w[state_idx],
        v_val = v[state_idx],
        p_action = p_action,
        reward = reward,
        td_error = td_error,
        noise_level = curr_noise,
        alpha = curr_alpha,
        beta = curr_beta,
        GuRuScore = row$GuRuScore,
        PopulationFrequency = row$PopulationFrequency,
        ClinVarPathogenicity = row$ClinVarPathogenicity
      ))
      
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

# save data ----
# After generating results, save them to an RDS file
saveRDS(results, file = "./data/results.rds")
saveRDS(all_epoch_metrics, file = "./data/all_epoch_metrics.rds")
saveRDS(all_train_history, file = "./data/all_train_history.rds")
saveRDS(auc_results, file = "./data/auc_results.rds")


# reload data ----
results <- readRDS("./data/results.rds")
all_epoch_metrics <- readRDS("./data/all_epoch_metrics.rds")
all_train_history <- readRDS("./data/all_train_history.rds")
auc_results <- readRDS("./data/auc_results.rds")
