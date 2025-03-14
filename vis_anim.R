# reload data ----
results <- readRDS("./data/results.rds")
all_epoch_metrics <- readRDS("./data/all_epoch_metrics.rds")
all_train_history <- readRDS("./data/all_train_history.rds")
auc_results <- readRDS("./data/auc_results.rds")

# Animations
# gif 1 ----
library(ggplot2)
library(gganimate)
library(dplyr)

all_train_history <- all_train_history %>%
  group_by(noise_level, alpha, beta) %>%
  mutate(cum_avg_reward = cumsum(reward) / row_number())

p_anim <- ggplot(all_train_history, aes(x = sample, y = cum_avg_reward, color = factor(noise_level))) +
  geom_line(size = 1) +
  labs(title = 'Learning Process: Epoch {closest_state}', 
       x = 'Training Sample Index', 
       y = 'Cumulative Average Reward', 
       color = 'Noise Level') +
  transition_states(epoch, transition_length = 1, state_length = 0, wrap = FALSE) +
  ease_aes('cubic-in-out')

anim <- animate(p_anim, nframes = 200, fps = 20,  width = 600, height = 500, res = 120)
anim_save(filename = "./figures/gif_genetic_rl_learning.gif", animation = anim)

# gif 2 ----
# This one is not so informative. It shows color for pathogenicity but we don't see any learning-related change. 

# p_scatter <- ggplot(all_train_history, aes(x = GuRuScore, y = -log(PopulationFrequency), color = factor(ClinVarPathogenicity))) +
#   geom_point(size = 3, alpha = 0.7) +
#   labs(title = 'Learning Process on Variant Features: Epoch {closest_state}',
#        x = 'GuRu Score', 
#        y = 'Population Frequency', 
#        color = 'ClinVar Pathogenicity') +
#   transition_states(epoch, transition_length = 1, state_length = 0, wrap = FALSE) +
#   ease_aes('cubic-in-out')
# 
# anim_scatter <- animate(p_scatter, nframes = 200, fps = 20)
# anim_save("./figures/gif_genetic_rl_scatter.gif", anim_scatter)

# gif 3 ----
# Compute average p_action per epoch from the training history

# all_train_history_1 <- all_train_history |> filter(noise_level == 0.1)
# 
# avg_p <- all_train_history_1 %>%
#   group_by(epoch) %>%
#   summarize(avg_p_action = mean(p_action, na.rm = TRUE))
# 
# # Create the animated scatter plot with annotation at x = 0, y = -log(0.01)
# p_scatter_pact <- ggplot(all_train_history_1, aes(x = GuRuScore, y = -log(PopulationFrequency), color = p_action)) +
#   geom_point(size = 3, alpha = 0.7) +
#   scale_color_gradient(low = "blue", high = "red") +
#   labs(title = 'Evolving Predicition (p action)\non Variant Features: Epoch {closest_state}',
#        subtitle = 'Subset noise level 0.1',
#        x = 'GuRu Score', 
#        y = '-log(Population Frequency)', 
#        color = 'Predicted p_action') +
#   # Add annotation using the precomputed averages
#   geom_text(data = avg_p, 
#             aes(x = 5, y = -log(0.01), label = paste0("Avg p action:\n", sprintf("%.3f", avg_p_action)),
#             inherit.aes = FALSE, size = 6, color = "black") +
#   transition_states(epoch, transition_length = 1, state_length = 0, wrap = FALSE) +
#   ease_aes('cubic-in-out')

# Subset data for noise_level 0.1
all_train_history_1 <- all_train_history %>% filter(noise_level == 0.1)

# Compute average cumulative reward per epoch
avg_cum_reward <- all_train_history_1 %>%
  group_by(epoch) %>%
  summarize(avg_cum_reward = mean(cum_avg_reward, na.rm = TRUE))

# Create the animated scatter plot with annotation at x = 5, y = -log(0.01)
p_scatter_pact <- ggplot(all_train_history_1, aes(x = GuRuScore, y = -log(PopulationFrequency), color = p_action)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = 'Evolving Prediction (p_action) on Variant Features: Epoch {closest_state}',
       subtitle = 'Subset: noise level 0.1',
       x = 'GuRu Score', 
       y = '-log(Population Frequency)', 
       color = 'Predicted p action') +
  # Add annotation using the precomputed average cumulative reward
  geom_text(data = avg_cum_reward, 
            aes(x = 5, y = -log(0.01), label = paste0("Avg Cum Reward:\n", sprintf("%.3f", avg_cum_reward))),
            inherit.aes = FALSE, size = 6, color = "black") +
  transition_states(epoch, transition_length = 1, state_length = 0, wrap = FALSE) +
  ease_aes('cubic-in-out')

anim_scatter <- animate(p_scatter_pact, nframes = 200, fps = 20,  width = 800, height = 600, res = 120)
anim_save("./figures/gif_genetic_rl_scatter_pact.gif", anim_scatter)

# gif 5 ----
# Ensure cumulative average reward is computed
all_train_history <- all_train_history %>%
  group_by(noise_level, alpha, beta) %>%
  mutate(cum_avg_reward = cumsum(reward) / row_number())

# Animated scatter plot with GeneNumber on the x-axis
p_scatter_gene <- 
  # ggplot(all_train_history, aes(x = GeneNumber, y = cum_avg_reward, color = factor(noise_level))) +
  # geom_point(shape = 21, color = "black", size = 3, alpha = 0.6) +
  
  ggplot(all_train_history, aes(x = GeneNumber, y = cum_avg_reward)) +
  geom_jitter(shape = 21, color = "black", size = 3, aes(fill = factor(noise_level))) +
  
  labs(title = 'Learning Process: Epoch {closest_state}',
       x = 'Gene Number', 
       y = 'Cumulative Average Reward',
       fill = 'Noise Level') +
  transition_states(epoch, transition_length = 1, state_length = 0, wrap = FALSE) +
  ease_aes('cubic-in-out')

# anim_scatter_gene <- animate(p_scatter_gene, nframes = 200, fps = 20)
anim_scatter_gene <- animate(p_scatter_gene, nframes = 200, fps = 20,  width = 600, height = 500, res = 120)
anim_save("./figures/gif_genetic_rl_scatter_gene.gif", anim_scatter_gene)
