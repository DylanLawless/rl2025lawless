# RL2025Lawless: An Actor-Critic Reinforcement Learning Framework for Genetic Variant Pathogenicity

## Overview

This project implements a reinforcement learning (RL) framework designed to estimate the probability of observing genetic variants in disease. Instead of directly predicting variant pathogenicity, the framework quantifies the cumulative evidence supporting a variant’s clinical observability within a Bayesian context. The approach integrates established genomic metrics—including the GuRu score (formerly known as the ACMGuru score), gene risk priors, and population frequency—to learn which evidence correctly supports known variant labels. This adaptive learning step is intended to form the basis for subsequent Bayesian integration, yielding nuanced probability estimates for variant-disease associations.

## Key Features

- **Reinforcement Learning Framework:** Uses an actor-critic algorithm to learn from simulated genomic data with variable label noise.
- **Bayesian Evidence Integration:** Focuses on estimating the probability of observing a variant in disease, rather than a direct pathogenicity prediction.
- **Simulated Data Generation:** Creates synthetic datasets with 2,000 variants (50% known, 50% unknown) and varying noise levels (10%, 20%, 30%) to mimic real-world genomic variability.
- **Performance Evaluation:** Assesses model performance using ROC curves, AUC, calibration plots, cumulative learning curves, and temporal-difference error metrics.
- **Parallel Computing:** Employs `doParallel` and `foreach` for efficient hyperparameter grid search over noise levels and learning rates.

## Dependencies

- **R** (tested on version X.X.X)
- R packages: `ggplot2`, `dplyr`, `caret`, `pROC`, `patchwork`, `doParallel`, `foreach`, `reshape2`, `knitr`

## Usage

1. **Run Simulations:** Execute `sim.R` or `sim_loop.R` to generate synthetic genomic data, train the RL model across various parameter settings, and evaluate performance.
2. **View Figures:** Output plots (e.g., ROC curves, calibration plots, learning curves) are saved in the `figures/` directory.
3. **Manuscript Preparation:** LaTeX source files in the `latex/` directory are used to compile the final manuscript PDF.



