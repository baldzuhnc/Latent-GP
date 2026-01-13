![plots/cover.png]

# Latent GPIRT: Gaussian Process Item Response Theory

This project implements a Gaussian Process Latent Trajectory Model within an Item Response Theory (IRT) framework. It is designed to estimate individual-level trajectories of attitudes or traits over multiple survey waves.

## Workflow Overview

### 1. Data Preparation
- Place your raw survey data (CSV format) into the `data/` folder.
- The pipeline currently expects `overall_extended_germany_data.csv` or `overall_extended_us_data.csv`.

### 2. Standardize Data (`1_dataloader.R`)
- Run this script to clean the raw data, standardize ID columns to `pid`, and save the data in a versioned `.rds` format.
- Output is stored in `data/de/` or `data/us/`.

### 3. Fit the Model (`2_pipeline.R`)
This script runs the Stan model. You can use it in two ways:
- **Interactive:** Open the script and set variables like `country`, `varset_name`, and `waves_from:waves_to` at the top.
- **CLI (Batch):** Run from the terminal for automated processing:
  ```bash
  Rscript 2_pipeline.R de trust 55 65
  ```
- **Outputs:** Stan fits are saved in `fitted/`, and summaries are saved in `gp_summary/`.

### 4. Analyze Results (`3_analysis.R`)
- Use this interactive script to visualize results.
- **PPC:** Posterior Predictive Checks to validate model fit.
- **Trajectories:** Visualize latent trait ($\theta$) change over time for individuals or the population.
- **IRT Params:** Inspect item thresholds ($\tau$) and GP parameters (length-scale, amplitude).

## Directory Structure

- `data/`: Raw CSV files and processed `.rds` data.
- `model/`: Stan source code (`grm_gp.stan`).
- `fitted/`: Saved model objects (serialized with `qs`).
- `gp_summary/`: Extracted posterior summaries for fast plotting.
- `plots/`: Exported visualizations.
- `pipeline_funs.R`: Core utility functions for data prep and model interaction.
- `analysis_funs.R`: Visualization and diagnostic functions.

## Requirements
- R with `tidyverse`, `cmdstanr`, `qs`, and `tidybayes`.
- Recent version of CmdStan (installed via `cmdstanr::install_cmdstan()`).
