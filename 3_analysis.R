# ========== Gaussian Process IRT Post-Fitting Analysis ==========
source("analysis_funs.R")

# 1. Selection & Loading ####

# Available results in fitted folder:
# list.files(fitted_path)

# Manual selection or pick latest
title <- "de_trust_55:65_nas:2_length:4_Fri Oct 31 14:11:25.rds" 

if (!file.exists(file.path(fitted_path, title))) {
  stop("Selected fit file not found: ", title)
}

message("Loading results for: ", title)

# Load all components
fit        <- qs::qread(file.path(fitted_path, title))
gp_summary <- readRDS(file.path(gp_summary_path, title))
yrep       <- readRDS(file.path(yrep_path, title))

# Data and mappings
data_bundle <- readRDS(file.path(data_path, title))
long_data   <- data_bundle$long_data
pid_map     <- data_bundle$pid_map
wave_map    <- data_bundle$wave_map

# 2. Posterior Predictive Checks (PPC) ####
# Compare observed distributions (bars) with model replications (dots)
ppc(long_data, dim = "wave", yrep = yrep, observed_only = TRUE)
ppc(long_data, dim = "pid",  yrep = yrep, observed_only = TRUE, n_pids = 16)
ppc(long_data, dim = "item", yrep = yrep, observed_only = TRUE)

# 3. Parameter Distributions ####

# Item Thresholds (Tau)
plot_irt_params(param = "item_thresh", model_fit = fit)

# Global GP Parameters (Length Scale / Amplitude)
# plot_param_distr(gp_summary, param = "length", bw = 0.5)
plot_param_distr(gp_summary, param = "amplitude", bw = 0.1)

# Diagnostic: Posterior scatter (check for correlations or funnels)
# bayesplot::mcmc_scatter(fit$draws(), pars = c("item_thresh[1,3]", "gp_amplitude_ind[8]"))

# 4. Latent Trajectories ####

# All-individual trend (Theta)
plot_theta_trajectory(gp_summary, plotly = FALSE, title = "Population Latent Trajectories")

# Individual comparison: Observed Mean Level (ML) vs. Estimated Theta
plot_ml_vs_theta(ldata = long_data, post_df = gp_summary, n_pids = 12)

# 5. Outlier Filtering ####
# Create a subset by removing extreme trajectories if necessary
# sub_gpsummary <- remove_outliers(gp_summary, theta_lower = -5, theta_upper = 5)
