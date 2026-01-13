# ========== Modeling attitudes in SOSEC via latent Gaussian Processes ==========
source("pipeline_funs.R")

# 1. Configuration & Parameters ####
# This script can be run via Rscript: Rscript pipeline.R <country> <varset> <waves_from> <waves_to>
# Or interactively (uses defaults below)

# Define available varsets
varsets <- list(
  trust      = c("F5A2_1", "F5A3_1", "F5A4_1", "F5A5_1", "F5A6_2", "F5A7_1"),
  discontent = c('F1A10_1', 'F1A16_1', 'F1A17_1', 'F1A19_1'),
  efficacy   = c('F1A5_2', 'F1A10_2'),
  anger      = c('F1A9_1', 'F1A15_1'),
  conspiracy = c('F1A16_2', 'F1A17_2', 'F5cA3_1')
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 4) {
  country         <- args[1]
  varset_name     <- args[2]
  waves_from      <- as.integer(args[3])
  waves_to        <- as.integer(args[4])
  consec_missings <- ifelse(length(args) >= 5, as.integer(args[5]), 2) #Optional: Set consecutive missings in command arguments, default: 2
  length_scale    <- ifelse(length(args) >= 6, as.numeric(args[6]), 4.0) #Optional: Set consecutive missings in command arguments, default: 4
} else {
  # Default settings for interactive use
  # Default settings for interactive use
  # Default settings for interactive use
  country         <- "de"
  varset_name     <- "anger"
  waves_from      <- 55
  waves_to        <- 65
  consec_missings <- 2
  length_scale    <- 4.0
}

if (!varset_name %in% names(varsets)) {
  stop(sprintf("Unknown varset '%s'. Available: %s", varset_name, paste(names(varsets), collapse = ", ")))
}

vars <- list(title = varset_name, vars = varsets[[varset_name]])

data <- get_data(country, waves_from, waves_to, consec_missings, vars$vars, vars$title, length_scale)

# Wrangle data ####

# Subset variables and indices, remove consecutive missings
subsetted_data <- subset_data(data)

if (nrow(subsetted_data) == 0) {
  stop("Subsetting returned no data. Check waves or ID selection.")
}

cat("Observations per wave:\n")
print(table(subsetted_data$wave))
cat("\nUnique individuals in timeframe: ",
    length(unique(subsetted_data$pid)), "\n")
cat("====================================\n\n")


# Tidyformat: Make long and fill missing indices ####
long_data <- subset_to_long_df(subsetted_data)

# Plot data ####
plot_data(long_data, dim = "wave")
plot_data(long_data, dim = "pid", n_pids = 40) #n_pids necessary for computational efficiency
plot_data(long_data, dim = "item")
plot_ml_trajectory(long_data) #scaled ML trajectories per individual

# Prepare for stan; Extract stan data and participant/time mapping ####
stan_ready <- prep_for_stan(long_data, length_scale)
model_data <- stan_ready$model_data

pid_map <- stan_ready$pid_map
wave_map <- stan_ready$wave_map



# ========== Run model & save ============
fit <- stan_sample(model_data, show_exceptions = FALSE)
save_fit(fit, title = data$title)

cat("\n============ Diagnostic summary ===============\n")
print(fit$diagnostic_summary())
cat("===============================================\n")

# ========== Extract & save ============
cat("\nExtracting and saving results...\n")
# Data
save_data(long_data, pid_map, wave_map, data$title)

#GP summary
gp_summary <- extract_posterior_summary(fit, pid_map, wave_map)
save_gp_summary(gp_summary, data$title)

#Yrep
yrep <- extract_yrep_sample(fit, pid_map, wave_map, n_samples = 100)
save_yrep(yrep, data$title)
## ==========================================