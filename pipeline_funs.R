#Functions for streamlining the latent variable pipeline

#0. Paths and parameters ####
model_path <- "model/grm_gp.stan"
gp_summary_path <- "gp_summary/"
yrep_path <- "yrep/"
fitted_path <- "fitted/"
data_path <- "subsetted_data/"

# Ensure output directories exist (paths defined in pipeline_funs.R)
for (path in c(gp_summary_path, yrep_path, fitted_path, data_path)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

#1. Imports ####
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(tidyverse)
p_load(cmdstanr)
p_load(tidybayes)
p_load(posterior)
p_load(peakRAM)

#3. Initialization of stan model ####
options(mc.cores = parallel::detectCores())
model <- cmdstan_model(model_path)

#4. Pipeline functions ####
get_latest <- function(country){
  # Return: path to latest wave file for `country`. Errors if none found.
  files <- list.files(paste0("data/", country), pattern = paste0("^W\\d+\\_",country, ".rds$"), full.names = TRUE)
  if(length(files) == 0){
    stop(paste0("No data files found for country: ", country))
  }
  wave_numbers <- as.numeric(gsub("\\D", "", basename(files)))
  return(files[which.max(wave_numbers)])
}

check_vars <- function(data, vars) {
  # Check if all items have the same scale length
  max_values <- data$df |>
    summarise(across(all_of(vars$vars), ~ max(.x, na.rm = TRUE))) |>
    as_vector()
  
  if (length(unique(max_values)) != 1) {
    message("Items have different scale lengths. Stopping.")
    cat("Max values per item:\n")
    print(max_values)
    stop("Scale length mismatch detected")
  }
}

get_data <- function(country, waves_from, waves_to, consec_missings, vars, title = "", length_scale){
  # Short: build a data descriptor list for later pipeline steps.
  # Args: country, waves_from/waves_to (integers), consec_missings (int), vars (character vector)
  # Returns: list with `title`, `df`, `vars`, `waves`, `missings_allowed`.
  
  data <- readRDS(get_latest(country)) %>% 
    dplyr::select(wave, pid, all_of(vars))

  check_vars(data, vars)
  
  title <- paste0(
    country, "_", 
    title, "_", 
    waves_from, ":",
    waves_to, "_", 
    "nas:", consec_missings, 
    "_length:", length_scale, "_",
    format(Sys.time(), "%a %b %d %X"))

  
  ret <- list(
    title = title, 
    df = data, 
    vars = vars, 
    waves = c(waves_from:waves_to), 
    missings_allowed = consec_missings)
  
  return(ret)
}

subset_data <- function(data, cutoff_ids = NULL){
	
  vars <- data$vars
  waves <- data$waves
  missings_allowed <- data$missings_allowed
  datadf <- data$df

  if(is.null(cutoff_ids)){
    # select all if null
    cutoff_ids <- unique(datadf$pid)
  }else{
    cat("Subsetting participants. ")
  }
	
	#1. Subset cols, waves and ids
	sdata <- datadf %>%
    filter(wave %in% waves, pid %in% cutoff_ids) %>%
    mutate(across(-pid, as.integer)) #make all vars integers
  
  #2. Subset data that contains more than consecutive missings
  cat("\n====================================\n")
  cat("Excluding observations with more than", missings_allowed , "missings in a row.")
  cat("\n====================================\n")

	all_combinations <- expand_grid(
  	pid = unique(sdata$pid),
    wave = unique(sdata$wave))
  
  # Join with existing data
  compl_data <- right_join(sdata, all_combinations, by = c("wave", "pid")) %>% 
  	arrange(pid, wave)
  
  # Check for consecutive NAs using run length encoding
  # Determine whether an observation (pid x wave) is missing across ALL items
  consec_missing <- compl_data %>%
    group_by(pid) %>%
    #mutate(is_missing = ifelse(wave %in% data$wave & pid %in% data$pid, F, T)) %>% #doesnt work but would be more robust
  	mutate(is_missing = is.na(get(vars[1]))) %>%
    mutate(run_id = cumsum(c(0, diff(is_missing) != 0))) %>%
    group_by(pid, run_id, is_missing) %>%
    mutate(run_length = n()) %>%
    ungroup() %>%
    filter(is_missing == TRUE & run_length > missings_allowed) %>% 
    distinct(pid) %>%
    pull(pid)
  
  rdata <- sdata %>%
  	filter(!pid %in% consec_missing)
    
  return(rdata)
  
}

subset_to_long_df <- function(sdata){
  # Short: convert wide participant-wave-item data to long format
  # Args: sdata with columns `pid`, `wave`, item columns...
  # Returns: long dataframe with columns pid, wave, item, Y
	all_combinations <- expand_grid(
		pid = unique(sdata$pid),
    item = 1:(ncol(sdata)-2),
    wave = unique(sdata$wave))
	
	ldata <- sdata %>%
    pivot_longer(cols = -c(pid, wave), names_to = "item", values_to = "Y") %>%
    mutate(item = as.integer(factor(item, levels = colnames(sdata)[3:ncol(sdata)]))) %>%
    right_join(all_combinations, by = c("pid", "wave", "item"))
	
	return(ldata)
}

prep_for_stan <- function(long_data, length_scale){
  # Short: prepare data structures and maps for Stan model input
  # Args: long_data (pid, wave, item, Y), length_scale (numeric)
  # Returns: list with `model_data`, `pid_map`, `wave_map` ready for Stan
	
	pid_map <- tibble(
		pid = unique(long_data$pid),
		i = seq_along(unique(long_data$pid)))
	
	wave_map <- tibble(
		wave = unique(long_data$wave),
		t = seq_along(unique(long_data$wave)))
	
	stan_data <- long_data %>% 
		left_join(pid_map, by = "pid") %>%
		left_join(wave_map, by = "wave") %>%
    dplyr::select(-c("pid", "wave")) %>%
		mutate_all(as.integer) %>% 
		drop_na()
	
	model_data <- list(
		K = max(unique(stan_data$Y)), #K categories
  	N = length(stan_data$Y),      #N total observations
  	P = max(unique(stan_data$item)), #P items
  	tmax = max(unique(stan_data$t)), #length of timescale
  	I = length(unique(stan_data$i)), #I participants
  	participant = stan_data$i,
  	item = stan_data$item,
  	time = stan_data$t,
  	Y = stan_data$Y,
    length_scale = length_scale)
  
	return(list("model_data" = model_data, 
							"pid_map" = pid_map, 
							"wave_map" = wave_map))
}

stan_sample <- function(model_data, show_exceptions = T){
  # Short: run CmdStan sampling on prepared `model_data`
  # Returns: CmdStanMCMC fit object
	
  if(!show_exceptions){
    cat("Caution: Not printing informational messages from stan!")
  }
  
  model_fit <- model$sample(data = model_data, 
														chains = 4,
                            adapt_delta = 0.95,
														parallel_chains = 4,
                            show_exceptions = show_exceptions,
														refresh = 500,
														save_warmup = TRUE)

	return(model_fit)
}

save_fit <- function(fit, title){
  # Short: persist a CmdStan fit object to `fitted_path` using `qs::qsave`
  # Args: fit (CmdStan fit), title (string)
	# Load CmdStan output files into the fitted model object.
	fit$draws() # Load posterior draws into the object.
	try(fit$sampler_diagnostics(), silent = TRUE) # Load sampler diagnostics.
	try(fit$init(), silent = TRUE) # Load user-defined initial values.
	try(fit$profiles(), silent = TRUE) # Load profiling samples.

  path <- paste0(fitted_path, title, ".rds")

	# Save the object to a file.
	qs::qsave(x = fit, file = path)
}




#5. Posterior ####


extract_posterior_summary <- function(model_fit, pid_map, wave_map) {
  # Short: extract posterior summaries for theta and GP amplitude
  # Returns tidy dataframe of posterior means and 95% intervals joined to pid/wave
  #draws <- model_fit$draws(variable = c("theta", "gp_length_ind", "gp_amplitude_ind"))
  #draws_raw <- model_fit$draws(variable = "theta")
  mem_usage2 <- peakRAM(
    draws <- as_draws_df(model_fit),
      # Get summaries for theta and amplitude and length
    
      summary <- draws %>%
        #gather_rvars(theta[i, t], gp_length_ind[i], gp_amplitude_ind[i]) %>%
        gather_rvars(theta[i, t], gp_amplitude_ind[i]) %>%
        mean_qi(.value) %>%
        left_join(pid_map, "i") %>%
        left_join(wave_map, "t") %>%
        rename(parameter = .variable, mean = .value, lower95 = .lower, upper95 = .upper) %>%
        relocate(wave, pid, i, t, parameter) %>%
        dplyr::select(-c(.width, .point, .interval))

  )

  cat("\n====================================\n")
  cat("Peak RAM during posterior processing:",mem_usage2$Peak_RAM_Used_MiB[1],"MB\n")
  cat("====================================\n")



  return(summary)
}


extract_posterior_theta <- function(model_fit, n_draws = 200) {
  # Short: extract a sample of theta draws
  # Returns: tibble of theta draws (i,t,.value)
	theta4_estim_raw  <- model_fit %>% 
		spread_draws(theta[i,t], ndraws = n_draws) %>%
		mutate(i = as.factor(i),
					 t = as.factor(t))
  return(theta4_estim_raw)
}

extract_yrep_sample <- function(model_fit, pid_map, wave_map, n_samples = 100) {
  # Short: extract a small sample of Y_rep draws for diagnostics
  # Only extract a small sample of Y_rep for diagnostics
  draws <- model_fit$draws(variables = "Y_rep", inc_warmup = FALSE)
  
  # Get Y_rep draws but sample only n_samples
  yrep_draws <- draws %>%
    posterior::subset_draws(variable = "Y_rep", 
                           draw = sample(posterior::ndraws(draws), 
                           min(n_samples, posterior::ndraws(draws)))) %>%
  	spread_draws(Y_rep[i,j,t]) %>%
    left_join(pid_map, by = "i") %>%
    left_join(wave_map, by = "t") %>%
    ungroup() %>%
    rename("item" = "j") %>%
    dplyr::select(-c(i,t))
  
  return(yrep_draws)
}



save_gp_summary <- function(post_summary_df, title){
  # Short: save posterior summary dataframe to `gp_summary_path`
  cat(paste0("Saving GP posterior summary dataframe to ", gp_summary_path, title, ".rds", "\n"))
	write_rds(post_summary_df, file = paste0(gp_summary_path, title, ".rds"))
}

save_yrep <- function(yrep, title){
  # Short: save `yrep` draws to `yrep_path`
  cat(paste0("Saving yrep to ", yrep_path, title, ".rds", "\n"))
	write_rds(yrep, file = paste0(yrep_path, title, ".rds"))
}

save_data <- function(long_data, pid_map, wave_map, title){
  # Short: save processed long data and maps to `data_path`
  tosave <- list(long_data = long_data, pid_map = pid_map, wave_map = wave_map)
  
  cat(paste0("Saving data to ", data_path, title, ".rds", "\n"))
	write_rds(tosave, file = paste0(data_path, title, ".rds"))
}


# Descriptive plotting ####

plot_data <- function(ldata, dim, n_pids = NULL){

  # Short: simple diagnostic plot of counts by category across `dim`
  
  data <- ldata %>% drop_na() 

  if(!dim %in% names(data)){
    stop("Specified dimension not in data. Choose from 'pid', 'wave', 'item'")
  }

  if(!is.null(n_pids)){
    subset_ids <- sample(unique(data$pid), 20)
    data <- data %>% filter(pid %in% subset_ids)
  }

  title <- case_when(
    dim == "wave" ~ "Wave t",
    dim == "item" ~ "Item j",
    dim == "pid" ~ "Individual i"
  )

  ggplot(data = data %>% group_by(.data[[dim]], Y) %>% summarize(Y_count = n())) +
    geom_col(aes(x = Y, y =  Y_count, fill = factor(Y))) +
    facet_wrap(~.data[[dim]]) +
    ggsci::scale_fill_jco() +
    ggthemes::theme_few(base_size = 15) + 
    labs(title = paste(title)) +
    theme(legend.position = "none")
  
}

plot_ml_trajectory <- function(ldata){
  # Short: plot mean (ML) trajectory per individual across waves

  data <- ldata %>% rename("t" = "wave", "i" = "pid") %>% drop_na() 

  p <- ggplot() +
    geom_line(data = data %>% group_by(i, t) %>%
      summarise(ml = mean(Y), .groups = "drop") %>%
      mutate(ml = as.numeric(scale(ml))), 
      aes(x = t, y = ml, group = i, color = factor(i)),
      alpha = 0.8,
      linetype = "dashed") +
    facet_wrap(~i) +
    labs(y = "Value", title = "ML") +
    theme(legend.position = "none")
    
  return(p)
  
}




