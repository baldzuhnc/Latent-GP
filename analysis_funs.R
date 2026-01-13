# Functions for conducting post-fitting analyses

source("pipeline_funs.R")


# Posterior Predictive Check: compares observed data (bars) with model replications (dots)
ppc <- function(long_data, dim, yrep, observed_only = TRUE, n_pids = 20){
  
  if(!dim %in% names(yrep)){
    stop("Specified dimension not in data. Choose from 'pid', 'wave', 'item'")
  }
  
  if(dim == "pid"){
    selected_pids <- sample(unique(long_data$pid), n_pids)
    long_data <- long_data %>% filter(pid %in% selected_pids)
    yrep <- yrep %>% filter(pid %in% selected_pids)
  }

  p <- plot_data(long_data, dim = dim)
  
  # Optionally filter Y_rep to only observed cases
  if(observed_only) {
    # Create a flag for observed (non-missing) values in original data
    observed_cases <- long_data %>%
      drop_na(Y) %>%
      mutate(observed = TRUE) %>%
      dplyr::select(pid, wave, item, observed)
    
    # Filter Y_rep to only include cases that were observed in original data
    yrep <- yrep %>%
      left_join(observed_cases, by = c("pid", "wave", "item")) %>%
      filter(observed == TRUE) %>%
      dplyr::select(-observed)
  }
  
  # Compute counts
  yrep_summary_counts <- yrep %>%
    group_by(.draw, .data[[dim]], Y_rep) %>%
    summarise(count_rep = n(), .groups = "drop") %>%
    group_by(.data[[dim]], Y_rep) %>%
    summarise(
      lower = quantile(count_rep, 0.1),
      mean = mean(count_rep),
      upper = quantile(count_rep, 0.9),
      .groups = "drop"
    )

  return(p + geom_pointrange(data = yrep_summary_counts, 
                              aes(ymin=lower, ymax = upper, x = Y_rep, y = mean), 
                              color = "darkred"))
}


# Plots estimated IRT parameters (thresholds) for each item
plot_irt_params <- function(param, model_fit){
  
  # if(!param %in% c("item_thresh", "item_disc")){
  #  stop("Choose an IRT parameter. Either 'item_thresh' or 'item_disc'.")
  # }
  # y_var <- if(param == "item_thresh") "number" else "item"
  y_var <- "number"
  
  irt_params_estim <- model_fit %>% 
      # gather_rvars(item_thresh[item, number], item_disc[item]) %>% 
      gather_rvars(item_thresh[item, number]) %>% 
      mean_qi(.value) %>%
      dplyr::select(item, number, .variable, .value, .lower, .upper) %>%
      rename(mean_estim = .value, parameter = .variable) %>%
      filter(parameter == param)
  
  if(param == "item_thresh"){
    p <- ggplot() +
      ggstance::geom_pointrangeh(data = irt_params_estim, aes(xmin=.lower, xmax = .upper, y = factor(.data[[y_var]]) , x = mean_estim), color = "darkred") +
      geom_vline(xintercept = 0, color = "darkred") +
      facet_wrap(~ item, scales = "free_x") +
      labs(y = "", title = "Thresholds tau", x = "") +
      ggthemes::theme_few(base_size = 15) +
      theme(legend.position = "none")
  }

  # current model doesnt include discrimination parameter
  # if(param == "item_disc"){
  #  p <- ggplot() +
  #    ggstance::geom_pointrangeh(data = irt_params_estim, aes(xmin=.lower, xmax = .upper, y = factor(.data[[y_var]]) , x = mean_estim), color = "darkred") +
  #    geom_vline(xintercept = 1, color = "darkred") +
  #    labs(y = "Item j", title = "Discrimination alpha", x = "") +
  #    ggthemes::theme_few(base_size = 15) +
  #    theme(legend.position = "none")
  # }

  return(p)
}



# Plots posterior density of global GP parameters (length or amplitude)
plot_param_distr <- function(post_df, param = "length", bw = 4) {

  if (param == "length") {
    selected_param <- "gp_length_ind"
    plot_title <- "Posterior Density of GP Length Scale"
  } else if (param == "amplitude") {
    selected_param <- "gp_amplitude_ind"
    plot_title <- "Posterior Density of GP Amplitude"
  } else {
    stop("Unknown parameter. Use 'length' or 'amplitude'.")
  }
  
  post_df %>% 
    filter(.data$parameter == selected_param) %>%
    ggplot(aes(x = mean)) +
    geom_density(bw = bw, fill = "#3182bd", alpha = 0.6, color = "#08519c") +
    labs(
      title = plot_title,
      x = "Posterior Mean",
      y = "Density"
    ) +
    ggthemes::theme_few(base_size = 15) +
    theme(panel.grid.minor = element_blank())
}



# Plots estimated latent trajectories (theta) for all individuals
plot_theta_trajectory <- function(post_df, y_lim = NA, plotly = F, title = ""){
 pocket_df <- post_df %>% filter(parameter == "theta")
 pocket_p <- ggplot() +
    geom_line(data = pocket_df, aes(x = wave, y = mean, group = pid, text = paste("pid:", pid)), alpha = 0.1) +
    theme_minimal(base_size = 14) +
    labs(x = "Wave", y = "Theta (\u03b8)", title = title)

  if (!any(is.na(y_lim))) pocket_p <- pocket_p + ylim(y_lim)
  if(plotly) pocket_p <- plotly::ggplotly(pocket_p, tooltip = "text")

  return(pocket_p)
}


# Compares raw mean-level scores (standardized) with estimated latent theta for a sample of individuals
plot_ml_vs_theta <- function(ldata, post_df, n_pids = 20){
  ldata_nafree <- ldata %>% drop_na()

  # Subset to a random sample of PIDs for plotting
  selected_pids <- sample(unique(ldata_nafree$pid), n_pids)
  ldata_sub <- ldata_nafree %>% filter(pid %in% selected_pids)
  post_sub  <- post_df %>% filter(pid %in% selected_pids)

  p <- ggplot() +
    geom_line(data = ldata_sub %>% group_by(pid, wave) %>%
      summarise(latent = mean(Y), .groups = "drop") %>%
      mutate(latent = as.numeric(scale(latent))), 
      aes(x = factor(wave), y = latent, group = pid, color = factor(pid)),
      alpha = 0.5) +
    geom_point(data = ldata_sub %>% group_by(pid, wave) %>%
      summarise(latent = mean(Y), .groups = "drop") %>%
      mutate(latent = as.numeric(scale(latent))), 
      aes(x = factor(wave), y = latent, group = pid, color = factor(pid)),
      size = 0.5, alpha = 0.5) +
    geom_line(data = post_sub %>% filter(parameter == "theta"), 
      aes(x = factor(wave), y = mean, group = pid, color = factor(pid))) +
    facet_wrap(~pid) +
    geom_hline(yintercept = 0, color = "grey", alpha = 0.8) +
    ggthemes::theme_few(base_size = 15) +
    labs(x = "") +
    theme(legend.position = "none")

  return(p)
}



# Filters out individuals with extreme theta values
remove_outliers <- function(post_df, theta_lower = -10, theta_upper = 10) {
  before <- length(unique(post_df$pid))
  
  valid_pids <- post_df %>%
    filter(parameter == "theta") %>%
    group_by(pid) %>%
    summarise(valid = max(mean) < theta_upper & min(mean) > theta_lower, .groups = "drop") %>%
    filter(valid) %>%
    pull(pid)
  
  removed <- before - length(valid_pids)
  if(removed > 0){
    cat(sprintf("Note: Removed %d trajectories with extreme values. %d remain.\n", removed, length(valid_pids)))
  }

  return(post_df %>% filter(pid %in% valid_pids))
}


