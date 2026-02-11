# ============================================================================
# Simulation 2. 12.12.2025 Sample size 1000
# ============================================================================

## Setup ----
#setwd("~/Simulations/EUPROACT Simulations/proact-sim")
library(here)
library(arrow)
library(tidyverse)
library(furrr)
library(glue)
library(fs)
library(rms)
library(VGAM)
#library(devtools)
#install_github("scjohannes/markov.misc")
library(markov.misc)

## Configuration ----
N_PATIENTS <- 250000 # Size of superpopulation
N_SIMULATIONS <- 1000 # Number of simulation iterations
SAMPLE_SIZE <- c(1000) # Patients per sample. Run 1000 in a separate file.
N_BOOTSTRAP <- 1000 # Bootstrap iterations per sample     ### CURRENTLY NEEDS TO BE CHANGED MANUALLY IN THE FUNCTION SPECIFICATION (STEP 5)
N_CORES <- parallel::detectCores() - 1 # Parallel workers
SEED <- 123456


# File paths ----
OUTPUT_PATH <- here("scripts", "output", "sim301225")
DATA_PATH <- here("scripts", "data")
final_seed_path <- "seed_after_initial_run30_12_25_n1000.rds"
results_prefix <- "simulation_results_301225" # es_label and ss_label will be added

# Create directories if they don't exist
dir_create(OUTPUT_PATH)
dir_create(DATA_PATH)

# 0, 1.01, 1.50, 1.99 days
effect_sizes <- c(0, -0.007155, -0.0103, -0.01375)
  

for (i in effect_sizes) {
  ## ============================================================================
  ## STEP 1: Generate Superpopulation ----
  ## ============================================================================

  # Generate trajectories using Brownian motion model
  # This creates a 6-state model similar to VIOLET trial
  longitudinal_state_data <- sim_trajectories_brownian(
    n_patients = N_PATIENTS,
    follow_up_time = 60,
    treatment_prob = 0.5,
    absorbing_state = 6,
    seed = SEED,
    mu_treatment_effect = i
  )

  # plot sops
  plot_sops(longitudinal_state_data, geom = "line")

  # 1.28 more days at home over 60 days
  es <- calc_time_in_state_diff(longitudinal_state_data, target_state = 1)[1, 1]

  # Convert to different data format for other analyses

  t_data <- states_to_ttest(longitudinal_state_data, target_state = 1)

  drs_data <- states_to_drs(
    longitudinal_state_data,
    target_state = 1,
    follow_up_time = 60,
    covariates = NULL
  )

  count_data <- states_to_tte_old(
    longitudinal_state_data,
    covariates = NULL
    )

  # Prepare for proportional odds modeling
  # This removes absorbing states and converts variables to appropriate types
  markov_model_data <- prepare_markov_data(
    data = longitudinal_state_data,
    absorbing_state = 6,
    ordered_response = TRUE,
    factor_previous = TRUE # Convert yprev to factor
  )

  ## ============================================================================
  ## STEP 2: Save Superpopulations to Disk ----
  ## ============================================================================

  es_format <- sprintf("%.2f", round(es,2))
  
  # Full longitudinal states data
  long_states_path <- path(
    DATA_PATH,
    paste0("superpopulation_longitudinal_states_", es_format,
           ".parquet"))
  write_parquet(longitudinal_state_data, long_states_path)
  message(paste0("Saved to ", long_states_path))
  message(glue("File size: {file_size(long_states_path) |> as.character()}"))
  
  # Markov
  markov_path <- path(
    DATA_PATH,
    paste0("superpopulation_markov_", es_format, ".parquet")
  )
  write_parquet(markov_model_data, markov_path)
  message(paste0("Saved to ", markov_path))
  message(glue("File size: {file_size(markov_path) |> as.character()}"))
  
  # t-test
  t_path <- path(
    DATA_PATH,
    paste0("superpopulation_ttest_", es_format, ".parquet")
  )
  write_parquet(t_data, t_path)
  message(paste0("Saved to ", t_path))
  message(glue("File size: {file_size(t_path) |> as.character()}"))
  
  # DRS 60
  drs_path <- path(
    DATA_PATH,
    paste0("superpopulation_drs_", es_format, ".parquet")
  )
  write_parquet(drs_data, drs_path)
  message(paste0("Saved to ", drs_path))
  message(glue("File size: {file_size(drs_path) |> as.character()}"))
  
  # AG count data format
  count_path <- path(
    DATA_PATH,
    paste0("superpopulation_count_", es_format, ".parquet"))
  
  write_parquet(count_data, count_path)
  message(paste0("Saved to ", count_path))
  message(glue("File size: {file_size(count_path) |> as.character()}"))
  
}

## ============================================================================
## STEP 3: Identify Available Effect Sizes from Generated Files ----
## ============================================================================

# Get all markov files to determine available effect sizes
markov_files <- dir_ls(DATA_PATH, glob = "*superpopulation_markov_*.parquet")

# Extract effect sizes from filenames
available_effect_sizes <- markov_files |>
  path_file() |>
  str_extract("(?<=superpopulation_markov_)[0-9.-]+(?=\\.parquet)") |>
  as.numeric() |>
  sort()

message(glue(
  "Found {length(available_effect_sizes)} effect size scenarios: {paste(available_effect_sizes, collapse = ', ')}"
))

# Test only the required scenarios
available_effect_sizes <- available_effect_sizes[c(1, 3, 5, 7)]

## =============================================================================
## STEP 4: Define Analysis Functions ----
## =============================================================================

fit_markov_boot <- function(data, iter) {
  # Fit VGAM model
  fit <- VGAM::vglm(
    y ~ tx + rcs(time, 4) + factor(yprev),
    family = VGAM::cumulative(parallel = TRUE, reverse = TRUE),
    data = data
  )

  # Bootstrap within this function
  boot_coef <- bootstrap_model_coefs(
    fit,
    data = data,
    n_boot = 1000,
    id_var = "id"
  )
  
  # Tidy bootstrap results and add required columns
  tidy_coefs <- tidy_bootstrap_coefs(boot_coef, 
                                     probs = c(0.024, 0.976),
                                     estimate = "median") |>
    filter(term == "tx") |>
    mutate(
      iter = iter,
      analysis = "markov_coef",
      se_type = "boot",
      # Calculate bootstrap SE from confidence limits (approximate)
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      conf_low = lower,
      conf_high = upper,
      .before = 1
    ) |>
    select(iter, analysis, se_type, term, estimate, std_error,
           statistic, p_value, conf_low, conf_high)

  # Compute cluster robust standard errors
  robust_vcov_cl <- robcov_vglm(fit, cluster = data$id)
  
  robust_markov_ses <- tibble(
    iter = iter,
    analysis = "markov_coefs_rob",
    se_type = "robust",
    term = "tx",
    estimate = robust_vcov_cl$coefficients["tx"],
    std_error = robust_vcov_cl$se["tx"],
    statistic = robust_vcov_cl$z["tx"],
    p_value = robust_vcov_cl$pvalues["tx"],
    conf_low = robust_vcov_cl$coefficients["tx"] - 1.977 * robust_vcov_cl$se["tx"],
    conf_high = robust_vcov_cl$coefficients["tx"] + 1.977 * robust_vcov_cl$se["tx"]
    )
  
  out <- rbind(tidy_coefs, robust_markov_ses)
  
  # Also return all bootstrap iteration results to be saved and analysed later
  boot_results <- boot_coef
  naive_vglm <- VGAM::summary(fit)@coef3["tx", c("Estimate", "Pr(>|z|)")]

  return(list(out, list(boot_results, naive_vglm)))
}

# DRS analysis with rms::orm
fit_drs_orm <- function(data, iter) {
  # Set up datadist for rms
  # dd <- rms::datadist(data)
  # assign('dd', dd, envir = .GlobalEnv)
  # options(datadist = "dd", )

  # Fit ORM model
  fit <- rms::orm(drs ~ tx, data = data)

  # Tidy coefficients with one-sided test (greater (higher odds of better state))
  result <- tidy_po(fit, alternative = "two.sided", conf_level = 0.96) |>
    dplyr::filter(term == "tx") |>
    dplyr::mutate(
      iter = iter,
      analysis = "drs_orm",
      se_type = "robust",
      .before = 1
    )

  return(list(result))
}

# T-test analysis
fit_ttest <- function(data, iter) {
  # Run one-sided t-test
  # For some reason this tests HA that tx0 < tx1
  # Would have expected it to be the opposite and to have to use alternative = "greater"
  t_result <- t.test(
    y ~ tx,
    data = data,
    alternative = "two.sided",
    conf.level = 0.96
  )

  # Format results to match expected structure
  result <- tibble::tibble(
    iter = iter,
    analysis = "ttest",
    se_type = "naive",
    term = "tx",
    estimate = -diff(t_result$estimate), # tx1 - tx0
    std_error = t_result$stderr,
    statistic = t_result$statistic,
    p_value = t_result$p.value,
    conf_low = -t_result$conf.int[2],
    conf_high = -t_result$conf.int[1]
  )

  return(list(result))
}

# Time-to-event analysis: time-to-discharge and time-to-recovery
fit_time_to_event <- function(data, iter) {
  # Get rid of any time-series attributes
  data <- as.data.frame(unclass(as.list(data)))
  
  # Restructure data for time-to-recovery analysis -----------------------------
  data$lead_dur <- ave(data$stop-data$start, data$id, FUN = function(t) {
    c(t[-1], 1) }) # the interval length of the next interval
  data$recovery <- data$y == 1 & data$lead_dur > 3 # mark recovery (>3 days at home)
  split_data <- split(data, data$id)
  temp_ttr <- lapply(split_data, FUN = function(x) {
    # If there is one or more recoveries, select the first one...
    if (any(x$recovery == TRUE)) {
      x[x$recovery & x$stop == min(x$stop[x$recovery]),]
      # ... else select the last entry (death, or administrative censoring) 
    } else {x[nrow(x), ]}
  })
  data_ttr <- do.call(rbind, temp_ttr) # time to recovery per participant
  
  cox_model_ttr <- survival::coxph(
    survival::Surv(time = stop, event = recovery) ~ tx,
    #cluster = id,   # Check the type-I-error without robust SEs
    data = data_ttr
  )
  
  tidy_ttr <- tidy_po(cox_model_ttr, alternative = "two.sided") |>
    dplyr::filter(term == "tx") |>
    dplyr::mutate(
      iter = iter,
      analysis = "time_to_recovery",
      se_type = "naive",
      .before = 1
    )
  
  # Restructure data for time-to-discharge analysis ----------------------------
  temp_ttd <- lapply(split_data, FUN = function(x) {
    # If there is at least one discharge, select the first one...
    if (any(x$y == 1)) {
      x[x$y == 1 & x$stop == min(x$stop[x$y == 1]),]
      # ... else select the last entry (death, or administrative censoring) 
    } else {x[nrow(x), ]}
  })
  data_ttd <- do.call(rbind, temp_ttd) # time to recovery per participant
  
  # Creat event variable
  data_ttd$discharge <- data_ttd$y == 1
  
  cox_model_ttd <- survival::coxph(
    survival::Surv(time = stop, event = discharge) ~ tx,
    #cluster = id,   # Check the type-I-error without robust SEs
    data = data_ttd
  )
  
  tidy_ttd <- tidy_po(cox_model_ttd, alternative = "two.sided") |>
    dplyr::filter(term == "tx") |>
    dplyr::mutate(
      iter = iter,
      analysis = "time_to_discharge",
      se_type = "naive",
      .before = 1
    )
  
  return(list(rbind(
    tidy_ttr,
    tidy_ttd
  )))
}

# Test proportions of death using logistic regression
fit_binary_death <- function(data, iter) {
  # Collapse count data to one row per patient
  # We determine death if the patient EVER reached state 6 (death_status)
  analysis_data <- data |>
    group_by(id) |>
    summarise(
      died = as.integer(any(y == 6)),
      tx = first(tx),
      yprev = first(yprev)
    ) |>
    ungroup()
  
  # Fit Logistic Regression (GLM)
  fit <- stats::glm(
    died ~ tx + yprev,
    family = stats::binomial(link = "logit"),
    data = analysis_data
  )
  
  # Extract results
  # We rename columns to match the convention used in your fit_ttest
  bin_result <- tidy_po(fit, alternative = "two.sided") |>
    filter(term == "tx") |>
    mutate(
      iter = iter,
      analysis = "death_binary",
      se_type = "naive",
      .before = 1
    )
  
  return(list(bin_result))
}


# Define fitting functions (constant across all effect sizes)
fit_functions <- list(
  markov_boot = fit_markov_boot,
  drs_orm = fit_drs_orm,
  ttest = fit_ttest,
  time_to_event = fit_time_to_event,
  binary_death = fit_binary_death
)

## ============================================================================
## STEP 5: Run Power Simulation for Each Effect Size ----
## ============================================================================

# Set up parallel processing
library(furrr)
plan(future.callr::callr, workers = N_CORES)

# Set master seed
RNGkind("L'Ecuyer-CMRG")
set.seed(123456)

# Loop through each available effect size and sample size
for (j in seq_along(SAMPLE_SIZE)) {
  for (i in seq_along(available_effect_sizes)) {
    # Effect size
    effect_size <- available_effect_sizes[i]
    es_label <- sprintf("%.2f", effect_size)
    
    # Sample size
    sample_size <- SAMPLE_SIZE[j]
    ss_label <- sprintf("%.2f", sample_size)
    
    message(glue("Processing effect size scenario: {es_label}"))
    
    # Define data paths for this effect size
    data_paths <- list(
      markov_boot = path(
        DATA_PATH,
        paste0("superpopulation_markov_", es_label, ".parquet")
      ),
      drs_orm = path(
        DATA_PATH,
        paste0("superpopulation_drs_", es_label, ".parquet")
      ),
      ttest = path(
        DATA_PATH,
        paste0("superpopulation_ttest_", es_label, ".parquet")
      ),
      time_to_event = path(
        DATA_PATH, 
        paste0("superpopulation_count_", es_label, ".parquet")),
      binary_death = path(
        DATA_PATH, 
        paste0("superpopulation_count_", es_label, ".parquet"))
    )
    
    # Verify all files exist
    missing_files <- !file_exists(unlist(data_paths))
    if (any(missing_files)) {
      warning(glue(
        "Missing files for effect size {es_label}: {paste(names(data_paths)[missing_files], collapse = ', ')}"
      ))
      next
    }
  
    rerandomize <- effect_sizes[i] == 0
    message(paste0("Using rerandomize? ", rerandomize))
    
    # Set path to save details for each scenario
    DETAILS_PATH <- path(
      OUTPUT_PATH,
      paste0("scenario_", es_label, "_", ss_label)
    )
    
    start <- Sys.time() 
    # Run simulations in parallel
    message(glue("Starting {N_SIMULATIONS} simulation iterations..."))
    results <- future_map_dfr(
      1:N_SIMULATIONS,
      ~ assess_operating_characteristics(
        iter_num = .,
        data_paths = data_paths,
        output_path = DETAILS_PATH,
        fit_functions = fit_functions,
        sample_size = sample_size,
        allocation_ratio = 0.5,
        rerandomize = rerandomize,                           # ONLY TRUE FOR H0!
        seed = NULL,
      ),
      .options = furrr_options(
        seed = TRUE,
        packages = c("rms", "VGAM", "dplyr", "tibble", "markov.misc")
      ),
      .progress = TRUE
    )
    stop <- Sys.time()
    sim_duration <- stop-start
    print(sim_duration)
    
    # Add effect size identifier to results
    results <- results |>
      dplyr::mutate(scenario_label = paste("ES = ", es_label,
                                           "SIZE = ", ss_label), .before = 1)
    
    
    # # Summarize operating characteristics
    # summary_results <- summarize_oc_results(
    #   results,
    #   true_value = NULL, # Use the true effect from superpopulation
    #   alpha = 0.025 # One-sided alpha
    # ) |>
    #   dplyr::mutate(scenario_label = paste("ES = ", es_label,
    #                                        "SIZE = ", ss_label), .before = 1)
    # 
    # # Display summary
    # message(glue("\nSummary for effect size {es_label}:"))
    # print(summary_results)
    # SAVE RAW RESULTS INSTEAD:
    
    # Save results for this effect size and sample size
    # Set output path
    results_path <- path(
      OUTPUT_PATH,
      paste0(results_prefix, es_label, ss_label, ".parquet")
    )
    
    write_parquet(results, results_path)
    message(glue("Saved: {results_path}"))
    
    # summary_path <- path(
    #   OUTPUT_PATH,
    #   paste0("summary_results_", es_label, ".parquet")
    # )
    # write_parquet(summary_results, summary_path)
    # message(glue("Saved: {summary_path}"))
  }
}

final_seed <- .Random.seed
saveRDS(final_seed, final_seed_path)

# Add this at the beginning of the document in future extensions of these simulations
#.Random.seed <- readRDS("seed_after_initial_run12_12_25_n1000_es1.50.rds")

# Close parallel workers
plan("sequential")

## ============================================================================
## STEP 6: Load and Combine All Results ----
## ============================================================================
# Load all individual simulation results
results_files <- dir_ls(OUTPUT_PATH, glob = "*simulation_results_*.parquet")
all_results_combined <- map_dfr(results_files, read_parquet)

# # Load all individual summary results
# summary_files <- dir_ls(OUTPUT_PATH, glob = "*summary_results_*.parquet")
# all_summaries_combined <- map_dfr(summary_files, read_parquet)

# Save combined results
results_path_all <- path(OUTPUT_PATH, "simulation_results_all.parquet")
write_parquet(all_results_combined, results_path_all)

# summary_path_all <- path(OUTPUT_PATH, "summary_results_all.parquet")
# write_parquet(all_summaries_combined, summary_path_all)

# # Display final combined summary
# all_summaries_combined

# Compute summary over all results
