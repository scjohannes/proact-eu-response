# ============================================================================
# Simulation 1. 05.12.2025
# ============================================================================

## Setup ----
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
SAMPLE_SIZE <- c(250, 500, 700) # Patients per sample
N_BOOTSTRAP <- 500 # Bootstrap iterations per sample     ### CURRENTLY NEEDS TO BE CHANGED MANUALLY IN THE FUNCTION SPECIFICATION (STEP 5)
N_CORES <- parallel::detectCores() - 1 # Parallel workers
SEED <- 558712

OUTPUT_PATH <- here("scripts", "output")
DATA_PATH <- here("scripts", "data")

# Create directories if they don't exist
dir_create(OUTPUT_PATH)
dir_create(DATA_PATH)

# 0 Days, 0.6 Days,
effect_sizes <- c(0, 
                  -0.004, -0.008, 
                  -0.012
)


for (i in effect_sizes) {
  ## ============================================================================
  ## STEP 1: Generate Superpopulation ----
  ## ============================================================================

  # Generate trajectories using Brownian motion model with default parameters
  # This creates a 6-state model similar to VIOLET trial
  # Default: no treatment effect (mu_treatment_effect = 0)
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
    markov_data, 
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

  # Full logitudinal states data 
  long_states_path <- path(
    DATA_PATH,
    paste0("superpopulation_longitudinal_states", "_", round(es, digits = 2), 
           ".parquet"))
  write_parquet(longitudinal_state_data, long_states_path)
  message(paste0("Saved to ", long_states_path))
  message(glue("File size: {file_size(long_states_path) |> as.character()}"))
  
  # Markov
  markov_path <- path(
    DATA_PATH,
    paste0("superpopulation_markov", "_", round(es, digits = 2), ".parquet")
  )
  write_parquet(markov_model_data, markov_path)
  message(paste0("Saved to ", markov_path))
  message(glue("File size: {file_size(markov_path) |> as.character()}"))

  # t-test
  t_path <- path(
    DATA_PATH,
    paste0("superpopulation_ttest", "_", round(es, digits = 2), ".parquet")
  )
  write_parquet(t_data, t_path)
  message(paste0("Saved to ", t_path))
  message(glue("File size: {file_size(t_path) |> as.character()}"))

  # DRS 60
  drs_path <- path(
    DATA_PATH,
    paste0("superpopulation_drs", "_", round(es, digits = 2), ".parquet")
  )
  write_parquet(drs_data, drs_path)
  message(paste0("Saved to ", drs_path))
  message(glue("File size: {file_size(drs_path) |> as.character()}"))
  
  # AG count data format
  count_path <- path(
    DATA_PATH, 
    paste0("superpopulation_count_", round(es, digits = 2), ".parquet"))
  
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

# Test only H0 and strong effect
#available_effect_sizes <- available_effect_sizes

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
  boot_output <- bootstrap_standardized_sops(
    fit,
    data = data,
    n_boot = 500,
    parallel = FALSE
  )
  
  # Compute time in target state by treatment
  boot_tis <- time_in_state(boot_output[["sops"]], target_states = 1)
  median_tis <- median(boot_tis$delta)
  lower_tis <- quantile(boot_tis$delta, 0.025)
  upper_tis <- quantile(boot_tis$delta, 0.975)
  boot_sops <- as.data.frame(
    list(estimate = median_tis, 
         lower = lower_tis, 
         upper = upper_tis)
    ) |> mutate(term = "tx", .before = 1)
  
  # Get bootstrapped coefficients
  boot_coefs <- tidy_bootstrap_coefs(boot_output[["coefs"]], 
                                     probs = c(0.025, 0.975), 
                                     estimate = "median") |>
    filter(term == "tx") |> 
    select(!n_iter)
  
  # Rowbind to get one row per method and add required columns 
  out <- rbind(boot_sops, boot_coefs) |>
    mutate(
      iter = rep(iter, 2),
      analysis = c("markov_taooh", "markov_coefs"),
      se_type = rep("boot", 2),
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
  
  return(out)
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
  result <- tidy_po(fit, alternative = "two.sided", conf_level = 0.95) |>
    dplyr::filter(term == "tx") |>
    dplyr::mutate(
      iter = iter,
      analysis = "drs_orm",
      se_type = "robust",
      .before = 1
    )

  return(result)
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
    conf.level = 0.95
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

  return(result)
}

# Define fitting functions (constant across all effect sizes)
fit_functions <- list(
  markov_taooh = fit_markov_boot,
  drs_orm = fit_drs_orm,
  ttest = fit_ttest
)

## ============================================================================
## STEP 5: Run Power Simulation for Each Effect Size ----
## ============================================================================

# Set up parallel processing
library(furrr)
plan(future.callr::callr, workers = N_CORES)

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
      markov_taooh = path(
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
      )
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
    
    start <- Sys.time() 
    # Run simulations in parallel
    message(glue("Starting {N_SIMULATIONS} simulation iterations..."))
    results <- future_map_dfr(
      1:N_SIMULATIONS,
      ~ assess_operating_characteristics(
        iter_num = .,
        data_paths = data_paths,
        fit_functions = fit_functions,
        sample_size = sample_size,
        allocation_ratio = 0.5,
        rerandomize = rerandomize,               # ONLY TRUE FOR H0!
        seed = SEED,
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
    
    # SAVE RAW RESULTS INSTEAD
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
    
    # Save results for this effect size
    results_path <- path(
      OUTPUT_PATH,
      paste0("simulation_results_051225", es_label, ss_label, ".parquet")
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
