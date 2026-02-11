#setwd()
library(arrow)
library(here)
library(purrr)
library(dplyr)
library(markov.misc)
library(ggplot2)
library(fs)

# Set paths
OUTPUT_PATH1 <- here("..", "proact-sim", "scripts", "output", "sim121225")
OUTPUT_PATH2 <- here("..", "proact-sim", "scripts", "output", "sim121225_n1000")

# Get all results that are present
results_files1 <- dir_ls(OUTPUT_PATH1, glob = "*simulation_*results_*.parquet")
results_files2 <- dir_ls(OUTPUT_PATH2, glob = "*simulation_*results_*.parquet")
results_files <- c(results_files1, results_files2)
results_files <- results_files[!grepl("all", results_files)]
all_results_combined <- map_dfr(results_files, read_parquet)

# Save combined results
results_path_all <- path(OUTPUT_PATH1, "simulation_results_all.parquet")
write_parquet(all_results_combined, results_path_all)

# Check detailed output
details_taooh_iter3 <- readRDS("scripts/output/sim121225/scenario_1.01_600.00/details/markov_taooh/markov_taooh_iter_3.rds")
details_drs_iter3 <- readRDS("scripts/output/sim121225/scenario_1.01_600.00/details/drs_orm/drs_orm_iter_3.rds")

# Get summary 
summary_HA <- all_results_combined |>
  mutate(decision = case_when(
    analysis == "markov_coefs" ~ estimate < 0 & conf_high < 0,
    analysis == "markov_taooh" ~ estimate > 0 & conf_low > 0,
    analysis == "ttest" ~ p_value < 0.05 & estimate < 0,
    analysis == "drs_orm" ~ p_value < 0.05 & estimate > 0,
    analysis == "time_to_discharge" ~ p_value < 0.05 & estimate > 0,
    analysis == "time_to_recovery" ~ p_value < 0.05 & estimate > 0,
    analysis == "death_binary" ~ p_value < 0.05 & estimate < 0
  )) |>
  group_by(scenario_label, analysis) |>
  summarise(power = mean(decision)) |>
  mutate(
    sample_size = as.numeric(vapply(strsplit(
      scenario_label, split = "SIZE = ", fixed = TRUE), 
      \(x) x[2],
      FUN.VALUE = character(1))),
    effect_size = as.numeric(vapply(strsplit(
      substr(scenario_label, start = 7, stop = 100), split = " SIZE = ", fixed = TRUE), 
      \(x) x[1],
      FUN.VALUE = character(1))),
    analysis = factor(analysis),
    MCSE = sqrt(power*(1-power)/1000),
    MCSE.lower = power - 1.96*MCSE,
    MCSE.upper = power + 1.96*MCSE
    ) |>
  filter(effect_size != -0.07)
  

# Plot results 
plotHA <- ggplot(subset(summary_HA, effect_size != 0.05), 
                 aes(x = sample_size, y = power, color = analysis)) +
  geom_point(position = position_dodge(width = 10)) +
  geom_line(position = position_dodge(width = 10)) +
  facet_wrap(~effect_size) + 
  geom_errorbar(aes(ymin = MCSE.lower, ymax = MCSE.upper), width = 10,
                position = position_dodge(width = 10)) +
  geom_hline(yintercept = 0.9, color = "red", lty = "dashed") + 
  theme_minimal() + scale_color_viridis_d(option = "H") +
  ylim(c(0,1)) +
  labs(x = "Total sample size", y = "Estimated power")
plotHA  
