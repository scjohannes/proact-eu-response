library(arrow)
library(here)
library(purrr)
library(dplyr)
library(markov.misc)
library(ggplot2)
library(fs)

# Pretend that we are in the proact-sim project
setwd("/home/felix/SHARE/imssrv_home/Simulations/EUPROACT Simulations/proact-sim")
#OUTPUT_PATH <- "/home/felix/SHARE/imssrv_home/Simulations/EUPROACT Simulations/proact-sim/scripts/output/sim051225"
OUTPUT_PATH <- "/home/felix/Simulations/EUPROACT Simulations/proact-sim/scripts/output"

# # Get all results that are present
# results_files <- dir_ls(OUTPUT_PATH, glob = "*simulation_*results_*.parquet")
# results_files <- results_files[!grepl("taooh|all", results_files)]
# all_results_combined <- map_dfr(results_files, read_parquet)
# 
# # Save combined results
# results_path_all <- path(OUTPUT_PATH, "simulation_results_all.parquet")
# write_parquet(all_results_combined, results_path_all)
all_results_combined0512 <- read_parquet(path(OUTPUT_PATH, 
                                              "simulation_results_all.parquet"))

# Get summary 
summary_HA0512 <- all_results_combined0512 |>
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
plotHA0512 <- ggplot(subset(summary_HA0512, effect_size != 0.05), 
                     aes(x = sample_size, y = power, color = analysis)) +
  geom_point(position = position_dodge(width = 10)) +
  geom_line(position = position_dodge(width = 10)) +
  facet_wrap(~effect_size) + 
  geom_errorbar(aes(ymin = MCSE.lower, ymax = MCSE.upper), width = 10,
                position = position_dodge(width = 10)) +
  geom_hline(yintercept = 0.025, color = "black") + 
  theme_minimal() + scale_color_viridis_d(option = "H") + ylim(c(0,1)) +
  labs(x = "Total sample size", y = "Power estimate")
plotHA0512

plotH00512 <- plot_results(data = subset(summary_HA0512, effect_size == 0.05), 
                           power, 
                           group = analysis,
                           x = sample_size,
                           combine = FALSE)

plotH00512$single_plots[[1]] + facet_wrap(~effect_size) + 
  geom_ribbon(aes(ymin = 0.025 - 1.96*sqrt(0.025*0.975/1000),
                  ymax = 0.025 + 1.96*sqrt(0.025*0.975/1000)),
              alpha = 0.01, color = NA) +
  geom_hline(yintercept = 0.025, color = "black") +
  ylim(c(0, 0.05)) + scale_color_viridis_d(option = "H") +
  labs(x = "Total sample size", y = "Type-I-error estimate")
