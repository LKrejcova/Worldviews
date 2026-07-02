library(tidyverse)
library(igraph)
library(glue)
library(ggraph)
library(scales)
library(reshape2)
library(mlf)
library(infotheo)
library(furrr)
library(future)
library(ggridges)
library(extraDistr)
library(R.utils)

source("./Functions/00_helpers.R")
source("./Functions/01_graph_generation.R")
source("./Functions/02_agent_dynamics.R")
source("./Functions/03_session_dynamics.R")
source("./Functions/04_metrics.R")
source("./Functions/05_simulation.R")


# ── 1. Experiment grid ────────────────────────────────────────────────────────
# set.seed(43)  # makes the grid itself reproducible

N_REPS <- 50  # replications per condition — adjust as needed

experiment_grid <- tidyr::crossing(
  nD         = 2:5,
  nA         = 2:5,
  nP         = 2:5,
  nS         = 2:5,
  nI         = 2:5,
  nR         = 2:5,
  group_size = c(3, 5, 10, 15, 20, 25),
  tau        = c(0.5, 0.66)
) |>
  mutate(
    run_id  = row_number(),
    rho     = runif(n(), 0.5, 1),
    p_links = runif(n(), 0.25, 0.75),
    replicate     = 1:3
  )

experiment_grid <- tidyr::crossing(
  nD         = 2:5,
  nA         = 2:5,
  nP         = 2:5,
  nS         = 2:5,
  nI         = 2:5,
  nR         = 2:5,
  group_size = c(3, 5, 10, 15, 20, 25),
  tau        = c(0.5, 0.66)
) |>
  tidyr::crossing(replicate = 1:3) |>
  mutate(
    run_id  = row_number(),
    rho     = runif(n(), 0.5, 1),
    p_links = runif(n(), 0.25, 0.75)
  )

# 3 replicates per parameter setup
experiment_grid <- tidyr::crossing(
  nD         = 2:5,
  nA         = 2:5,
  nP         = 2:5,
  nS         = 2:5,
  nI         = 2:5,
  nR         = 2:5,
  group_size = c(3, 5, 10, 15, 20, 25),
  tau        = c(0.5, 0.66)
) |>
  mutate(
    rho     = runif(n(), 0.5, 1),
    p_links = runif(n(), 0.25, 0.75)
  ) |>
  tidyr::crossing(replicate = 1:3) |>
  mutate(run_id = row_number())


test <- pmap(experiment_grid[156, ], one_run)[[1]]
str(test)

any(sapply(test$shared_history, function(x) any(is.na(x))))

pilot_test <- pmap(experiment_grid[1:3, ], one_run)
graph_cols <- c("shared_graph", "original_graph")

# Data frame of scalar values
pilot_df <- map_dfr(pilot_test, ~{
  as_tibble(.x[setdiff(names(.x), graph_cols)])
})

# Separate list of graphs
pilot_graphs <- map(pilot_test, ~{
  .x[graph_cols]
})

pilot_graphs <- setNames(
  map(pilot_test, ~ .x[graph_cols]),
  map_chr(pilot_test, ~ as.character(.x$run_id))
)


# A few checks before we start the parallel run
set.seed(1)
pilot_rows <- experiment_grid[sample(nrow(experiment_grid), 100), ]

system.time({
  pilot_results <- pmap(pilot_rows, one_run) # apply one_run to each row of pilot_rows
  pilot_df <- bind_rows(pilot_results) 
  graphs <- pilot_results[]
})


# Sanity checks
summary(pilot_df)
pilot_df |> count(n_nodes_original)          # should vary, not be stuck at one value
pilot_df |> filter(is.na(shared_consensus_time)) |> nrow()  # how many never converge?
range(pilot_df$n_nodes_shared)

# ── 3. Run in parallel ────────────────────────────────────────────────────────
plan(multisession, workers = parallel::detectCores() - 1)

system.time({
  results_raw <- future_pmap(
    pilot_rows,
    one_run,
    .options = furrr_options(seed = FALSE),
    .progress = TRUE)
})

# ── 4. Flatten and save ───────────────────────────────────────────────────────
results_df <- bind_rows(results_raw)

# Data frame of scalar values
graph_cols <- c("shared_graph", "original_graph")

results_df <- map_dfr(results_raw, ~{
  as_tibble(.x[setdiff(names(.x), graph_cols)])
})

# Separate list of graphs
results_graphs <- map(results_raw, ~{
  .x[graph_cols]
})

saveRDS(results_df, "./Results/simulation_results.rds")
write_csv(results_df, "./Results/simulation_results.csv")

saveRDS(results_graphs, "./Results/Simulation_graphs.rds")
