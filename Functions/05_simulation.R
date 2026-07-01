# ── 1. Experiment grid ────────────────────────────────────────────────────────
set.seed(43)  # makes the grid itself reproducible

N_REPS <- 50  # replications per condition — adjust as needed

set.seed(42)

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
    rep     = NA  # no replication needed — graph structure is now exact
  )

# ── 2. Single-run function ────────────────────────────────────────────────────
# Wraps everything for one row of the grid.
# Returns a flat named list — becomes one row in your results table.

# In one_run(), iterate over ncomponent_grid rows like this:
one_run <- function(run_id, group_size, tau, rho, p_links,
                    nD, nA, nP, nS, nI, nR, rep, ...) {
  
  # E_candidates <- initial_graph(
  #   ncomponents = c(nD, nA, nP, nS, nI, nR),
  #   p_links     = p_links,
  #   signs       = c(.5,.5,.05,.5,.5,.05,.05,.5,.05),
  #   weighted    = FALSE,
  #   exact       = TRUE                               
  # )
  
  graph_out    <- initial_graph(
    ncomponents = c(nD, nA, nP, nS, nI, nR),
    p_links     = p_links,
    signs       = c(.5,.5,.05,.5,.5,.05,.05,.5,.05),
    weighted    = FALSE,
    exact       = TRUE
  )
  E_candidates <- graph_out$E_candidates
  V_nodes      <- graph_out$V
  feat_cols    <- graph_out$feat_cols
  
  # E_candidates <- initial_graph(...)
  cat("n_nodes at creation:", length(unique(E_candidates$u)), "\n")
  
  
  
  # --- Graph ------------------------------------------------------------------
  #  E_candidates <- initial_graph(
  #   ncomponents = c(3, 3, 3, 3, 3, 3),  # fixed for now
  #   p_links     = p_links,
  #   signs       = c(.5,.5,.05,.5,.5,.05,.05,.5,.05),
  #   weighted    = FALSE
  # )
  E_tbl <- moderate_E(E_candidates)
  
  # --- Agents -----------------------------------------------------------------
  agent_list <- generate_agents(group_size)
  agents <- lapply(seq_len(group_size), function(i)
    make_agent(agent_list$id[i], wvec(agent_list$w_mix[i]),
               E_candidates, alpha_default, feat_cols)
  )
  
  # --- Simulate ---------------------------------------------------------------
  rounds <- 500 * group_size  # scale rounds with group size
  res <- run_session(
    E_candidates, E_tbl, agents, alpha_default, feat_cols,
    tau = tau, rounds = rounds,
    lambda = 1.0, gamma = 1.0, rho = rho
  )
  
  # --- Summarise --------------------------------------------------------------
  
  shared_con <- detect_shared_consensus(res)
  shared_graph <- G_shared(res, E_candidates, V_nodes)
  n_nodes_shared <- sum(degree(shared_graph) > 0)
  
  shared_consensus_final <- shared_con$shared_consensus && n_nodes_shared > 0
  agent_con     <- detect_agent_consensus(res)
  wv_mix        <- characterize_worldview_mix(agent_list)
  evidence_met  <- edge_metrics(res, shared_graph, E_candidates, E_tbl)
  
  cat("n_nodes at return:", length(unique(E_candidates$u)), "\n")
  
  # --- Return a flat list (one row) -------------------------------------------
  list(
    run_id                 = run_id,
    group_size             = group_size,
    tau                    = tau,
    rho                    = rho,
    p_links                = p_links,
    n_nodes_original       = length(unique(E_candidates$u)),
    n_edges_original       = nrow(E_candidates),
    n_nodes_shared         = n_nodes_shared,
    n_edges_shared         = gsize(shared_graph),
    evidence_original      = mean(E_tbl$E),
    evidence_shared        = evidence_met$evidence_mean,
    DR_paths               = DAPSIR_path(shared_graph),
    shared_consensus      = shared_consensus_final,
    shared_consensus_time  = shared_con$shared_consensus_time,
    agent_consensus        = agent_con$agent_consensus,
    agent_consensus_time   = agent_con$agent_consensus_time,
    final_agreement        = tail(agent_shared_disagreement(res), 1),
    wv_diversity           = wv_mix$diversity,
    wv_tension             = wv_mix$tension,
    wv_dominant            = paste(wv_mix$dominant_worldviews, collapse = "&"),
    prop_H                 = wv_mix$proportions[1],
    prop_I                 = wv_mix$proportions[2],
    prop_F                 = wv_mix$proportions[3],
    prop_E                 = wv_mix$proportions[4],
    n_cycles               = length(FindCycles(shared_graph)),
    n_disconnected         = disconnected(shared_graph)
  )
}


test <- pmap(experiment_grid[51, ], one_run)[[1]]
str(test)

any(sapply(test$shared_history, function(x) any(is.na(x))))

pilot_test <- pmap(experiment_grid[94:97, ], one_run)
pilot_df <- bind_rows(pilot_test)
summary(pilot_df)


# A few checks before we start the parallel run
set.seed(1)
pilot_rows <- experiment_grid[sample(nrow(experiment_grid), 30), ]

system.time({
  pilot_results <- pmap(pilot_rows, one_run)
  pilot_df <- bind_rows(pilot_results)  
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
    .options = furrr_options(seed = TRUE))
})

# ── 4. Flatten and save ───────────────────────────────────────────────────────
results_df <- bind_rows(results_raw)

saveRDS(results_df, "simulation_results.rds")
write_csv(results_df, "simulation_results.csv")