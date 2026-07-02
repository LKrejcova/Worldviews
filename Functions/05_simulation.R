
# ── 1. Single-run function ────────────────────────────────────────────────────
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
  # cat("n_nodes at creation:", length(unique(E_candidates$u)), "\n")
  
  
  
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
  
  original_graph <- E_candidates %>% select(u,v) %>% graph_from_data_frame(vertices = tibble(name=V_nodes), directed=TRUE)
  
  shared_consensus_final <- shared_con$shared_consensus && n_nodes_shared > 0
  agent_con     <- detect_agent_consensus(res)
  wv_mix        <- characterize_worldview_mix(agent_list)
  evidence_met  <- edge_metrics(res, shared_graph, E_candidates, E_tbl)
  
  n_cycles = tryCatch(
    R.utils::withTimeout(length(FindCycles(shared_graph)), timeout = 10, onTimeout = "error"),
    error = function(e) NA_integer_
  )
  
  # cat("n_nodes at return:", length(unique(E_candidates$u)), "\n")
  
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
    mean_evidence_original = mean(E_tbl$E),
    median_evidence_original = median(E_tbl$E),
    mean_evidence_shared   = evidence_met$evidence_mean,
    median_evidence_shared = evidence_met$evidence_median,
    DR_paths               = DAPSIR_path(shared_graph),
    shared_consensus       = shared_consensus_final,
    shared_consensus_time  = shared_con$shared_consensus_time,
    agent_consensus        = agent_con$agent_consensus,
    agent_consensus_time   = agent_con$agent_consensus_time,
    final_agreement        = agent_shared_disagreement(res),
    wv_diversity           = wv_mix$diversity,
    wv_tension             = wv_mix$tension,
    wv_dominant            = paste(wv_mix$dominant_worldviews, collapse = "&"),
    prop_H                 = wv_mix$proportions[1],
    prop_I                 = wv_mix$proportions[2],
    prop_F                 = wv_mix$proportions[3],
    prop_E                 = wv_mix$proportions[4],
    n_cycles               = n_cycles,
    n_disconnected         = disconnected(shared_graph),
    shared_graph           = shared_graph,
    original_graph         = original_graph
    
  )
}

