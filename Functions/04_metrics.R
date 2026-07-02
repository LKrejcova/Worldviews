G_shared <- function(res, E_candidates, V){
  inc <- which(res$shared) # gets edge_id for the edges included in the shared graph
  E_candidates %>% filter(edge_id %in% inc) %>% select(u,v) %>% graph_from_data_frame(vertices = tibble(name=V), directed=TRUE)
}

DAPSIR_path <- function(graph){
  # Vertex names
  D <- grep("^D", V(graph)$name, value = TRUE)
  R <- grep("^R", V(graph)$name, value = TRUE)
  
  # Initialize matrix
  reachable <- matrix(
    0,
    nrow = length(D),
    ncol = length(R),
    dimnames = list(D, R)
  )
  
  # Fill matrix
  for (d in D) {
    for (r in R) {
      paths <- all_shortest_paths(graph, from = d, to = r)$res
      reachable[d, r] <- as.integer(length(paths) > 0)
    }
  }
  
  dimnames(reachable) <- list(D,R)
  reachability <- sum(reachable) 
  return(reachability)
}

# let's look at the number of nodes that have a degree > 0 (i.e. they are connected to the graph)
disconnected <- function(graph){
  sum(degree(graph)==0)
}

edge_metrics <- function(res, graph, E_candidates, E_score){
  df <- E_candidates %>% mutate(include = as.logical(res$shared)) %>% left_join(E_score, by="edge_id")
  tibble(
    n_edges = sum(df$include),
    n_nodes = sum(degree(graph) > 0),
    evidence_mean = mean(df$E[df$include]),
    evidence_median = median(df$E[df$include])
  )
}

############################ Edge density ###################################

# edge density between component types
get_edge_density <- function(results){
  shared_edge_features <- E_candidates %>%
    mutate(incl_shared = results$shared,# attach the included edges to the features dataframe
           from = substr(u,1,1), # remove the numbers from the components
           to = substr(v,1,1)
    )
  
  # need to do this so that unrealized edges are still visible
  allowed <- tibble(
    from = c("D","A","P","S","I","R","R","R","R"),
    to   = c("A","P","S","I","R","D","A","S","I")
  )
  
  # get the metrics
  result <- shared_edge_features %>%
    group_by(from, to) %>%
    summarise(
      density = mean(incl_shared),
      n_possible = n(),
      n_present = sum(incl_shared),
      .groups = "drop"
    ) %>%
    right_join(allowed, by = c("from", "to"))
}

# Now we look at the features
get_feature_metrics <- function(results){
  shared_edge_features <- E_candidates %>%
    mutate(incl_shared = results$shared)
  
  result <- shared_edge_features %>%
    group_by(u_feature, v_feature) %>%
    summarise(
      density = mean(incl_shared),
      n_possible = n(),
      n_present = sum(incl_shared),
      .groups = "drop"
    )
}

# Consensus detection ----------------------------------------------------------
detect_shared_consensus <- function(res) {
  shared_change <- sapply(2:length(res$shared_history), function(t) {
    sum(res$shared_history[[t]] != res$shared_history[[t - 1]])
  })
  is_zero <- shared_change == 0
  n <- length(is_zero)
  
  # Find the last position where a change occurred; consensus starts right after
  last_change <- max(which(!is_zero), 0)  # 0 if no changes at all
  consensus_time <- if (last_change < n) last_change + 1 else NA_integer_
  
  shared_consensus <- !is.na(consensus_time) &&
    consensus_time < 0.97 * length(res$shared_history)
  
  list(
    shared_consensus      = shared_consensus,
    shared_consensus_time = consensus_time
  )
}

# Jaccard index allows for better comparison between graphs of different size than the Hamming distance
agent_shared_disagreement <- function(results) {
  
  indiv  <- results$indiv
  shared <- results$shared
  
  mean(
    apply(indiv, 2, function(agent_graph) {
      jaccard(agent_graph, shared)
    })
  )
}


detect_agent_consensus <- function(res) {
  consensus_time <- which(
    sapply(seq_along(res$jaccard_trace), function(t) {
      all(res$jaccard_trace[t:length(res$jaccard_trace)] == 1)
    })
  )[1]
  list(
    agent_consensus      = !is.na(consensus_time),
    agent_consensus_time = consensus_time
  )
}

# Characterise the worldview mix of the agents ---------------------------------
characterize_worldview_mix <- function(agent_list) {
  
  # --- Proportions ---
  counts <- table(factor(agent_list$w_mix, levels = c("HIER", "INDIV", "FATAL", "EGAL")))
  p <- counts / sum(counts)
  
  pH <- p["HIER"]  # Hierarchist:   HiGrid-HiGroup
  pI <- p["INDIV"]  # Individualist: LoGrid-LoGroup
  pF <- p["FATAL"]  # Fatalist:      HiGrid-LoGroup
  pE <- p["EGAL"]  # Egalitarian:   LoGrid-HiGroup
  
  # --- Dominance ---
  p_max     <- max(p)
  top_label <- names(p)[which.max(p)]
  DR        <- p_max / 0.25 # DR = 1 => no dominance, DR = 4 => maximum dominant
  
  # Co-dominance: top two together > 70% and within 10pp of each other
  p_sorted    <- sort((p), decreasing = TRUE)
  co_dominant <- (p_sorted[1] + p_sorted[2] > 0.70) &
    (p_sorted[1] - p_sorted[2] < 0.10)
  
  top2_labels <- names(sort(p, decreasing = TRUE))[1:2]
  
  
  if (isTRUE(co_dominant)) {
    dominant_worldviews <- top2_labels
  } else {
    dominant_worldviews <- top_label
  }
  
  # no dominance: when all worldviews present have the same proportion (0.25, 0.33, or 0.5)
  
  no_dominance <- p_sorted[1] == p_sorted[2]
  
  if (isTRUE(no_dominance)) {
    dominant_worldviews <- "NoDom"
  }
  
  if(p_sorted[1] == 1) {
    dominant_worldviews <- "Unanim"
  }
  
  # --- Diversity (Simpson's Dominance) ---
  # Ok now we replace this with Simpson
  p_nonzero <- p[p > 0] # get non-zero proportions
  
  # Calculate D
  D <- sum(p_nonzero^2)
  
  
  # --- Tension (diagonal opposition) ---
  T_raw <- pH * pI + pF * pE
  T     <- T_raw / 0.5  # normalize to [0, 1]
  
  # --- Output ---
  list(
    proportions         = as.numeric(p),
    names               = names(p),
    dominant_worldviews = dominant_worldviews,
    dominance_ratio     = round(as.numeric(DR), 3),
    dominance           = round(as.numeric(D), 3),
    diversity           = round(as.numeric(1-D),3),
    tension             = round(as.numeric(T), 3)
  )
}