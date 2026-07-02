# Social update dynamics -------------------------------------------------------
aggregate_group <- function(indiv_mat, weights, tau){
  # indiv_mat: edges x N logical matrix
  p <- as.numeric(indiv_mat %*% normalize(weights))
  as.logical(p >= tau)
}

degroot_update <- function(indiv_mat, W, rho=1){ #rho changes
  # indiv_mat: edges x N (0/1); W: N x N row-stochastic
  M <- apply(indiv_mat, 2, as.numeric)       # edges x N retrieves the individ matrix
  M_new <- M %*% t(W)                        # edges x N multiplies by the influence weights (worldview x worldview)
  M_out <- (1-rho)*M + rho*M_new # add the social influences to the indiv matrices, weighted by 1-rho and rho, respectively
  (M_out >= 0.5)  # only keep the edges that are more than 50% "on"
}

# Run a session ----------------------------------------------------------------
# run_session <- function(E_candidates, E_score, agents, alpha_list, feat_cols,
#                         tau=0.6, rounds=200, lambda=1, gamma=1, rho=0.51, W=NULL, save_snapshots = TRUE){
#   # set up row and column number for future matrices
#   N <- length(agents)
#   edges <- nrow(E_candidates)
#   
#   # make the worldview matrix
#   wmix_mat <- do.call(rbind, lapply(agents, function(a) as.numeric(a$w_mix)))
#   rownames(wmix_mat) <- sapply(agents, `[[`, "id")
#   # W <- make_W(wmix_mat)
#   W <- make_W_comp(wmix_mat, alpha = 4)
#   
#   # make the matrix of individual graphs: edges x agents 
#   indiv <- do.call(cbind, lapply(agents, function(a) as.numeric(a$include)))
#   indiv <- matrix(indiv, nrow=edges, ncol=N)
#   weights <- rep(1, N)
#   
#   # make the shared graph based on the individual graphs
#   shared <- aggregate_group(indiv, weights, tau)
#   
#   # retrieve edge features: features, evidence scores, and proportion of agents who include the edge
#   psi_list <- lapply(1:edges, function(i) psi_of(E_candidates[i,], feat_cols))
#   E_vec <- E_score$E[match(1:edges, E_score$edge_id)]
#   social <- rowMeans(indiv)
#   
#   # Only allocate snapshot storage if needed
#   round_snapshots <- if (save_snapshots) vector("list", rounds) else NULL
#   
#   # Lightweight tracking that works regardless of save_snapshots
#   shared_history     <- vector("list", rounds)   # just the shared logical vector each round
#   jaccard_trace    <- numeric(rounds) # mean Jaccard between each agent and the shared graph
#   
#   
#   for (t in 1:rounds){
#     speaker <- ((t - 1) %% N) + 1 # select one speaker
#     probs <- sapply(1:edges, function(ei){ # calculate the probability of selecting a given edge => which edge do I pay most attention to?
#       proposal_prob(agents[[speaker]]$w_mix, psi_list[[ei]], alpha_list,
#                     E_vec[ei], social[ei], lambda, gamma)
#     })
#     probs <- probs/sum(probs) # normalize the probability
#     e_star <- sample(1:edges, 1, prob=probs) # sample one edge following the probability distribution
#     indiv[e_star, speaker] <- TRUE
#     # if the speaker included e_star originally: indiv[e_star, speaker] == 1 => as.logical => TRUE, !as.logical => FALSE
#     # if the speaker didn't include e_star: indiv[e_star, speaker] == 0 => as.logical => FALSE, !as.logical => TRUE
#     # basically the speaker toggles the edge on/off
#     indiv <- degroot_update(indiv, W, rho) # agents all update their graph
#     shared <- aggregate_group(indiv, weights, tau) # update shared graph
#     social <- rowMeans(indiv) # recalculate the proportion of agents who include the edge
#     
#     # Always track these lightweight summaries
#     shared_history[[t]] <- shared
#     jaccard_trace[t] <- mean(apply(indiv, 2, function(a) jaccard(a, shared)))
#     
#     # Only store full snapshots if requested
#     if (save_snapshots){
#       round_snapshots[[t]] <- list(indiv=indiv, shared=shared)
#     }
#     
#   }
#   
#   
#   list(
#     shared               = shared,
#     indiv                = indiv,
#     W                    = W,
#     snapshots          = round_snapshots,       # NULL if save_snapshots=FALSE
#     shared_history     = shared_history,        # always available, cheap to store
#     jaccard_trace = jaccard_trace      # always
#   )
# }

run_session <- function(E_candidates, E_score, agents, alpha_list, feat_cols,
                        tau=0.6, rounds=200, lambda=1, gamma=1, rho=0.51, W=NULL, 
                        save_snapshots = FALSE){
  
  N     <- length(agents)
  edges <- nrow(E_candidates)
  
  # Worldview matrix and influence weights
  wmix_mat <- do.call(rbind, lapply(agents, function(a) as.numeric(a$w_mix))) # get the worldviews of each agent
  rownames(wmix_mat) <- sapply(agents, `[[`, "id")
  W <- make_W_comp(wmix_mat, alpha = 4) # get the worldview x competence matrix
  
  # Individual graph matrix: edges x agents
  indiv   <- do.call(cbind, lapply(agents, function(a) as.numeric(a$include))) # get the initial graph the agents build
  indiv   <- matrix(indiv, nrow = edges, ncol = N) # save it as a matrix
  weights <- rep(1, N) # change if you want the agents to  have different weights
  
  # Shared graph and social signal
  shared <- aggregate_group(indiv, weights, tau) # make the shared graph
  social <- rowMeans(indiv) # get the proportion of agents who include each edge
  
  # Evidence vector
  E_vec <- E_score$E[match(1:edges, E_score$edge_id)] # get the evidence score for each edge
  
  # ── Precompute once before the loop ────────────────────────────────────────
  
  # psi_matrix: edges x (1 + n_features)
  # replaces psi_list + the per-edge unlist/setNames inside the loop
  psi_matrix <- cbind(
    `(Intercept)` = 1,
    as.matrix(E_candidates[, feat_cols])
  ) # extracts the features of each edge and puts them into a matrix (where is the sign though?)
  
  # Weighted alpha vector per agent: collapses the sapply over worldviews
  # into a single named numeric vector per agent
  alpha_per_agent <- lapply(agents, function(a) {
    Reduce("+", lapply(names(a$w_mix), function(k) {
      a$w_mix[k] * alpha_list[[k]]
    }))
  })  # extracts 
  
  # Evidence probabilities don't change between rounds
  p_evid <- sigmoid(lambda * E_vec)
  
  # ── Storage ────────────────────────────────────────────────────────────────
  round_snapshots <- if (save_snapshots) vector("list", rounds) else NULL
  shared_history  <- vector("list", rounds)
  jaccard_trace   <- numeric(rounds)
  
  # ── Main loop ──────────────────────────────────────────────────────────────
  for (t in 1:rounds) {
    speaker <- ((t - 1) %% N) + 1
    
    # Single matrix multiply replaces sapply(1:edges, ...) + edge_prior_logit
    log_odds <- as.numeric(psi_matrix %*% alpha_per_agent[[speaker]]) # multiply the agent coefficients with the edge features
    p_prior  <- sigmoid(log_odds)
    p_soc    <- exp(gamma * social)   # social updates every round
    
    probs  <- p_prior * p_evid * p_soc
    probs  <- probs / sum(probs)
    e_star <- sample(1:edges, 1, prob = probs)
    
    indiv[e_star, speaker] <- TRUE
    indiv  <- degroot_update(indiv, W, rho)
    shared <- aggregate_group(indiv, weights, tau)
    social <- rowMeans(indiv)
    
    # Lightweight tracking
    shared_history[[t]] <- shared
    jaccard_trace[t]    <- mean(apply(indiv, 2, function(a) jaccard(a, shared)))
    
    if (save_snapshots) {
      round_snapshots[[t]] <- list(indiv = indiv, shared = shared)
    }
  }
  
  list(
    shared         = shared,
    indiv          = indiv,
    W              = W,
    snapshots      = round_snapshots,
    shared_history = shared_history,
    jaccard_trace  = jaccard_trace
  )
}
