# This file contains the functions used to create agents

# Agent structure --------------------------------------------------------------
make_agent <- function(id, w_mix, E_candidates, alpha_list, feat_cols){
  # Initial map sampled from worldview-weighted priors (no evidence yet)
  pri <- vapply(seq_len(nrow(E_candidates)), function(i){
    row_vals <- unlist(E_candidates[i, feat_cols, drop = FALSE], use.names = FALSE)
    psi <- c(`(Intercept)`=1, setNames(as.numeric(row_vals), feat_cols))
    lp  <- edge_prior_logit(w_mix, psi, alpha_list)
    sigmoid(lp)
  }, numeric(1))
  
  include0 <- rbinom(length(pri), 1, pri)==1
  list(id=id, w_mix=w_mix, include=include0)
}

# Generate agents --------------------------------------------------------------
# function to generate agent IDs and their corresponding worldview
generate_agents <- function(n_agents, probs = c(HIER = 0.25, EGAL = 0.25, INDIV = 0.25, FATAL = 0.25)) {
  
  worldviews <- c("HIER", "EGAL", "INDIV", "FATAL")
  
  # Validate probs
  if (!all(names(probs) %in% worldviews) || length(probs) != 4) {
    stop("probs must be a named vector with names: HIER, EGAL, INDIV, FATAL")
  }
  if (abs(sum(probs) - 1) > 1e-6) {
    stop("probs must sum to 1")
  }
  
  # Sample worldviews
  assigned <- sample(worldviews, size = n_agents, replace = TRUE, prob = probs[worldviews])
  
  # Build IDs: track counter per worldview prefix
  prefix_map <- c(HIER = "H", EGAL = "E", INDIV = "I", FATAL = "F")
  counters <- c(H = 0L, E = 0L, I = 0L, F = 0L)
  
  ids <- vapply(assigned, function(wv) {
    prefix <- prefix_map[wv]
    counters[prefix] <<- counters[prefix] + 1L
    paste0(prefix, counters[prefix])
  }, character(1))
  
  data.frame(id = ids, w_mix = assigned, row.names = NULL)
}

# Proposal probability for an edge e by agent i at time t ----------------------
proposal_prob <- function(w_mix, psi_vec, alpha_list, E_e, social_sig, lambda=1, gamma=1){
  p_prior <- sigmoid(edge_prior_logit(w_mix, psi_vec, alpha_list))
  p_evid  <- sigmoid(lambda * E_e)
  p_soc   <- exp(gamma * social_sig)
  p_prior * p_evid * p_soc
}

# Create worldview mixtures (one-hot for simplicity) ---------------------------
wvec <- function(k){
  setNames(rep(0, length(worldviews)), worldviews) %>% {.[k] <- 1; .}
}

# Get worldview priors ---------------------------------------------------------
worldviews <- c("HIER","EGAL","INDIV","FATAL")

# Default coefficients for each worldview over psi(e) features (tune/fit with seed data - informed judgement)
alpha_default <- list(
  HIER = c(`(Intercept)`= -1.0, sign = 0.2, delay = -0.1, feat_authority= 1.0,
           feat_equity=-0.1, feat_market= 0.1, feat_risk=-0.2),
  EGAL = c(`(Intercept)`= -1.0, sign = 0.1, delay = 0.1,  feat_authority=-0.2,
           feat_equity= 1.1, feat_market=-0.3, feat_risk= 0.0),
  INDIV= c(`(Intercept)`= -1.0, sign = 0.2, delay = 0.0,  feat_authority=-0.1,
           feat_equity=-0.1, feat_market= 1.0, feat_risk=-0.1),
  FATAL= c(`(Intercept)`= -1.2, sign =-0.1, delay = 0.3,  feat_authority=-0.1,
           feat_equity= 0.0, feat_market=-0.1, feat_risk= 0.9)
)

# Compute prior log-odds for an edge given a worldview mix, the characteristics of that edge, and the coefficients for that worldview
edge_prior_logit <- function(w_mix, psi_vec, alpha_list){
  # w_mix: named vector over worldviews summing to 1
  # psi_vec: named vector of features with intercept
  sum(sapply(names(w_mix), function(k){
    ak <- alpha_list[[k]]
    # align names
    ak <- ak[names(psi_vec)]
    w_mix[k] * sum(ak * psi_vec)
  }))
}

# Social matrix W (row-stochastic) to inform influence -------------------------
# Here: homophily by worldview distance, we can set quite a few different models here depending on the social influence we want to consider
# options: homophily by worldviews, one dominant individual, faction formation, socnet with varying modularity

# homophily by worldviews
make_W <- function(wmix_mat, alpha=4){
  # wmix_mat: N x K matrix of worldview mixtures
  N <- nrow(wmix_mat)
  D <- as.matrix(dist(wmix_mat, method="manhattan"))
  W <- exp(-alpha*D)
  diag(W) <- 0
  W <- W/rowSums(W + 1e-12)
  W
}

## competence
# assign competence to each agent


# competence AND homophily
make_W_comp <- function(wmix_mat, alpha = 4){
  N <- nrow(wmix_mat)
  # randomly assign competence
  competence <- rbeta(N, shape1 = 2, shape2 = 3)
  
  # define competence preferences
  comp_pref <- c(
    "HIER" = 1,
    "EGAL" = 0.1,
    "INDIV" = 0.5,
    "FATAL" = 0.1
  )
  
  # determine what the observer's competence preference is (if worldview mix)
  obs_pref <- as.numeric(wmix_mat %*% comp_pref)
  
  # define homophily
  D <- as.matrix(dist(wmix_mat, method="manhattan"))
  W <- exp(-alpha*D)
  diag(W) <- 0
  
  # temper homophily by competence
  for(i in 1:N){
    W[i, ] <- W[i, ] * (1 + obs_pref[i] * competence) #baseline influence is determined by worldview rather than competence
  }
  # normalize
  W <- W/rowSums(W + 1e-12)
  W
}

