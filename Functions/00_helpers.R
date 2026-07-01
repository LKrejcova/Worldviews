
sigmoid <- function(x) 1/(1+exp(-x))
normalize <- function(x) x/sum(x)

# Jaccard similarity between two logical vectors
jaccard <- function(a, b) {
  a <- as.logical(a); b <- as.logical(b)
  intersection <- sum(a & b)
  union   <- sum(a | b)
  if (union == 0) return(NA_real_)
  intersection/union
}

# Sorensen-Dice index
sorensen_dice <- function(a,b){
  J <- jaccard(a,b)
  2 * J / (1 + J)
}

# Simple cycle counter up to length L (approx; for diagnostics)
count_cycles <- function(g, max_len=4){
  cycles <- list()
  for (k in 2:max_len){
    cyc_k <- suppressWarnings(igraph::cycles(g, maxlen=k))
    if (!is.null(cyc_k)) cycles[[as.character(k)]] <- length(cyc_k)
  }
  tibble(k = names(cycles), n = unlist(cycles))
}

# finds all cycles (could be an issue if the graph is too big)
FindCycles = function(g) {
  Cycles = NULL
  for(v1 in V(g)) {
    if(degree(g, v1, mode="in") == 0) { next }
    GoodNeighbors = neighbors(g, v1, mode="out")
    GoodNeighbors = GoodNeighbors[GoodNeighbors > v1]
    for(v2 in GoodNeighbors) {
      TempCyc = lapply(all_simple_paths(g, v2,v1, mode="out"), function(p) c(v1,p))
      TempCyc = TempCyc[which(sapply(TempCyc, length) > 3)]
      TempCyc = TempCyc[sapply(TempCyc, min) == sapply(TempCyc, `[`, 1)]
      Cycles  = c(Cycles, TempCyc)
    }
  }
  Cycles
}