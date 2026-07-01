# These functions serve to generate the initial DAPSIR graph the agents look at

## Make DAPSIR graph -----------------------------------------------------------
makeDAPSIR <- function(ncomponents=c(3,6,4,4,4,3), p_links=0.25, 
                       signs=c(.75,.95,.05,.95,.95,.05,.05,.95,.05),
                       weighted=FALSE, exact=FALSE) {  # add exact argument
  require(extraDistr)
  
  if (exact) {
    # Use ncomponents as exact node counts
    nD <- ncomponents[1]
    nA <- ncomponents[2]
    nP <- ncomponents[3]
    nS <- ncomponents[4]
    nI <- ncomponents[5]
    nR <- ncomponents[6]
  } else {
    # Original behaviour: treat ncomponents as Poisson means
    nD <- rtpois(1, ncomponents[1], a=0)
    nA <- rtpois(1, ncomponents[2], a=0)
    nP <- rtpois(1, ncomponents[3], a=0)
    nS <- rtpois(1, ncomponents[4], a=0)
    nI <- rtpois(1, ncomponents[5], a=0)
    nR <- rtpois(1, ncomponents[6], a=0)
  }
  
  nstates<-nD+nA+nP+nS+nI+nR # total number of components
  DPSIR<-matrix(0,nstates,nstates) # create null matrix with proper dimensions
  rownames(DPSIR)<-c(rep("D",nD),rep("A",nA),rep("P",nP),rep("S",nS),rep("I",nI),rep("R",nR))
  colnames(DPSIR)<-c(rep("D",nD),rep("A",nA),rep("P",nP),rep("S",nS),rep("I",nI),rep("R",nR))
  
  #row feeds to columns
  
  #network topology
  # here we control the network density
  #p_links #a quarter of possible links are linked 
  #because then the p can be applied independently to all links, makes complex matrix subsetting easier
  # note we can change to make p component specific if we want
  #possible non-zero matrix elements
  ##### addition we need to ensure that all drivers are connected to at least one pressure
  ### all pressures are connected to at least one state
  ### all states are connected to at least one impact
  # rowsums >0
  # we have to sample elements instead to ensure at least one link per row AND one link per column
  
  #########################
  ## careful here when p_links becomes small that the while loop do not go to infinite!!
  ## update 7 Mar, got ride of naughty while loops, but dealing with (a more appropriate)
  ## restricted manipulation of p_links. that is because we need to have a minimal connection set
  ## at least one driver for all acitivities and all drivers connected to some activity, etc
  ########################
  
  # while (all(rowSums(as.matrix(DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A")]))!=0)==FALSE) {
  # DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A")]<-rbinom(prod(dim(as.matrix(DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A")]))),1,p_links)
  # }
  for (i in 1:nD) {
    DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A")][i,]<-replace(DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A")][i,],c(sample(seq(1:nA),ceiling(nA*p_links),replace=F)),1)
  }
  if (length(which(colSums(DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A"),drop=FALSE])==0))>0) {
    grab<-length(which(colSums(DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A"),drop=FALSE])==0))
    cols<-which(colSums(DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A"),drop=FALSE])==0)
    DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A")]<-replace(DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A")],cbind(sample(seq(1:nD),grab,replace=T),cols),1)
  }
  
  for (i in 1:nA) {
    DPSIR[which(rownames(DPSIR)=="A"),which(colnames(DPSIR)=="P")][i,]<-replace(DPSIR[which(rownames(DPSIR)=="A"),which(colnames(DPSIR)=="P")][i,],c(sample(seq(1:nP),ceiling(nP*p_links),replace=F)),1)
  }
  
  if (length(which(colSums(DPSIR[which(rownames(DPSIR)=="A"),which(colnames(DPSIR)=="P"),drop=FALSE])==0))>0) { # detects empty columns
    grab<-length(which(colSums(DPSIR[which(rownames(DPSIR)=="A"),which(colnames(DPSIR)=="P"),drop=FALSE])==0)) # gets the number of empty columns
    cols<-which(colSums(DPSIR[which(rownames(DPSIR)=="A"),which(colnames(DPSIR)=="P"),drop=FALSE])==0) # gets the indices of the empty columns
    DPSIR[which(rownames(DPSIR)=="A"),which(colnames(DPSIR)=="P")]<-replace(DPSIR[which(rownames(DPSIR)=="A"),which(colnames(DPSIR)=="P")],cbind(sample(seq(1:nA),grab,replace=T),cols),1)	}
  
  for (i in 1:nP) {
    DPSIR[which(rownames(DPSIR)=="P"),which(colnames(DPSIR)=="S")][i,]<-replace(DPSIR[which(rownames(DPSIR)=="P"),which(colnames(DPSIR)=="S")][i,],c(sample(seq(1:nS),ceiling(nS*p_links),replace=F)),1)
  }
  if (length(which(colSums(DPSIR[which(rownames(DPSIR)=="P"),which(colnames(DPSIR)=="S"),drop=FALSE])==0))>0) {
    grab<-length(which(colSums(DPSIR[which(rownames(DPSIR)=="P"),which(colnames(DPSIR)=="S"),drop=FALSE])==0))
    cols<-which(colSums(DPSIR[which(rownames(DPSIR)=="P"),which(colnames(DPSIR)=="S"),drop=FALSE])==0)
    DPSIR[which(rownames(DPSIR)=="P"),which(colnames(DPSIR)=="S")]<-replace(DPSIR[which(rownames(DPSIR)=="P"),which(colnames(DPSIR)=="S")],cbind(sample(seq(1:nP),grab,replace=T),cols),1)
  }
  
  for (i in 1:nS) {
    DPSIR[which(rownames(DPSIR)=="S"),which(colnames(DPSIR)=="I")][i,]<-replace(DPSIR[which(rownames(DPSIR)=="S"),which(colnames(DPSIR)=="I")][i,],c(sample(seq(1:nI),ceiling(nI*p_links),replace=F)),1)
  }
  if (length(which(colSums(DPSIR[which(rownames(DPSIR)=="S"),which(colnames(DPSIR)=="I"),drop=FALSE])==0))>0) {
    grab<-length(which(colSums(DPSIR[which(rownames(DPSIR)=="S"),which(colnames(DPSIR)=="I"),drop=FALSE])==0))
    cols<-which(colSums(DPSIR[which(rownames(DPSIR)=="S"),which(colnames(DPSIR)=="I"),drop=FALSE])==0)
    DPSIR[which(rownames(DPSIR)=="S"),which(colnames(DPSIR)=="I")]<-replace(DPSIR[which(rownames(DPSIR)=="S"),which(colnames(DPSIR)=="I")],cbind(sample(seq(1:nS),grab,replace=T),cols),1)}
  
  for (i in 1:nI) {
    DPSIR[which(rownames(DPSIR)=="I"),which(colnames(DPSIR)=="R")][i,]<-replace(DPSIR[which(rownames(DPSIR)=="I"),which(colnames(DPSIR)=="R")][i,],c(sample(seq(1:nR),ceiling(nR*p_links),replace=F)),1)
  }
  if (length(which(colSums(DPSIR[which(rownames(DPSIR)=="I"),which(colnames(DPSIR)=="R"),drop=FALSE])==0))>0) {
    grab<-length(which(colSums(DPSIR[which(rownames(DPSIR)=="I"),which(colnames(DPSIR)=="R"),drop=FALSE])==0))
    cols<-which(colSums(DPSIR[which(rownames(DPSIR)=="I"),which(colnames(DPSIR)=="R"),drop=FALSE])==0)
    DPSIR[which(rownames(DPSIR)=="I"),which(colnames(DPSIR)=="R")]<-replace(DPSIR[which(rownames(DPSIR)=="I"),which(colnames(DPSIR)=="R")],cbind(sample(seq(1:nI),grab,replace=T),cols),1)	}
  
  for (i in 1:nR) {
    DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)!="R")][i,]<-replace(DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)!="R")][i,],c(sample(seq(1:(nstates-nR)),ceiling((nstates-nR)*p_links),replace=F)),1)
  }
  
  # with that last one we don't constrain where we apply responses, we can change that of course
  
  # #### make sure all components are attached to the previous components
  # unthetered<-which(colSums(DPSIR[1:(nstates-nR),])[(nD+1):nstates]==0)
  # levels<-c("D","A","P","S","I","R")
  # if (length(unthetered)>0) {
  # for (i in 1:length(unthetered)) {
  # from.l<-levels[which(levels==names(unthetered)[i])-1]
  # col.l<-nD+as.numeric(unthetered[i])
  # DPSIR[sample(which(rownames(DPSIR)==from.l),replace=T,1),col.l]<-1
  # }
  # }
  
  
  #network signs
  #here we indirectly control network stability
  #p positive
  Sda<-signs[1] #not quite sure of this one
  Sap<-signs[2]
  Sps<-signs[3]
  Ssi<-signs[4]
  Sir<-signs[5]
  Srd<-signs[6]
  Srp<-signs[7]
  Srs<-signs[8]
  Sri<-signs[9]
  
  nz<-which(DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A")]==1)
  DPSIR[which(rownames(DPSIR)=="D"),which(colnames(DPSIR)=="A")][nz]<-sample(c(-1,1),length(nz),replace=TRUE,prob=c((1-Sda),Sda))
  
  nz<-which(DPSIR[which(rownames(DPSIR)=="A"),which(colnames(DPSIR)=="P")]==1)
  DPSIR[which(rownames(DPSIR)=="A"),which(colnames(DPSIR)=="P")][nz]<-sample(c(-1,1),length(nz),replace=TRUE,prob=c((1-Sap),Sap))
  
  nz<-which(DPSIR[which(rownames(DPSIR)=="P"),which(colnames(DPSIR)=="S")]==1)
  DPSIR[which(rownames(DPSIR)=="P"),which(colnames(DPSIR)=="S")][nz]<-sample(c(-1,1),length(nz),replace=TRUE,prob=c((1-Sps),Sps))
  
  nz<-which(DPSIR[which(rownames(DPSIR)=="S"),which(colnames(DPSIR)=="I")]==1)
  DPSIR[which(rownames(DPSIR)=="S"),which(colnames(DPSIR)=="I")][nz]<-sample(c(-1,1),length(nz),replace=TRUE,prob=c((1-Ssi),Ssi))
  
  nz<-which(DPSIR[which(rownames(DPSIR)=="I"),which(colnames(DPSIR)=="R")]==1)
  DPSIR[which(rownames(DPSIR)=="I"),which(colnames(DPSIR)=="R")][nz]<-sample(c(-1,1),length(nz),replace=TRUE,prob=c((1-Sir),Sir))
  
  nz<-which(DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)=="D")]==1)
  DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)=="D")][nz]<-sample(c(-1,1),length(nz),replace=TRUE,prob=c((1-Srd),Srd))
  
  nz<-which(DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)=="P")]==1)
  DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)=="P")][nz]<-sample(c(-1,1),length(nz),replace=TRUE,prob=c((1-Srp),Srp))
  
  nz<-which(DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)=="S")]==1)
  DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)=="S")][nz]<-sample(c(-1,1),length(nz),replace=TRUE,prob=c((1-Srs),Srs))
  
  nz<-which(DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)=="I")]==1)
  DPSIR[which(rownames(DPSIR)=="R"),which(colnames(DPSIR)=="I")][nz]<-sample(c(-1,1),length(nz),replace=TRUE,prob=c((1-Sri),Sri))
  
  if (weighted) {
    #network weights
    # here we "complexify" the network representation
    #weights can be functions to simplify the Jacobian do we stick to polynomials? 
    # P~aD+b | P~aD^2+bD+c
    # we can still include a limiting behaviour with polynomials:
    #dP/dt~aP^2(1-bD-P)+c or something like that
    
    #### phase I - constant (linear relationships)
    
    DPSIR[DPSIR!=0]<-DPSIR[DPSIR!=0]*runif(length(DPSIR[DPSIR!=0]),1e-12,1)
  }
  
  rownames(DPSIR)<-c(paste0(rep("D",nD),1:nD),paste0(rep("A",nA),1:nA),paste0(rep("P",nP),1:nP),paste0(rep("S",nS),1:nS),paste0(rep("I",nI),1:nI),paste0(rep("R",nR),1:nR))
  colnames(DPSIR)<-c(paste0(rep("D",nD),1:nD),paste0(rep("A",nA),1:nA),paste0(rep("P",nP),1:nP),paste0(rep("S",nS),1:nS),paste0(rep("I",nI),1:nI),paste0(rep("R",nR),1:nR))
  
  
  return(DPSIR)
}

# Draw the initial graph, assign features to each node -------------------------

initial_graph <- function(ncomponents, p_links, signs, weighted, exact){
  startDAPSIR <- makeDAPSIR(
    ncomponents = ncomponents,
    p_links     = p_links,
    signs       = signs,
    weighted    = FALSE,
    exact       = TRUE
  )
  
  # Convert matrix to edge list
  idx <- which(startDAPSIR != 0, arr.ind = TRUE)
  
  E_candidates <- data.frame(
    u    = rownames(startDAPSIR)[idx[,1]],
    v    = colnames(startDAPSIR)[idx[,2]],
    sign = startDAPSIR[idx]
  )
  
  # Assign random features to components
  components <- unique(c(E_candidates$u, E_candidates$v))
  features   <- c("Authority", "Equity", "Market", "Risk")
  component_features <- data.frame(
    component = components,
    feature   = sample(features, length(components), replace = TRUE)
  )
  
  # Add u and v features
  E_candidates <- E_candidates %>%
    left_join(component_features, by = c("u" = "component")) %>%
    rename(u_feature = feature) %>%
    left_join(component_features, by = c("v" = "component")) %>%
    rename(v_feature = feature)
  
  feat_cols <- c("sign","delay","feat_authority","feat_equity","feat_market","feat_risk")
  
  # Add edge-level feature flags
  E_candidates <- E_candidates %>%
    mutate(
      feat_authority = if_else(u_feature == "Authority" | v_feature == "Authority", 1, 0),
      feat_equity    = if_else(u_feature == "Equity"    | v_feature == "Equity",    1, 0),
      feat_market    = if_else(u_feature == "Market"    | v_feature == "Market",    1, 0),
      feat_risk      = if_else(u_feature == "Risk"      | v_feature == "Risk",      1, 0),
      delay          = rbinom(n(), 1, 0.2),
      edge_id        = row_number()
    ) %>%
    mutate(across(all_of(feat_cols), as.numeric))
  
  # Return both E_candidates and V so callers don't depend on global V
  list(
    E_candidates = E_candidates,
    V            = components, 
    feat_cols    = feat_cols
  )
}

# Get the evidence scores: three different types of evidence distribution ------
weak_E <- function(E_candidates){
  tibble(edge_id = E_candidates$edge_id, E = rnorm(nrow(E_candidates), -1, 1)) %>%
    mutate(E = pmax(pmin(E, 2), -2))
}

moderate_E <- function(E_candidates){
  tibble(edge_id = E_candidates$edge_id, E = rnorm(nrow(E_candidates), 0, 1)) %>%
    mutate(E = pmax(pmin(E, 2), -2))
}

strong_E <- function(E_candidates){
  tibble(edge_id = E_candidates$edge_id, E = rnorm(nrow(E_candidates), 1, 1)) %>%
    mutate(E = pmax(pmin(E, 2), -2))
}

E_score <- moderate_E

# Build psi vector with intercept for an edge row ------------------------------
# this gets the characteristics of each edge so the agents can react to them
psi_of <- function(row, feat_cols){
  vals <- unlist(row[, feat_cols, drop = FALSE], use.names = FALSE)
  out <-c(`(Intercept)` = 1, setNames(as.numeric(vals), feat_cols))
  out
}