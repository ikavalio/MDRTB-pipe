library(igraph)

# E1(x) = |x - p|
frn.phi1 <- function(x, p) abs(rep(x, length(p)) - p) # p - probability of relevance

# E2(x, y) = l * |cor(x, Y)| * I(x != y)
frn.phi2 <- function(r, lmbda) 2 * lmbda * abs(r)

frn.find.significant <- function(mutations, corr.matrix, lambda = 0.2, corr.threshold = 0.7) {
  R <- corr.matrix
  R[abs(R) < corr.threshold] <- 0
  
  E1 <- frn.phi1(1, mutations$relevance) # energy for X = 1
  E0 <- frn.phi1(0, mutations$relevance) # energy for X = 0
  E2 <- frn.phi2(R, lambda)
  
  adjM.rows <- c(rownames(R), "s", "t")
  adjM.cols <- c(rownames(R), "s", "t")
  adjM <- matrix(rep(0, length(adjM.rows) * length(adjM.cols)), 
                 nrow = length(adjM.rows), ncol = length(adjM.cols), 
                 dimnames = list(adjM.rows, adjM.cols))
  
  energComparison <- E0 < E1
  energValue <- abs(E1 - E0)
  
  # connect significant mutations with 's' and not significant with 't'
  adjM["s",] <- c(ifelse(energComparison, energValue, 0), 0, 0) # last 2 els are "s" and "t"
  adjM[,"t"] <- c(ifelse(energComparison, 0, energValue), 0, 0) # last 2 els are "s" and "t"
  #E2[lower.tri(E2, diag = TRUE)] <- 0 # use only one edge between any two variables
  diag(E2) <- 0
  adjM[1:ncol(snps.cor), 1:ncol(snps.cor)] <- E2
  
  frn <- graph.adjacency(adjM, weighted = TRUE, mode = "directed")
  #cut <- stMincuts(frn, source = "s", target = "t")
  E(frn)$capacity <- E(frn)$weight
  flow <- graph.maxflow(frn, "s", "t")
  setdiff(V(frn)$name[flow$partition2], "t")
}