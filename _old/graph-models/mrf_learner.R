# simulation
# variables init Markov chain: P(v1=1)=0.7 v1 [0.7]-> v2 [0.2]-> v3 [0.1]-> v4 [0.6]-> v5

train.init.p <- 0.7
train.prob <- c(0.7, 0.2, 0.1, 0.6)
train.dsize <- 1000

theta.beta <- 0.7
theta.noise <- 0.2
theta.thresh <- 1e-3
data.dim <- length(chain.prob) + 1
params.tot <- data.dim ^ 2
cd.its <- 100
mcmc.n <- 5
learn.rate <- 0.1
persistent <- TRUE

init.chain.sample <- function(p.init.one, chain.prob) {
  init.state <- rbinom(1, 1, p.init.one)
  Reduce(function(l, r) c(l, ifelse(rbinom(1, 1, r) == 1, l, 1 - l)), chain.prob, init.state)
}

input.data <- as.matrix(
  Reduce(function(l, r) 
    rbind(l, init.chain.sample(train.init.p, train.prob)), 
    1:train.dsize, 
    data.frame()
  ))

theta.init <- function(n, beta = 0.7, noise.sd = 0.2) {
  stopifnot(beta != 1)
  src <- rnorm(n, mean = log(beta / (1 - beta)), sd = noise.sd)
}

theta.graph.adjust <- function(theta, ext, edge.bidirectional = FALSE) {
  theta <- matrix(theta, ncol = ext, nrow = ext)
  diag(theta) <- 0 # no edges to itself
  if(!edge.bidirectional) {
    theta[lower.tri(theta, FALSE)] <- t(theta)[lower.tri(theta, FALSE)]
  }
  as.numeric(theta)
}

theta.fix <- function(theta, threshold) {
  theta[abs(theta) < threshold] <- 0
  theta
}

f <- function(x, theta) {
  exp(theta %*% x)
}

phi <- function(x) {
  as.integer(sapply(x, function(el) el == x))
}

calc.alpha <- function(x.prop, x.curr, theta) {
  exp(theta %*% (phi(x.prop) - phi(x.curr)))
}

e.estimate <- function(smpls) {
  phi.m <- apply(smpls, 1, phi)
  apply(phi.m, 1, mean)
}

.next.mutate <- function(x) {
  t.index <- sample(1:length(x), 1)
  x[t.index] <- 1 - x[t.index]
  x
}

.next.swap <- function(x) {
  t.indexes <- sample(1:length(x), 2, replace = FALSE)
  x[rev(t.indexes)] <- x[t.indexes]
  x
}

.next.reverse <- function(x) {
  rev(x)
}

.next.reverse.halves <- function(x) {
  sz <- length(x)
  hv <- floor(sz / 2) + 1
  c(x[hv:sz], x[1:hv - 1])
}

next.sample <- function(x) {
  # randomly select binary coordinate and flip it
  actions <- c("mutate", "swap.els", "swap.halves", "reverse")
  f <- switch(sample(actions, 1),
              mutate = .next.mutate,
              swap.els = .next.swap,
              swap.halves = .next.reverse.halves,
              reverse = .next.reverse 
  )
  f(x)
}

mcmc.sample <- function(x, theta, n = 10) {
  x.hat <- x
  for(i in 1:n) {
    next.accepted <- FALSE
    try.i <- 0
    next.proposal <- next.sample(x.hat)
    alpha <- calc.alpha(next.proposal, x.hat, theta)
    if( as.logical(ifelse(alpha < 1, rbinom(1, 1, alpha), alpha)) ) x.hat <- next.proposal 
  }
  x.hat
}

cd.run <- function(input.data,
                   params.tot,
                   steps = 10000,
                   mcmc.n = 5, 
                   beta = 0.7,
                   theta.noise = 0.2,
                   theta.thresh = 1e-3,
                   learn.rate = 0.2,
                   persistent = TRUE) {
  theta.curr <- theta.graph.adjust(theta.init(params.tot, beta, theta.noise), ncol(input.data))
  e.sample <- e.estimate(input.data)
  for(i in 1:steps) {
    x.chain.source <- if(persistent) x.upd else input.data
    x.upd <- t(apply(x.chain.source, 1, mcmc.sample, theta = theta.curr, n = mcmc.n))
    e.upd <- e.estimate(x.upd)
    theta.curr <- theta.fix(theta.curr + learn.rate * (e.sample - e.upd), theta.thresh)
    # adjust learning rate
    if(i %% 10 == 0) cat(sprintf("Iteration %d done...\n", i))
  }
  theta.curr
}

theta.hat <- cd.run(input.data,
                    params.tot,
                    steps = cd.its,
                    mcmc.n = mcmc.n,
                    beta = theta.beta,
                    theta.noise = theta.noise,
                    theta.thresh = theta.thresh,
                    learn.rate = learn.rate,
                    persistent = persistent)
