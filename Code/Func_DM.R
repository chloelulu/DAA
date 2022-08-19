require(MASS)
require(dirmult)

rdirichlet.m <- function (alpha) {
  Gam <- matrix(rgamma(length(alpha), shape = alpha), nrow(alpha), ncol(alpha)) 
  t(t(Gam) / colSums(Gam))
}

EstParaDM <- function (otu.tab, ref.pct = 0.05) {
  # mu = log(pi_k / pi_K+1) = beta_k0 + sum(beta_kj X_j) 
  # mu on the logit scale, phi - overdipsersion parameter
  
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste0('OTU', 1 : nrow(otu.tab))
  }
  
  otu.tab <- t(otu.tab)
  
  dirmult.obj <- dirmult(otu.tab)
  phi <- dirmult.obj$theta
  pi <- dirmult.obj$pi
  
  ord <- order(pi, decreasing = TRUE)
  otu.tab <- otu.tab[, ord]
  pi <- pi[ord]
  
  ref.ind <- max(1, floor(ref.pct * ncol(otu.tab)))
  ref.otu <- names(pi)[ref.ind]
  
  pi <- pi[-ref.ind]
  ref <- pi[ref.ind]
  
  mu <- log(pi / ref)
  otu.tab <- t(cbind(otu.tab[, -ref.ind], otu.tab[, ref.ind]))
  rownames(otu.tab) <- c(names(pi), ref.otu)
  
  return(list(mu = mu, V = phi, ref.otu = ref.otu, otu.tab = otu.tab))
}




SimulateMSeq <- function (
  # Input data
  otu.tab, otu.tree = NULL, 
  
  # Count generation model
  nSam = ncol(otu.tab), nOTU = nrow(otu.tab), 
  model.paras = NULL, 
  
  # Signal structure
  diff.otu.pct = 0,  diff.otu.direct = c('balanced', 'unbalanced'), 
  diff.otu.mode = c('random', 'abundant', 'rare', 'tree'), tree.prior = 0.5,
  include.top.otu = FALSE, k.top.otu = 2, 
  
  # Covariate effect
  covariate.type = c("continuous", "binary"),  grp.ratio = 1,  #covariate.eff.mean = 1,
  covariate.eff.mean = 0,
  covariate.eff.sd = 0, weight.abund.func = function (x) 0,
  
  # grp.ratio = 1;covariate.eff.mean = 0;covariate.eff.sd = 0; weight.abund.func = function (x) 0;
  # conf.cov.cor = 0.6;
  # conf.diff.otu.pct = 0.05; conf.nondiff.otu.pct = 0.1; confounder.eff.mean = 0;
  # confounder.eff.sd = 0
  # depth.theta = 5;  depth.conf.factor = 0;
  # outlier.prop = 0; otulier.change.func = function (x) sqrt(x);
  # zinfl.otu.pct = 0; zinfl.pi.low = 0; zinfl.pi.high = 0.8
  
  # Confounder effect
  confounder.type = c('none', 'continuous', 'binary', 'both'), conf.cov.cor = 0.6,
  conf.diff.otu.pct = 0.05, conf.nondiff.otu.pct = 0.1, confounder.eff.mean = 0, 
  confounder.eff.sd = 0,
  
  # Sequence depth
  depth.mu = 10000, depth.theta = 5,  depth.conf.factor = 0,
  
  # Outlier 
  outlier.prop = 0, otulier.change.func = function (x) sqrt(x),
  
  # Zero-inflation
  zinfl.otu.pct = 0, zinfl.pi.low = 0, zinfl.pi.high = 0.8
) {
  
  diff.otu.mode <- match.arg(diff.otu.mode)
  diff.otu.direct <- match.arg(diff.otu.direct)
  covariate.type  <- match.arg(covariate.type)
  confounder.type <- match.arg(confounder.type)
  
  if (is.null(model.paras)) {
    if (nOTU > nrow(otu.tab)) {
      stop('The OTU number requested exceeds the maximum!\n')
    } 
    otu.tab <- t(otu.tab)
    otu.tab.p <- otu.tab / rowSums(otu.tab)
    otu.tab <- t(otu.tab[, order(colMeans(otu.tab.p), decreasing = TRUE)[1 : nOTU]])
    model.paras <- EstParaDM(otu.tab)

  } else {	
    if (nOTU > length(model.paras$mu) + 1) {
      stop('The OTU number requested exceeds the maximum!\n')
    } 
    
    if (nOTU <  length(model.paras$mu) + 1) {
      warning('Some rare OTUs will be dropped to achieve the requested OTU number!\n')
      model.paras$mu <- model.paras$mu[1 : (nOTU - 1)]
    }
  }
  
  otu.names <- c(names(model.paras$mu), model.paras$ref.otu)
  
  if (!is.null(otu.tree)) {
    if (sum(!(otu.names %in% otu.tree$tip.label)) != 0) {
      stop("The OTU table contains unknown OTUs! OTU names\n\t\t\t\t\tin the OTU table and the tree should match!")
    }
    absent <- otu.tree$tip.label[!(otu.tree$tip.label %in% otu.names)]
    otu.tree <- drop.tip(otu.tree, absent)
  } else {
    if (diff.otu.mode == 'tree') stop('Selecting OTUs by tree structure requires an OTU tree!\n')
  }
  
  # Generate confounder
  if (confounder.type == "none")  {
    confounder.type <- "continuous"
    confounder.eff.mean <- 0
  } 
  if (confounder.type == "continuous") Z <- cbind(rnorm(nSam))
  if (confounder.type == "binary") Z <- cbind(c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2)))
  if (confounder.type == "both") Z <- cbind(rnorm(nSam), c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2)))
  
  # Generate covariate - assume equal confounding effects for "both"
  rho <- sqrt(conf.cov.cor ^ 2 / (1 - conf.cov.cor ^ 2))
  
  if (covariate.type == 'continuous') {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z)))  + rnorm(nSam)
  }
  
  if (covariate.type == "binary") { 
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z)))  + rnorm(nSam)
    # Will not use logistic model
    X <- cbind(ifelse(X <= quantile(X, grp.ratio / (1 + grp.ratio)), 0, 1))
  }
  
  # Generate log fold change for covariate effect
  if (diff.otu.direct == "balanced") {
    eta.diff <- sample(c(rnorm(floor(nOTU / 2), mean = -covariate.eff.mean, sd = covariate.eff.sd), 
                         rnorm(nOTU - floor(nOTU / 2) - 1, mean = covariate.eff.mean, sd = covariate.eff.sd)))  %*% t(scale(X))
    
  }
  
  if (diff.otu.direct == "unbalanced") {
    eta.diff <- rnorm(nOTU - 1, mean = covariate.eff.mean, sd = covariate.eff.sd) %*% t(scale(X))		
  }
  
  # Generate log fold change for confounder effect - assume balanced change
  eta.conf <-  sample(c(rnorm(floor(nOTU / 2), mean = -confounder.eff.mean, sd = confounder.eff.sd), 
                        rnorm(nOTU - floor(nOTU / 2) - 1, mean = confounder.eff.mean, sd = confounder.eff.sd)))  %*% t(scale(scale(Z) %*% rep(1, ncol(Z)))) 
  
  # Get the differential OTU id
  otu.ord <- 1 : (nOTU - 1)
  diff.otu.ind <- NULL
  diff.otu.num <- round(diff.otu.pct * nOTU)
  if (include.top.otu == TRUE) {
    # Create strong compositional effecy by including k most abundant OTUs
    diff.otu.ind <- c(diff.otu.ind, 1 : k.top.otu)
    otu.ord <- setdiff(otu.ord, diff.otu.ind)
    diff.otu.num <- diff.otu.num - k.top.otu
  }
  
  if (diff.otu.mode == 'random') diff.otu.ind <- c(diff.otu.ind, sample(otu.ord, diff.otu.num))
  if (diff.otu.mode == 'abundant') diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[1 : round(length(otu.ord) / 4)], diff.otu.num))
  if (diff.otu.mode == 'rare') diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[round(3 * length(otu.ord) / 4) : length(otu.ord)], diff.otu.num))
  
  # Sampling according to the tree, may be revised in the future
  if (diff.otu.mode == 'tree') {
    D <- cophenetic(otu.tree)
    D <- D[otu.names, otu.names]
    otu.tree.ord <- setdiff(order(abs(mvrnorm(1, mu = rep(0, ncol(D)), Sigma = exp(-D * tree.prior))), decreasing = TRUE), c(diff.otu.ind, nOTU))
    diff.otu.ind <- c(diff.otu.ind, otu.tree.ord[1 : diff.otu.num])
  }
  
  # Get the confounder OTU id
  if(covariate.eff.mean!=0){
    conf.otu.ind <- c(sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct)), 
                      sample(setdiff(1 : (nOTU - 1), diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))
    eta.conf[setdiff(1 : (nOTU - 1), conf.otu.ind), ] <- 0
  }
  eta.diff[setdiff(1 : (nOTU - 1), diff.otu.ind), ] <- 0
  
  # Generate zero-index for zeroinflated model
  zero.ind.mat <- matrix(1, nSam, nOTU - 1)
  zero.col.ind <- sample(1 : (nOTU - 1), round(zinfl.otu.pct * nOTU))
  zero.col.pct <- NULL
  for (ind in zero.col.ind) {
    zero.pct <- runif(1, zinfl.pi.low, zinfl.pi.high)
    zero.num <- round(zero.pct * nSam)
    zero.col.pct <- c(zero.col.pct, zero.pct)
    zero.ind.mat[sample(1 : nSam, zero.num), ind] <- 0
  }
  
  # Generate the proportion
  eta.exp <- exp(t(model.paras$mu + eta.diff + eta.conf))
  eta.exp <- eta.exp * zero.ind.mat
  eta.exp.sum <- rowSums(eta.exp) + 1
  otu.tab.prop <- cbind(eta.exp / eta.exp.sum, 1 / eta.exp.sum)
  otu.tab.gamma <- otu.tab.prop * (1 - model.paras$V) / model.paras$V
  otu.tab.prop <- t(rdirichlet.m(t(otu.tab.gamma)))

  # Add outliers
  outlier.ind <- sample(which(otu.tab.prop != 0), floor(outlier.prop * sum(otu.tab.prop != 0)))
  otu.tab.prop[outlier.ind] <- otulier.change.func(otu.tab.prop[outlier.ind])
  
  # Renormalize
  otu.tab.prop <- otu.tab.prop / rowSums(otu.tab.prop)
  otu.tab.prop <- t(otu.tab.prop)
  
  # Generate the sequence depth
  nSeq <- rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta)
  
  otu.tab.sim <- sapply(1:ncol(otu.tab.prop), function (i) rmultinom(1, nSeq[i], otu.tab.prop[, i]))
  colnames(otu.tab.sim) <- paste0('S', 1 : ncol(otu.tab.sim))
  rownames(otu.tab.sim) <- otu.names
  
  return(list(call = match.call(), otu.tree = otu.tree, model.paras = model.paras, otu.tab.sim = otu.tab.sim,
              nSeq = nSeq, otu.tab.prop = otu.tab.prop, otu.names = otu.names, X = X, Z = Z,
              #diff.otu.ind = (1 : nOTU) %in% diff.otu.ind, conf.otu.ind = (1 : nOTU) %in% conf.otu.ind,
              zeroinfl.ind = zero.col.ind, zeroinfl.pct = zero.col.pct,
              outlier.ind = outlier.ind))
}

