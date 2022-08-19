rdirichlet.m <- function (alpha) {
  Gam <- matrix(rgamma(length(alpha), shape = alpha), nrow(alpha), ncol(alpha))
  t(t(Gam) / colSums(Gam))
}
EstPara0 <- function (otu.tab) {
  
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste0('OTU', 1 : nrow(otu.tab))
  } # otu * sample
  samplenames = colnames(otu.tab)
  taxnames = rownames(otu.tab)
  
  dirmult.paras <- dirmult::dirmult(t(otu.tab))
  
  gamma = dirmult.paras$gamma
  names(gamma) = names(dirmult.paras$pi)
  
  # Add pseduo count(each OTU add gamma estimated from dirmult)
  otu.tab = sapply(1:ncol(otu.tab), function (i) gamma + otu.tab[,i]) # C_ij otu * sample
  
  # back to dirchlet, calculate the true proportion
  otu.tab.p <- rdirichlet.m(otu.tab) # P_ij nOTU*nSam
  colnames(otu.tab.p) = samplenames
  rownames(otu.tab.p) = taxnames
  
  # order OTUs by mean OTU proportion, for later selection
  ord = order(rowMeans(otu.tab.p), decreasing = TRUE)
  otu.tab.p =  otu.tab.p[ord,]
  
  # apply size factor
  Si = exp(rnorm(ncol(otu.tab.p)))
  otu.tab0 = t(t(otu.tab.p)*Si)
  colnames(otu.tab0) = colnames(otu.tab.p)
  return(list(mu = otu.tab.p, otu.tab = otu.tab0))
}


EstPara <- function (otu.tab, dirmult.paras) {
  
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste0('OTU', 1 : nrow(otu.tab))
  } # otu * sample
  
  samplenames = colnames(otu.tab)
  taxnames = rownames(otu.tab)
  
  # For saving time, here will use the parameters already estimated
  gamma = dirmult.paras$gamma
  names(gamma) = names(dirmult.paras$pi)
  gamma = gamma[taxnames]
  
  # Add pseduo count(each OTU add gamma estimated from dirmult)
  otu.tab = sapply(1:ncol(otu.tab), function (i) gamma + otu.tab[,i]) # C_ij otu * sample
  
  # back to dirchlet, calculate the true proportion
  otu.tab.p <- rdirichlet.m(otu.tab) # P_ij nOTU*nSam
  colnames(otu.tab.p) = samplenames
  rownames(otu.tab.p) = taxnames
  
  # order OTUs by mean OTU proportion, for later selection
  ord = order(rowMeans(otu.tab.p), decreasing = TRUE)
  otu.tab.p =  otu.tab.p[ord,]
  
  # apply size factor
  Si = exp(rnorm(ncol(otu.tab.p)))
  otu.tab0 = t(t(otu.tab.p)*Si)
  colnames(otu.tab0) = colnames(otu.tab.p)
  
  return(list(mu = otu.tab.p, otu.tab = otu.tab0))
}

## Little different from the general simulation function in cluster, since we want to make the nOTU and nSam is the same as the nrow/ncol of the input otu.tab
SimulateSeq <- function (
  # Input data
  otu.tab, model.paras = NULL, 
  # OTU filtering 
  # prev = 0.05,  # Model selection
  model = c('loglinear','nonmonotone','heterogeneity'),
  # Count generation modeL
  # nSam = ncol(otu.tab), # to make sure the same sample size as the raw data in NULL comparison
  nSam =  c('nSam_L1','nSam_L2', 'nSam_L3','nSam_L4','nSam_L5'),
  nOTU = c('nOTU_L1','nOTU_L2','nOTU_L3','nOTU_L4','nOTU_L5'), 
  nSub = c('Sub_L1','Sub_L2','Sub_L3'),
  # Signal structure
  diff.otu.pct = c('low', 'medium', 'high'),
  diff.otu.direct = c('balanced', 'unbalanced'), 
  diff.otu.mode = c('mix', 'abundant', 'rare'), 
  include.top.otu = FALSE, 
  k.top.otu = 5, 
  # Covariate effect
  covariate.type = c("continuous", "binary"),  
  covariate.eff.mean =c('none','L1','L2','L3','L4','L5'),
  grp.ratio = 1,covariate.eff.sd = 0, weight.abund.func = function (x) 0,
  # Confounder effect
  confounder.type = c('none', 'continuous', 'binary', 'both'), 
  conf.cov.cor = 0.6, conf.diff.otu.pct = 0.05, conf.nondiff.otu.pct = 0.1, confounder.eff.mean = 0, confounder.eff.sd = 0,
  # Sequence depth
  depth.mu = c('D1','D2','D3','D4','D5'), depth.theta = 5,  depth.conf.factor = c('none','DL1','DL2','DL3','DL4','DL5')
  
){
  
  model <- match.arg(model)
  nSam <- match.arg(nSam)
  nOTU <- match.arg(nOTU)
  nSub <- match.arg(nSub)
  diff.otu.pct <- match.arg(diff.otu.pct)
  diff.otu.direct <- match.arg(diff.otu.direct)
  diff.otu.mode <- match.arg(diff.otu.mode)
  covariate.type  <- match.arg(covariate.type)
  covariate.eff.mean <- match.arg(covariate.eff.mean)
  confounder.type <- match.arg(confounder.type)
  depth.mu <- match.arg(depth.mu)
  depth.conf.factor <- match.arg(depth.conf.factor)
  
  # Adjusted for simulation purpose
  nSams = c(50, 100, 150, 200, 250)
  names(nSams) = c('nSam_L1','nSam_L2', 'nSam_L3','nSam_L4','nSam_L5')
  
  nOTUs = c(50, 100, 200, 300, 500)
  # nOTUs = c(10, 30, 50, 100, 500) 
  names(nOTUs) = c('nOTU_L1','nOTU_L2','nOTU_L3','nOTU_L4','nOTU_L5')
  
  nSubs = seq(0.1, 0.9, by = 0.3)
  names(nSubs) = c('Sub_L1','Sub_L2','Sub_L3')
  
  covariate.eff.means = c(0, 0.4, 0.8, 1.2, 1.6, 2) #c(0, 0.6, 1.3, 1.7, 2.2, 2.6)#
  names(covariate.eff.means) = c('none','L1','L2','L3','L4','L5') 
  
  diff.otu.pcts = c(0, 0.05, 0.1, 0.2)
  names(diff.otu.pcts) = c('none','low','medium','high')
  
  depth.conf.factors = c(0, 0.7, 0.9, 1.1, 1.3, 1.5)
  names(depth.conf.factors) = c('none','DL1','DL2','DL3','DL4','DL5')
  
  depth.mus = c(500, 1000, 10000, 50000, 100000)
  names(depth.mus) = c('D1','D2','D3','D4','D5')
  
  adj.nSam <- c('nSam_L1' = 1, 'nSam_L2' = 1, 'nSam_L3' = 1, 'nSam_L4' = 1, 'nSam_L5' = 1)
  adj.density <- c('none' = 1, 'low' = 1, 'medium' = 1, 'high' = 1)
  adj.abundance = c('rare'= 1, 'mix'= 1, 'abundant'= 1)
  
  
  # main part
  nOTU <- nOTUs[nOTU]
  nSam <- nSams[nSam]
  nSub <- nSubs[nSub]
  covariate.eff.mean = covariate.eff.means[covariate.eff.mean]
  diff.otu.pct = diff.otu.pcts[diff.otu.pct]
  depth.conf.factor = depth.conf.factors[depth.conf.factor]
  depth.mu = depth.mus[depth.mu]
  
  #load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/VaginalIntroitus_V35_dirmult.RData')
  #load('~/Documents/Mayo_Research/SemiSimulation/DAMethodsEvaluation/data/SimulationEvaluation/VaginalIntroitus_V35.RData')
  otu.tab0 = otu.tab

  ## Estimated parameters
  if (is.null(model.paras)) {
    model.paras <- EstPara(otu.tab, dirmult.paras = dirmult.paras)
  } else {
    model.paras <- EstPara0(otu.tab = otu.tab) # full otu.tab + already estimated parameters
  }
  
  
  ## ADD 07/08/2021, rare OTU does not have power, debug here
  idx.sample <- sample(colnames(model.paras$otu.tab), nSam, replace = T)
  otu.tab <- model.paras$otu.tab[,idx.sample]
  idx.otu <- names(rowMeans(otu.tab)[order(-rowMeans(otu.tab))][1:nOTU])
  otu.tab <- otu.tab[idx.otu,]
  non.otu <- rownames(model.paras$otu.tab)[!(rownames(model.paras$otu.tab) %in% idx.otu)]
  non.sam <- colnames(model.paras$otu.tab)[!(colnames(model.paras$otu.tab) %in% idx.sample)]
  otu.tab.unselect = otu.tab0[non.otu,][,non.sam] # For overfitting test
  
  
  # User select random nSam and top nOTU
  # otu.names <- colnames(model.paras$otu.tab)
  # otu.tab = model.paras$otu.tab[(1 : (nOTU)),] # otu * sample
  # idx.otu = rownames(otu.tab)
  # idx.sample <- sample(otu.names, nSam)
  # idx.nonsample <- colnames(otu.tab)[!(colnames(otu.tab) %in% idx.sample)]
  # otu.tab = otu.tab[,idx.sample] # otu * sample, !!! This is the final selected OTU tab by User
  # otu.tab.unselect = otu.tab0[c(1 : (nOTU)),][,idx.nonsample] # For overfitting test
  
  # Generate subgroup 
  if(model =='heterogeneity'){
    nOTU0 = nOTU
    nOTU = nOTU * nSub
    idx.subgrp = sample(rownames(otu.tab), nOTU)
    otu.tab.sub = otu.tab[idx.subgrp,]
    otu.tab.nonsub = otu.tab[setdiff(rownames(otu.tab), idx.subgrp),]
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
    X <- cbind(ifelse(X <= quantile(X, grp.ratio / (1 + grp.ratio)), 0, 1))
  }
  
  # Generate log fold change for covariate effect
  covariate.eff.mean1 = covariate.eff.mean 
  covariate.eff.mean2 = covariate.eff.mean *  adj.density[names(diff.otu.pct)] * adj.abundance[diff.otu.mode] * adj.nSam[names(nSam)]
  
  if (diff.otu.direct == "balanced") {
    if (diff.otu.mode == 'abundant'){
      eta.diff <- sample(c(rnorm(floor(nOTU / 2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                           rnorm(nOTU - floor(nOTU / 2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))  %*% t(scale(X)) # multiple means * xb, eta = xb
    }else if(diff.otu.mode == 'rare'){
      eta.diff <- sample(c(rnorm(floor(nOTU / 2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                           rnorm(nOTU - floor(nOTU / 2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))  %*% t(scale(X))
    }else{
      eta.diff <- c(sample(c(rnorm(floor(nOTU / 4), mean = -covariate.eff.mean1, sd = covariate.eff.sd), 
                             rnorm(floor(nOTU / 2) - floor(nOTU / 4), mean = covariate.eff.mean1, sd = covariate.eff.sd))),
                    sample(c(rnorm(floor((nOTU -floor(nOTU / 2))/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd), 
                             rnorm(nOTU -floor(nOTU / 2) - floor((nOTU -floor(nOTU / 2))/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*% t(scale(X))
    }
  }
  
  
  if (diff.otu.direct == "unbalanced") {
    if (diff.otu.mode == 'abundant'){
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% t(scale(X)) 
    }else if(diff.otu.mode == 'rare'){
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*% t(scale(X)) 
    }else{
      eta.diff <- c(sample(c(rnorm(floor(nOTU / 2), mean = covariate.eff.mean1, sd = covariate.eff.sd))), 
                    sample(c(rnorm(nOTU - floor(nOTU / 2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*% t(scale(X))
    }
  }
  
  # Generate log fold change for confounder effect - assume balanced change
  eta.conf <-  sample(c(rnorm(floor(nOTU / 2), mean = -confounder.eff.mean, sd = confounder.eff.sd), 
                        rnorm(nOTU - floor(nOTU / 2), mean = confounder.eff.mean, sd = confounder.eff.sd)))  %*% t(scale(scale(Z) %*% rep(1, ncol(Z)))) 
  
  # Get the differential OTU id
  otu.ord <- 1 : (nOTU)
  diff.otu.ind <- NULL
  diff.otu.num <- round(diff.otu.pct * nOTU)
  if (include.top.otu == TRUE) {
    # Create strong compositional effecy by including k most abundant OTUs
    diff.otu.ind <- c(diff.otu.ind, 1 : k.top.otu)
    otu.ord <- setdiff(otu.ord, diff.otu.ind)
    diff.otu.num <- diff.otu.num - k.top.otu
    eta.diff[diff.otu.ind,] = abs(eta.diff[diff.otu.ind,]) # Add 11/18/2020 for same direction of top OTUs
  }
  
  if (diff.otu.mode == 'mix') diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[1 : round(length(otu.ord) / 4)], floor(diff.otu.num/2)),
                                                sample(otu.ord[round(3 * length(otu.ord) / 4) : length(otu.ord)], diff.otu.num-floor(diff.otu.num/2)))
  if (diff.otu.mode == 'abundant') diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[1 : round(length(otu.ord) / 4)], diff.otu.num))
  if (diff.otu.mode == 'rare') diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[round(3 * length(otu.ord) / 4) : length(otu.ord)], diff.otu.num))
  
  # Get the confounder OTU id
  conf.otu.ind <- c(sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct)), sample(setdiff(1 : (nOTU), diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))
  
  eta.diff[setdiff(1 : (nOTU), diff.otu.ind), ] <- 0
  eta.conf[setdiff(1 : (nOTU), conf.otu.ind), ] <- 0
  
  
  # model selection
  if(model =='loglinear'){
    # eta.exp <- exp(t(mu + eta.diff + eta.conf))
    eta.exp <- exp(t(eta.diff + eta.conf))
    eta.exp <- eta.exp * t(otu.tab)
  }
  
  if(model =='nonmonotone'){
    eta = eta.diff + eta.conf
    eta.exp <- dnorm(t(eta))
    eta.exp <- 3 * eta.exp
    # eta.exp[,setdiff(1 : (nOTU), diff.otu.ind)] = dnorm(0) *3
    # eta.exp[,setdiff(1 : (nOTU), conf.otu.ind)] = dnorm(0) *3
    eta.exp <-  eta.exp * t(otu.tab)
    # eta.exp <- dnorm(t(eta))
    # eta.exp <- t(mu) * eta.exp * 2 # 5000* 100 * t(100* 5000)
    # eta.exp <- 1/((t(eta))^2 +0.05)
  }
  
  if(model =='heterogeneity'){
    eta.exp <- exp(t(eta.diff + eta.conf)) # nSam * notu
    eta.exp <- eta.exp/rowSums(eta.exp)
    eta.exp <- eta.exp * colSums(otu.tab.sub) # 400*100 no rownames + colnames
    colnames(eta.exp) = rownames(otu.tab.sub)
    rownames(eta.exp) = colnames(otu.tab.sub)
    eta.exp <- cbind(eta.exp, t(otu.tab.nonsub)) # cbind(100*400, 100* 600) = 100*1000
    eta.exp <- eta.exp[,rownames(otu.tab)]
  }
  
  # Renormalize
  # make sure the input is sam * otu
  otu.tab.prop <- eta.exp / rowSums(eta.exp)
  otu.tab.prop <- t(otu.tab.prop) # otu * sample
  
  # Generate the sequence depth
  nSeq <- rnegbin(nSam, mu = depth.mu * exp(scale(X) * depth.conf.factor), theta = depth.theta)
  otu.tab.sim <- sapply(1:ncol(otu.tab.prop), function (i) rmultinom(1, nSeq[i], otu.tab.prop[,i]))
  
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- colnames(eta.exp)
  
  if(model =='heterogeneity'){
    diff.otu.names = rownames(otu.tab.sub)[diff.otu.ind]
    diff.otu.ind = which(rownames(otu.tab.sim) %in% diff.otu.names)
    conf.otu.names = rownames(otu.tab.sub)[conf.otu.ind]
    conf.otu.ind = which(rownames(otu.tab.sim) %in% conf.otu.names)
    diff.otu.ind = (1 : nOTU0) %in% diff.otu.ind
    conf.otu.ind = (1 : nOTU0) %in% conf.otu.ind
  }else{
    diff.otu.ind = (1 : nOTU) %in% diff.otu.ind
    conf.otu.ind = (1 : nOTU) %in% conf.otu.ind
  }
  
  return(list(call = match.call(), model.paras = model.paras, 
              otu.names = idx.otu, idx.sample = idx.sample, 
              otu.tab.sim = otu.tab.sim, otu.tab.unselect = otu.tab.unselect,
              nSeq = nSeq, otu.tab.prop = otu.tab.prop, X = X,  Z = Z, 
              conf.otu.ind = conf.otu.ind, diff.otu.ind = diff.otu.ind))
}
