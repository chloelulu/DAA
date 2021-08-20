# 2019_02_20
# Beta-mixture to address sequence depth confounding
# Permutation-based FDR control
# Multiple-stage normalization
# Variance modulation by LIMMA - hasn't been finished
# Expectation of F statistics based on posterior sampling
# Omnibus test for combining diffrent effects
# Version 0.9 [PermuteDAA3.2]

require(matrixStats)
require(permute)
require(nlme)
require(vegan)
require(ape)
require(statmod)


# Please indicate the code is copied from vegan package for portability
getPermuteMatrix <- function (perm, N, strata = NULL) {
	if (length(perm) == 1) {
		perm <- how(nperm = perm)
	}
	if (!missing(strata) && !is.null(strata)) {
		if (inherits(perm, "how") && is.null(getBlocks(perm))) 
			setBlocks(perm) <- strata
	}
	if (inherits(perm, "how")) 
		perm <- shuffleSet(N, control = perm)
	if (is.null(attr(perm, "control"))) 
		attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"), 
						nperm = nrow(perm)), class = "how")
	perm
}

# Permutation-based FDR control
perm_fdr_adj <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	F0 <- F0[ord]
	perm.no <- ncol(Fp)
	Fp <- as.vector(Fp)
	Fp <- Fp[!is.na(Fp)]
	Fp <- sort(c(Fp, F0), decreasing = F)
	n <- length(Fp)
	m <- length(F0)
	FPN <- (n + 1) - match(F0, Fp) - 1:m
	p.adj.fdr <- FPN / perm.no / (1:m)
#		p.adj.fdr <- sapply(F0, function(x) sum(Fp >= 
#									x, na.rm=TRUE) / perm.no)/(1:length(F0))
	p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
}

# Permutation-based FWER control
# Conservative
perm_fwer_adj <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	F0 <- F0[ord]
	col.max <- colMaxs(Fp, na.rm=TRUE)
	p.adj.fwer <- sapply(F0, function(x) mean(col.max >= x))[order(ord)]
}


# Westfall young too slow
perm_fwer_adj2 <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	m <- length(F0)
	F0 <- F0[ord]
	Fp <- Fp[ord, , drop=FALSE]
	col.max <- Fp[m, ]
	p.adj.fwer <- sapply(m:1, function(i) {
				x <- F0[i]
				y <- Fp[i, ]
				col.max <<- ifelse(y > col.max, y, col.max)
				mean(col.max >= x)
			})
	p.adj.fwer <- rev(p.adj.fwer)
	p.adj.fwer <- pmin(1, rev(cummin(rev(p.adj.fwer))))[order(ord)]
}

# Pad back the NA values
na.pad <- function (vec, ind) {
	vec0 <- numeric(length(ind))
	vec0[!ind] <- vec
	vec0[ind] <- NA
	vec0
}

# Pad back the NA values
one.pad <- function (vec, ind) {
	vec0 <- numeric(length(ind))
	vec0[!ind] <- vec
	vec0[ind] <- 1
	vec0
}


find.ref <- function(otu.tab, ref.pct = 0.50) {
	p <- nrow(otu.tab)
	otu.tab <- otu.tab + 1
	res <- matrix(0, p, p)
	for (i in 1 : (p - 1)) {
		
		for (j in (i + 1) : p) {
			res[i, j] <- res[j, i] <- sd(log(otu.tab[i, ] / otu.tab[j, ]))
		}
	}
	ref.ind <- order(rowMedians(res))[1 : round(ref.pct * p)]
	return(ref.ind)
}

fastLM <- function(Y, M) {
	Y <- as.matrix(Y)
	XXI <- solve(t(M) %*% M)
	dof <- ncol(Y) - ncol(M)
	est <- XXI %*% t(M) %*% t(Y)
	resid <- t(Y) - M %*% est
	sigma <- sqrt(colSums(resid^2)/dof)
	pval <- 2 * pt(-abs(t(est/(sqrt(diag(XXI)))) / sigma), dof)
	return(list(pval = pval, sigma = sigma, dof = dof, est = est))
}

find.ref.z1 <- function(otu.tab, M,  ref.pct = 0.50) {
	p <- nrow(otu.tab)
	n <- ncol(otu.tab)
	otu.tab <- otu.tab + 1
	otu.tab <- t(otu.tab)
	res <- matrix(NA, p, p)
	
	XXI <- solve(t(M) %*% M)
	dof <- n - ncol(M)
	XXIM <- XXI %*% t(M)
	
	for (i in 1 : (p - 1)) {
		Y <- (as.matrix(log(otu.tab[, i] / otu.tab[, (i + 1) : p])))

		est <- XXIM %*% Y
		resid <- Y - M %*% est
		sigma <- sqrt(colSums(resid^2) / dof)
		
		res[i, (i + 1) : p] <- res[(i + 1) : p, i] <- sigma
	}
	ref.ind <- order(rowMedians(res, na.rm = TRUE))[1 : round(ref.pct * p)]
	return(ref.ind)
}

find.ref.z2 <- function(otu.tab, Q0, Q01, ref.pct = 0.50) {
	p <- nrow(otu.tab)
	n <- ncol(otu.tab)
	otu.tab <- otu.tab + 1
	otu.tab <- t(otu.tab)
	res <- matrix(0, p, p)

	for (i in 1 : (p - 1)) {
		Y <- t(as.matrix(log(otu.tab[, i] / otu.tab[, (i + 1) : p])))
		TSS <- rowSums(Y^2)
		MSS01 <- rowSums((Y %*% Q01)^2)
		MSS0 <- rowSums((Y %*% Q0)^2)
		
		MSS <- (MSS01 - MSS0) 
		RSS <- (TSS - MSS01) 
		
		res[i, (i + 1) : p] <- res[(i + 1) : p, i] <- MSS / RSS
	}
	ref.ind <- order(rowMedians(res))[1 : round(ref.pct * p)]
	return(ref.ind)
}

find.ref.z3 <- function(otu.tab, M,  ref.pct = 0.50) {
	p <- nrow(otu.tab)
	n <- ncol(otu.tab)
	otu.tab <- t(otu.tab)
	res <- matrix(0, p, p)
	
	XXI <- solve(t(M) %*% M)
	dof <- n - ncol(M)
	XXIM <- XXI %*% t(M)
	
	for (i in 1 : (p - 1)) {
		Y <- (as.matrix(log(otu.tab[, i] / otu.tab[, (i + 1) : p])))
		
		est <- XXIM %*% Y
		resid <- Y - M %*% est
		sigma <- sqrt(colSums(resid^2) / dof)
		
		res[i, (i + 1) : p] <- res[(i + 1) : p, i] <- sigma
	}
	ref.ind <- order(rowMeans(res))[1 : round(ref.pct * p)]
	return(ref.ind)
}

# Permutation-based differential abundance analysis
ZicoSeq.plus <- function (
		meta.dat, comm, grp.name, adj.name = NULL, 
		prev.filter = 0.1, abund.filter = 10, nonzero.filter = 5, min.prop = 0, 
		is.winsor = TRUE, winsor.qt = 0.97,
		is.prior = TRUE,  post.sample.no = 25, 
		link.func = list(function (x) x^0.5),
		perm.no = 99,  strata = NULL, stats.combine.func = max, 
		ref.pct = 0.50, stage.no = 5, excl.pct = 0.50,
		is.fwer = FALSE, verbose = TRUE, return.comm = FALSE, return.perm.F = FALSE, ...) {
	
	# Args:
	#   meta.dat: a data frame containing the sample information
	#   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
	#   grp.name: a character, variable of interest; it could be numeric or categorical; should be in "meta.dat"
	#   adj.name: a character vector, variable(s) to be adjusted; they could be numeric or categorical; should be in "meta.dat"
	#   prev.filter: features with prevalence (i.e., nonzero proportion) less than "prev.cutoff" will be filtered
	#   abund.filter: features with a total counts less than "abund.cutoff" will be filtered 
	#   nonzero.filter: features with nonzero OTUs less than "nonzero.filter" will be filtered
	#	is.winsor: a logical value indicating whether winsorization should be performed to replace outliers. The default is TRUE.		   
	#   winsor.qt: the winsorization quantile, above which the counts will be replaced
	#   is.prior: a logical value indicating whether to perform posterior inference based on some prior distribution on the proportion data
	#   prior.dist: prior distribution, either two-component beta-binomial mixture ("BetaMix") or zeroinflated beta-binomial ("ZIBB")
	#   post.method: method for posterior inference, either based on posterior sampling ("sample") or approximate posterior mean ("mean")
	#   post.sample.no:  the number of posterior samples if posterior sampling is used
	#   link.func:  a list of functions that connects the ratios to the covariates
	#   link.d.func: a list of the derivative function of "link.func"; only need to specifiy when "post.method" is "mean"
	#   variance.EB: a logical value indicating whehter to perform empirical Bayes based variance shrinkage
	#   df.prior: the degree of freedom of the prior inverse gamma distribution for variance shrinkage
	#   perm.no: the number of permutations; If the raw p values are of the major interest, set "perm.no" to at least 999
	#   strata:  a factor indicating the permutation strata; permutation will be confined to each stratum 
	#   stats.combine.func: function to combine the F-statistic for the omnibus test
	#   stage.no: the number of stages if multiple-stage ratio stategy is used
	#   topK: the number of dominant features that will be excluded in the initial stage ratio calculation
	#   stage.fdr: the fdr cutoff below which the features will be excluded for calculating the ratio
	#   stage.max.pct: the maximum percentage of features that will be excluded
	#   is.fwer: a logical value indicating whether the family-wise error rate control (West-Young) should be performed
	#   is.tree.fdr: a logical value indicating whether tree-based false discovery rate shuold be carried out
	#   tree:  a class of "phylo", the tree relats all the OTUs, and should have the same names in "comm"
	#   verbose: a logical value indicating whether the trace information should be printed out
	#   return.comm: a logical value indicating whether the wisorized, filtered "comm" matrix should be returned
	#   ...:  arguments passing to tree-based fdr control
	
	#
	# Returns:
	#   a list that contains:
	#      call: the call
	#      comm: the wisorized, filtered "comm" matrix
	#      filter.ind: a vector of logical values indicating which features are tested
	#      R2: a matrix of percent explained variance (number of features by number of functions)
	#      F0: a matrix of F-statistics (number of features by number of functions)
	#      RSS: a matrix of residual sum squares (number of features by number of functions)
	#      df.model, df.residual:  degree of freedoms for the model and residual space
	#      p.raw : the raw p-values based on permutations (not accurate if "perm.no" is small)
	#      p.adj.fdr: permutation-based FDR-adjusted p-values
	#      p.adj.tree.fdr: permutation-based tree FDR-adjusted p-values 
	#      p.adj.fwer: permutation-based FWER-adjusted (West-Young) p-values
	#      tree.fdr.obj: the object returned by the "TreeFDR"
	
	
	this.call <- match.call()
	###############################################################################
	# Winsorization to reduce the influence of outlier counts
	if (is.winsor == TRUE) {
		depth <- colSums(comm)
		comm <- apply(comm, 1, function (x) {
					x <- x / depth
					qt <- quantile(x, winsor.qt)
					x[x >= qt] <- qt
					round(x * depth)
				})
		comm <- t(comm)
		
	}
	###############################################################################
	# Filter to remove very rare taxa
	filter.ind <- rowMeans(comm != 0) >= prev.filter & rowSums(comm) >= abund.filter & rowSums(comm != 0) >= nonzero.filter
	names(filter.ind) <- rownames(comm)
	
	if (verbose)  cat(sum(!filter.ind), ' features are filtered!\n')
	
	comm <- comm[filter.ind, ]
	
	sample.no <- ncol(comm)
	otu.no <- nrow(comm)
	row.names <- rownames(comm)
	depth <- colSums(comm)
	
	if (verbose) cat('The data has ', sample.no, ' samples and ', otu.no, ' features will be tested!\n' )

	###############################################################################
	# Generate samples from posterior distribution (stacking it)
	if (is.prior == TRUE) {
		
		if (verbose) cat('Fitting beta mixture ...\n')
		comm.p <- apply(comm, 1, function (x) {			
					err1 <- try(res <- bbmix.fit.MM(x, depth), silent = TRUE)
					
					# Handle error
					if (class(err1) != 'try-error') {
						prop1.1 <- rbeta(sample.no * post.sample.no, shape1 = x + res$shape1.1, shape2 = res$shape1.2 + depth - x)
						prop1.2 <- rbeta(sample.no * post.sample.no, shape1 = x + res$shape2.1, shape2 = res$shape2.2 + depth - x)
						prop <- ifelse(runif(sample.no * post.sample.no) <= res$q1, prop1.1, prop1.2)
					} else {
						prop <- x / depth
						v <- var(prop)
						m <- mean(prop)
						
						a1 <- ((1 - m) / v - 1 / m) * m ^ 2
						a2 <- a1 * (1 / m - 1)
						
						if (is.na(a1) | a1 < 0) {
							# uniform prior
							prop <- rbeta(sample.no * post.sample.no, shape1 = x + 1, shape2 = otu.no + depth - x)
						} else {
							# beta prior
							prop <- rbeta(sample.no * post.sample.no, shape1 = x + a1, shape2 = a2 + depth - x)
						}	
					}
					return(prop)
				})
		
		comm.p <- t(comm.p)
		comm.p.list <- list()
		st <- 1
		end <- sample.no
		for (i in 1 : post.sample.no) {
			comm.p.list[[i]] <- comm.p[, st : end]
		# Normalization 
		#    comm.p.list[[i]] <- t(t(comm.p[, st : end]) / colSums(comm.p[, st : end]))
			st <- st + sample.no
			end <- end + sample.no
		}
	} else {
		comm.p.list <- list()
		comm.p.list[[1]] <- t(t(comm) / depth)
		post.sample.no <- 1
	}
	
	# Replace zeros or extremely small values for log calculation
	for (i in 1 : post.sample.no) {
		temp <- comm.p.list[[i]] 
		temp[temp <= min.prop] <- min.prop
		comm.p.list[[i]] <- temp
	}
	
	###############################################################################
	# Covariate space (including intercept)
	if (!is.null(strata)) {
		strata <- factor(strata)
	}
	
	if (is.null(adj.name)) {
		M0 <- model.matrix(~ 1, meta.dat) 
	} else {
		data0 <- meta.dat[, c(adj.name), drop = FALSE]
		M0 <- model.matrix( ~ ., data0) 
	}
	
	data1 <- meta.dat[, c(grp.name), drop = FALSE]
	M1 <-  model.matrix( ~ ., data1)[, -1, drop = FALSE]  # No intercept
	
	M01 <- cbind(M0, M1)
	
	# QR decompostion
	qrX0 <- qr(M0, tol = 1e-07)
	Q0 <- qr.Q(qrX0)
	Q0 <- Q0[, 1:qrX0$rank, drop = FALSE]
	H0 <- (Q0 %*% t(Q0))
	
	qrX1 <- qr(M1, tol = 1e-07)
	Q1 <- qr.Q(qrX1)
	Q1 <- Q1[, 1:qrX1$rank, drop = FALSE] 
	
	qrX01 <- qr(M01, tol = 1e-07)
	Q01 <- qr.Q(qrX01)
	Q01 <- Q01[, 1:qrX01$rank, drop = FALSE] 
	
	R0 <- as.matrix(resid(lm(Q1 ~ Q0 - 1)))
	
	pX0 <- ncol(Q0)
	pX1 <- ncol(Q1)
	pX01 <- ncol(Q01)
	
	df.model <- pX01 - pX0
	df.residual <- sample.no - pX01
	
	func.no <- length(link.func)
	
	###############################################################################
	if (verbose)  cat('Finding the references ...\n')
	
	ref.ind <- find.ref.z1(comm, M0, ref.pct = ref.pct)
#	rabund <- colMeans(t(comm) / depth)
#	ref.ind <- setdiff(ref.ind, order(rabund, decreasing = TRUE)[1 : 20])
	
	rabund <- colMeans(t(comm) / colSums(comm))
	size.factor <- colSums(comm[ref.ind, ])
	###############################################################################
	# Perform multiple stage normalization
	if (verbose) cat('Permutation testing ...\n')
	norm.ind <- NULL
	
	for (i in 1 : stage.no) {
		# Reference proportion
		divisor <- size.factor / depth
		
		# Create the giant Y matrix
		Y <- matrix(NA, sample.no, func.no * otu.no * post.sample.no)
		
		# Change order - func.no * otu.no
		for (k in 1 : post.sample.no) {
			for (j in 1 : func.no) {
				func <- link.func[[j]]
				comm.p <- comm.p.list[[k]]
				
				Y[, (k - 1) * func.no * otu.no + func.no * (0 : (otu.no - 1)) + j] <-
						func(t(comm.p) / divisor)  # No scaling first
				
			}
		}
		
		Y <- t(Y)
		TSS <- rowSums(Y^2)
		MSS01 <- rowSums((Y %*% Q01)^2)
		MSS0 <- rowSums((Y %*% Q0)^2)
		
		MSS <- (MSS01 - MSS0) 
		RSS <- (TSS - MSS01) 
		
		perm.ind <- getPermuteMatrix(perm.no, sample.no, strata = strata)
		perm.no <- nrow(perm.ind)
		
		MRSSp <- sapply(1 : perm.no, function(ii) {
					if (verbose) {
						if (ii %% 10 == 0) cat('.')
					}
					
					Rp <- R0[perm.ind[ii, ], , drop = FALSE]
					# Project to the reisdual space
					Rp <- Rp - H0 %*% Rp
					
					qrRp <- qr(Rp, tol = 1e-07)
					Q1p <- qr.Q(qrRp)
					Q1p <- Q1p[, 1:qrRp$rank, drop = FALSE]
					
					MSS01p <- MSS0 + rowSums((Y %*% Q1p)^2) 
					
					MSSp <- (MSS01p - MSS0) 
					RSSp <- (TSS - MSS01p) 
					
					c(MSSp, RSSp)
					
				})
		
		unit <- func.no * otu.no * post.sample.no
		MSSp <- MRSSp[1 : unit, ]
		RSSp <- MRSSp[(unit + 1) : (2 * unit), ]
		
		# EB is based on the aggregated RSS 
		RSS.m <- array(RSS, c(func.no, otu.no,  post.sample.no)) 
		RSS.m <- t(apply(RSS.m, c(1, 2), mean))  # otu.no by func.no
		
		F0.m <- array((MSS / df.model) / (RSS / df.residual), c(func.no, otu.no,  post.sample.no))
		F0.m <-  t(apply(F0.m, c(1, 2), mean))  # otu.no by func.no
		
		R2.m <- array(MSS / TSS, c(func.no, otu.no, post.sample.no))
		R2.m <- t(apply(R2.m, c(1, 2), mean))  # otu.no by func.no
		
		F0 <- (MSS / df.model)  /  (RSS  / df.residual)
		Fp <- (MSSp / df.model)  /  (RSSp  / df.residual)
		
		# Expectation of F0 and Fp
		F0 <- array(F0, c(func.no, otu.no, post.sample.no)) 
		Fp <- array(Fp, c(func.no, otu.no, post.sample.no, perm.no)) 
		
		F0 <- apply(F0, c(1, 2), mean)    # func.no * otu.no
		Fp <- apply(Fp, c(1, 2, 4), mean) # func.no * otu.no * perm.no
		
		
		###############################################################################
		# Omnibus test by taking maximum
		F0 <- apply(F0, 2, stats.combine.func)
		Fp <- apply(Fp, c(2, 3), stats.combine.func)  # otu.no by perm.no
		
		if (verbose) cat('\n')
		
		if (mean(is.na(F0)) >= 0.1) {
			warning('More than 10% observed F stats have NA! Please check! \n')
		}
		
		if (mean(is.na(Fp)) >= 0.1) {
			warning('More than 10% permuted F stats have NA! Please check! \n')
		}
		
		na.ind <- is.na(F0)
		F0 <- F0[!na.ind]
		Fp <- Fp[!na.ind, ]

		which.nan.ind <- which(!na.ind)
		###############################################################################
		p.raw <- rowMeans(cbind(Fp, F0) >= F0)
		p.adj.fdr <- perm_fdr_adj(F0, Fp)
		
		if (i == stage.no) {
			break
		} else {
#			if (mean(p.raw <= 0.05) >= excl.pct) {
#				ind <- order(p.raw)[1 : round(length(p.raw) * excl.pct)]
#			} else {
#				ind <- which(p.raw <= 0.05)
#			}

			ind <- order(p.adj.fdr)[1 : round(length(p.adj.fdr) * excl.pct)]

			
#			size.factor <- colSums(comm[setdiff(ref.ind, which.nan.ind[ind]), ])
	        size.factor <- colSums(comm[setdiff(1:nrow(comm), which.nan.ind[ind]), ])
			norm.ind <- cbind(norm.ind, !(1:nrow(comm) %in% setdiff(ref.ind, which.nan.ind[ind])))
		}
		
	}

	p.adj.fdr <- perm_fdr_adj(F0, Fp)

	p.raw <- na.pad(p.raw, na.ind)
	p.adj.fdr <- na.pad(p.adj.fdr, na.ind)
	
	names(p.raw) <- names(p.adj.fdr) <- rownames(R2.m) <- rownames(RSS.m) <- rownames(F0.m) <- row.names
	colnames(R2.m) <- colnames(F0.m) <- colnames(RSS.m) <- paste0('Func', 1 : func.no)
	
	if (is.fwer) {
		p.adj.fwer <- perm_fwer_adj2(F0, Fp)
		p.adj.fwer <- na.pad(p.adj.fwer, na.ind)
		names(p.adj.fwer)  <- row.names
	} else {
		p.adj.fwer <- NULL
	}
	
	if (!return.comm) {
		comm <- NULL
	}
	
	if (!return.perm.F) {
		Fp <- NULL
	}
	
	if (verbose) cat('Completed!\n')
	
	return(list(call = this.call, comm = comm, filter.ind = filter.ind,
					R2 = R2.m, F0 = F0.m, Fp = Fp, RSS = RSS.m, df.model = df.model, df.residual = df.residual,
					p.raw = p.raw, p.adj.fdr = p.adj.fdr,  p.adj.fwer = p.adj.fwer))
	
}


