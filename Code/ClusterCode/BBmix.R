# Estimate the posterior using beta mixture

require(rmutil)
###############################################################################################
bb.loglik <- function (paras, ct, dep, q1) {
	
	m1 <- exp(paras[1]) / (1 + exp(paras[1]))
	s1 <- exp(paras[2])
	
	-sum(dbetabinom(ct, dep, m1, s1, log = TRUE))
	
}

bb.fit  <- function (ct, dep,  winsor.qt = 1.0) {
	# Initialization
	if (mean(ct == 0) >= 0.95 | sum(ct != 0) < 5) {
		stop('The number of nonzero is too small to fit the model! Consider removing it from testing!\n')
	}
	
	# Initialization
	prop0 <- ct / dep
	qt <- quantile(prop0, winsor.qt)
	prop0[prop0 >= qt] <- qt
	var0 <- var(prop0)
	mean0 <- mean(prop0) 
	
	shape1 <- ((1 - mean0) / var0 - 1 / mean0) * mean0 ^ 2
	shape2 <- shape1 * (1 / mean0 - 1)
	
	m1 <- shape1 / (shape1 + shape2)
	s1 <- shape1 + shape2
	
	if (s1 < 0) {
		s1 <- 1
	}
	
	
	obj <- nlm(bb.loglik, c(log(m1 / (1 - m1)), log(s1)), ct = ct, dep = dep, q1 = q1)
	
	m1 <- exp(obj$estimate[1]) / (1 + exp(obj$estimate[1]))
	s1 <- exp(obj$estimate[2])
	
	shape1.1  <- m1 * s1
	shape1.2  <- s1 - shape1.1
	
	return(list(shape1.1 = shape1.1, shape1.2 = shape1.2))
	
}

bb.fit.MM <-  function (ct, dep,  winsor.qt = 1.0) {

	# Initialization
	if (mean(ct == 0) >= 0.95 | sum(ct != 0) < 5) {
		stop('The number of nonzero is too small to fit the model! Consider removing it from testing!\n')
	}
	
	# Initialization
	prop0 <- ct / dep
	qt <- quantile(prop0, winsor.qt)
	prop0[prop0 >= qt] <- qt
	var0 <- var(prop0)
	mean0 <- mean(prop0) 
	
	shape1 <- ((1 - mean0) / var0 - 1 / mean0) * mean0 ^ 2
	shape2 <- shape1 * (1 / mean0 - 1)
	
	
	return(list(shape1.1 = shape1, shape1.2 = shape2))
}

###############################################################################################
bbmix.fit.MM <- function (ct, dep,  nIter = 10, winsor.qt = 1.0) {
	
	if (mean(ct == 0) >= 0.95 | sum(ct != 0) < 5) {
		stop('The number of nonzero is too small to fit the model! Consider removing it from testing!\n')
	}
	
	# Initialization
	prop0 <- ct / dep
	qt <- quantile(prop0, winsor.qt)
	prop0[prop0 >= qt] <- qt
	
	var1 <- var(prop0)
	mean1 <- mean(prop0) / 2
	
	var2 <- var1
	mean2 <- 3 * mean(prop0) / 2
	
	pi <- 0.5
	for (i in 1 : nIter) {
		
		
		shape1.1 <- ((1 - mean1) / var1 - 1 / mean1) * mean1 ^ 2
		shape1.2 <- shape1.1 * (1 / mean1 - 1)
		
		shape2.1 <- ((1 - mean2) / var2 - 1 / mean2) * mean2 ^ 2
		shape2.2 <- shape2.1 * (1 / mean2 - 1)
		
		m1 <- shape1.1 / (shape1.1 + shape1.2)
		s1 <- shape1.1 + shape1.2
		
		m2 <- shape2.1 / (shape2.1 + shape2.2)
		s2 <- shape2.1 + shape2.2
		
		f1 <- dbetabinom(ct, dep, m1, s1)
		f2 <- dbetabinom(ct, dep, m2, s2)
		
		q1 <-  pi * f1 /  (pi * f1 + (1 - pi) * f2)
		q2 <-  1 - q1
		
		pi <- mean(q1)
		
		# Rough estimation
		mean1 <- sum(prop0 * q1) / sum(q1)
		var1 <- sum((prop0 - mean1)^2 * q1) / sum(q1)
		
		mean2 <- sum(prop0 * q2) / sum(q2)
		var2 <- sum((prop0 - mean2)^2 * q2) / sum(q2)
	}
	
	return(list(shape1.1 = shape1.1, shape1.2 = shape1.2, shape2.1 = shape2.1, shape2.2 = shape2.2, pi = pi, q1 = q1))
	
}

###############################################################################################
bbmix.loglik <- function (paras, ct, dep, q1) {
	
	m1 <- exp(paras[1]) / (1 + exp(paras[1]))
	s1 <- exp(paras[2])
	
	m2 <- exp(paras[3])/ (1 + exp(paras[3]))
	s2 <- exp(paras[4])
	
	-sum(q1 * dbetabinom(ct, dep, m1, s1, log = TRUE) + (1 - q1) * 
					dbetabinom(ct, dep, m2, s2, log = TRUE))
	
}

bbmix.fit <- function (ct, dep,  nIter = 10, winsor.qt = 0.97) {
	
	# Initialization
	if (mean(ct == 0) >= 0.95 | sum(ct != 0) < 5) {
		stop('The number of nonzero is too small to fit the model! Consider removing it from testing!\n')
	}
	
	# Initialization
	prop0 <- ct / dep
    qt <- quantile(prop0, winsor.qt)
	prop0[prop0 >= qt] <- qt
	var0 <- var(prop0)
	mean0 <- mean(prop0) 
	
	shape1 <- ((1 - mean0) / var0 - 1 / mean0) * mean0 ^ 2
	shape2 <- shape1 * (1 / mean0 - 1)
	
	m0 <- shape1 / (shape1 + shape2)
	s0 <- shape1 + shape2
	
	if (s0 < 0) {
		s0 <- 1
	}
	
	m1 <- 0.75 * m0
	m2 <- 1.25 * m0
	s1 <- s0
	s2 <- s0
    pi <- 0.5
	
	for (i in 1 : nIter) {
		
		# Expectation step
		f1 <- dbetabinom(ct, dep, m1, s1)
		f2 <- dbetabinom(ct, dep, m2, s2)
		
		q1 <-  pi * f1 /  (pi * f1 + (1 - pi) * f2)
		q2 <-  1 - q1
		
		# Maximization step
		pi <- mean(q1)
		
		obj <- nlm(bbmix.loglik, c(log(m1 / (1 - m1)), log(s1), log(m2 / (1 - m2)), log(s2)), ct = ct, dep = dep, q1 = q1)
		
		m1 <- exp(obj$estimate[1]) / (1 + exp(obj$estimate[1]))
		s1 <- exp(obj$estimate[2])
		m2 <- exp(obj$estimate[3]) / (1 + exp(obj$estimate[3]))
		s2 <- exp(obj$estimate[4])

	}
	shape1.1  <- m1 * s1
	shape1.2  <- s1 - shape1.1
	
	shape2.1  <- m2 * s2
	shape2.2  <- s2 - shape2.1
	

	return(list(shape1.1 = shape1.1, shape1.2 = shape1.2, shape2.1 = shape2.1, shape2.2 = shape2.2, pi = pi, q1 = q1))
	
}
###############################################################################################
zibb.loglik <- function (paras, ct, dep, q1) {
	
	m2 <- exp(paras[1]) / (1 + exp(paras[1]))
	s2 <- exp(paras[2])
	
	-sum((1 - q1) * dbetabinom(ct, dep, m2, s2, log = TRUE))
	
}

zibb.fit <- function (ct, dep,  nIter = 10, winsor.qt = 0.97) {
	
	# Initialization
	if (mean(ct == 0) >= 0.95 | sum(ct != 0) < 5) {
		stop('The number of nonzero is too small to fit the model! Consider removing it from testing!\n')
	}
	
	ct1 <- ct[ct != 0]
	dep1 <- dep[ct != 0]
	# Initialization
	prop0 <- ct1 / dep1
	qt <- quantile(prop0, winsor.qt)
	prop0[prop0 >= qt] <- qt
	var0 <- var(prop0)
	mean0 <- mean(prop0) 
	
	shape1 <- ((1 - mean0) / var0 - 1 / mean0) * mean0 ^ 2
	shape2 <- shape1 * (1 / mean0 - 1)
	
	m0 <- shape1 / (shape1 + shape2)
	s0 <- shape1 + shape2
	
	if (s0 < 0) {
		s0 <- 1
	}
	
	m2 <- 0.75 * m0
	s2 <- s0
	pi <- 0.75 * mean(ct == 0)
	
	for (i in 1 : nIter) {
		
		# Expectation step
		f1 <- ifelse(ct == 0, 1, 0)
		f2 <- dbetabinom(ct, dep, m2, s2)
		
		q1 <-  pi * f1 /  (pi * f1 + (1 - pi) * f2)
		q2 <-  1 - q1
		
		# Maximization step
		pi <- mean(q1)
		
		obj <- nlm(zibb.loglik, c(log(m2 / (1 - m2)), log(s2)), ct = ct, dep = dep, q1 = q1)
		
		m2 <- exp(obj$estimate[1]) / (1 + exp(obj$estimate[1]))
		s2 <- exp(obj$estimate[2])

	}
	shape2.1  <- m2 * s2
	shape2.2  <- s2 - shape2.1
	
	return(list(shape2.1 = shape2.1, shape2.2 = shape2.2,  pi = pi, q1 = q1))
	
}

###############################################################################################
zibb.fit.MM <- function (ct, dep,  nIter = 10, winsor.qt = 1.0) {
	
	# Initialization
	if (mean(ct == 0) >= 0.95 | sum(ct != 0) < 5) {
		stop('The number of nonzero is too small to fit the model! Consider removing it from testing!\n')
	}
	
	ct1 <- ct[ct != 0]
	dep1 <- dep[ct != 0]
	# Initialization
	prop0 <- ct1 / dep1
	qt <- quantile(prop0, winsor.qt)
	prop0[prop0 >= qt] <- qt
	var0 <- var(prop0)
	mean0 <- mean(prop0) 
	
	shape1 <- ((1 - mean0) / var0 - 1 / mean0) * mean0 ^ 2
	shape2 <- shape1 * (1 / mean0 - 1)
	
	m0 <- shape1 / (shape1 + shape2)
	s0 <- shape1 + shape2
	
	if (s0 < 0) {
		s0 <- 1
	}
	
	m2 <- 0.75 * m0
	s2 <- s0
	pi <- 0.75 * mean(ct == 0)
	
	##########################
	prop0 <- ct / dep
	qt <- quantile(prop0, winsor.qt)
	prop0[prop0 >= qt] <- qt
	##########################
	
	for (i in 1 : nIter) {
		
		# Expectation step
		f1 <- ifelse(ct == 0, 1, 0)
		f2 <- dbetabinom(ct, dep, m2, s2)
		
		q1 <-  pi * f1 /  (pi * f1 + (1 - pi) * f2)
		q2 <-  1 - q1
		
		# Maximization step
		pi <- mean(q1)
		
		mean2 <- sum(prop0 * q2) / sum(q2)
		var2 <- sum((prop0 - mean2)^2 * q2) / sum(q2)
		
		shape2.1 <- ((1 - mean2) / var2 - 1 / mean2) * mean2 ^ 2
		shape2.2 <- shape2.1 * (1 / mean2 - 1)
		
		m2 <- shape2.1 / (shape2.1 + shape2.2)
		s2 <- shape2.1 + shape2.2
		
	}
	shape2.1  <- m2 * s2
	shape2.2  <- s2 - shape2.1
	
	return(list(shape2.1 = shape2.1, shape2.2 = shape2.2,  pi = pi, q1 = q1))
	
}



###############################################################################################
bbmix.func <- function (x, res) {
	res$pi * dbeta(x, res$shape1.1, res$shape1.2) + (1 - res$pi) * dbeta(x, res$shape2.1, res$shape2.2)
}

bb.func <- function (x, res) {
	dbeta(x, res$shape1.1, res$shape1.2) 
}

zibb.func <- function (x, res) {
	(1 - res$pi) * dbeta(x, res$shape2.1, res$shape2.2) 
}

###############################################################################################

# Test
#temp <- load('~/Dropbox/Workspace/MayoClinic/Microbiome_data/Data/combo.RData')
#dim(otu.tab)
#
##otu.tab <- otu.tab[, sample(ncol(otu.tab), 20)]
#otu.tab <- otu.tab[rowMeans(otu.tab != 0) >= 0.15, ]
#dim(otu.tab)
#
#pdf('~/Dropbox/Workspace/MayoClinic/Methodology/2017_08_25_PermuteDAA/DepthConfounding/bbmix.fit.combo.pdf')
#par(mfrow=c(2, 2))
#dep <- colSums(otu.tab)
# 
#for (i in 1:nrow(otu.tab)) {
#	if (i %% 100 == 0) cat('.')
#	ct <- otu.tab[i, ]
#	
#	prop0 <- ct / dep
##	try(res0 <- zibb.fit.MM(ct, dep, winsor.qt = 1.0))
#	try(res0 <- bbmix.fit.MM(ct, dep, winsor.qt = 0.97))
#
#	hist(prop0, prob=TRUE, breaks = 20, main = 
#					paste0('MM:', paste0(paste0(c('a1=', 'b1='), formatC(unlist(res0[1:2]), digits=3)), collapse=','),
#							 '\n', paste0(paste0(c('a2=', 'b2=', 'p='), formatC(unlist(res0[3:5]), digits=3)), collapse=',')))
#	x <- seq(0.00001, max(prop0), len = 1000)
##	lines(x, zibb.func(x, res0),  col = 'red')
#	lines(x, bbmix.func(x, res0),   col = 'red')
#	
#	propE <- res0$q1 * (ct + res0$shape1.1) / (dep + res0$shape1.1 + res0$shape1.2) +
#			(1 - res0$q1) * (ct + res0$shape2.1) / (dep + res0$shape2.1 + res0$shape2.2)
#	
##	propE <- (ct + res0$shape1.1) / (dep + res0$shape1.1 + res0$shape1.2)
#	
##   propE <- (1 - res0$q1) * (ct + res0$shape2.1) / (dep + res0$shape2.1 + res0$shape2.2)
#	
#	plot(prop0, propE,  main = 'Observed vs posterior mean', cex = log10(dep) - 2.5, pch = 19, col = rgb(0.5, 0.5, 0.5, 0.5))
#	abline(0, 1)
##	
#}
#dev.off()
#
#
#
#require(GUniFrac)
#data(throat.otu.tab)
#otu.tab <- t(throat.otu.tab)
#otu.tab <- otu.tab[rowMeans(otu.tab != 0) > 0.15, ]
#
#pdf('~/Dropbox/Workspace/MayoClinic/Methodology/2017_08_25_PermuteDAA/DepthConfounding/ZIBB.fit.throat.pdf')
#par(mfrow=c(2, 2))
#
#dep <- colSums(otu.tab)
#
#for (i in 1:nrow(otu.tab)) {
#	if (i %% 100 == 0) cat('.')
#	ct <- otu.tab[i, ]
#	
#	prop0 <- ct / dep
#	try(res0 <- zibb.fit.MM(ct, dep, winsor.qt = 1.0))
##	try(res1 <- bbmix.fit(ct, dep, winsor.qt = 0.97))
#	
#	hist(prop0, prob=TRUE, breaks = 20, main = 
#					paste0('MM:', paste0(paste0(c('a1=', 'b1='), formatC(unlist(res0[1:2]), digits=3)), collapse=',')))
#	#		, '\n', paste0(paste0(c('a2=', 'b2=', 'p='), formatC(unlist(res0[3:5]), digits=3)), collapse=',')))
#	x <- seq(0.00001, max(prop0), len = 1000)
#	lines(x, zibb.func(x, res0),  col = 'red')
#	
##	propE <- res0$q1 * (ct + res0$shape1.1) / (dep + res0$shape1.1 + res0$shape1.2) +
##			(1 - res0$q1) * (ct + res0$shape2.1) / (dep + res0$shape2.1 + res0$shape2.2)
#	
##	propE <- (ct + res0$shape1.1) / (dep + res0$shape1.1 + res0$shape1.2)
#	
#	propE <- (1 - res0$q1) * (ct + res0$shape2.1) / (dep + res0$shape2.1 + res0$shape2.2)
#	
#	plot(prop0, propE,  main = 'Observed vs posterior mean', cex = log10(dep) - 2, pch = 19, col = rgb(0, 0, 1, 0.5))
#	abline(0, 1)
##	lines(x, bbmix.func(x, res1),   col = 'red')
#}
#dev.off()

