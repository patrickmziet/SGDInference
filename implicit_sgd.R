## Compact code for running Implicit SGD .
## Panos Toulis, panos.toulis@chicagobooth.edu
require(Rcpp)
sourceCpp("isgd.cpp")
source("lib.R")
library(mvtnorm)
library(expm)
require(BH)
require(parallel)

isgd <- function(data, C=NA, gamma=1, npass=1, seed=runif(1) * 1e6, use_permutation=T) {
                                        # wrapper function
    if(npass>1) use_permutation=TRUE
    as.numeric(implicit_sgd(data$X, data$Y, C=C, model_name=data$model, gamma =gamma, n_pass=npass, seed=seed,
                            use_permutation=use_permutation))
}

isgd_0 <- function(theta0, data, C=NA, gamma=1, npass=1, seed=runif(1) * 1e6, use_permutation=T) {
                                        # wrapper function
    if(npass>1) use_permutation=TRUE
    as.numeric(implicit_sgd_0(theta0, data$X, data$Y, C=C, 
                              model_name=data$model, gamma =gamma, n_pass=npass, seed=seed,
                              use_permutation=use_permutation))
}# returns last iterate 

isgd_const <- function(data, C=NA, gamma=1, npass=1, seed=runif(1) * 1e6, use_permutation=T) {
                                        # wrapper function
    if(npass>1) use_permutation=TRUE
    as.numeric(implicit_sgd_const(data$X, data$Y, C=C, model_name=data$model, gamma=gamma, n_pass=npass, seed=seed,
                                  use_permutation=use_permutation))
}

aisgd <- function(data, C=NA, gamma=1, npass=1, seed=runif(1) * 1e6, use_permutation=T) {
                                        # wrapper function
    if(npass>1) use_permutation=TRUE
    as.numeric(implicit_sgd_avg(data$X, data$Y, C=C, model_name=data$model, gamma=gamma, n_pass=npass, seed=seed,
                                use_permutation=use_permutation))
}

aisgd_0 <- function(theta0, data, C=NA, gamma=1, npass=1, seed=runif(1) * 1e6, use_permutation=T) {
                                        # wrapper function
    if(npass>1) use_permutation=TRUE
    as.numeric(implicit_sgd_avg_0(theta0, data$X, data$Y, C=C, model_name=data$model, gamma=gamma, n_pass=npass, seed=seed,
                                  use_permutation=use_permutation))
}

isgd_adagrad <- function(data, seed=runif(1) * 1e6) {
                                        # wrapper function
    implicit_Adagrad(data$X, data$Y, model_name=data$model, seed=seed)
}

isgd_R <- function(data, theta0, gamma=1, verbose=TRUE) {
    N = nrow(data$X)
    p = ncol(data$X)
                                        # output
    theta_all = matrix(0, nrow=p, ncol=N)
    theta_prev = theta0
    
    for(i in 1:N) {
        xi = data$X[i, ]
        yi = data$Y[i]
        
        predi = sum(theta_prev * xi) 
        xi2 = sum(xi^2)
        
        gi = gamma / i
        
        implicitBound <- function(ksi) {
                                        # this returns the value  yn - h(theta_{n-1}' xn + xn^2 ξ)  -- for a GLM
                                        # this is a scalar.
            return(gi * (yi - data$glm_link(predi + xi2 * ksi)))
        }
        
                                        # 1. Define the search interval
        Bi = c(min(implicitBound(0), 0), max(implicitBound(0), 0))
                                        #Bounds2_n= c(min(implicitBound2(0), 0), max(implicitBound2(0), 0))
        
        implicit_fn <- function(u) {
            u  - implicitBound(u)
        }
        
                                        # 2. Solve implicit equation
        ksi_star = 0
        if(Bi[1] != Bi[2]) {
            ksi_star = uniroot(implicit_fn, interval=Bi)$root
        }
        ## Main update.
        theta_new = theta_prev + ksi_star * xi
        ## Store information.
        
        theta_all[, i] <- theta_new
                                        # dist0 = c(dist0, L2(theta_new, theta_all[, 1]))
                                        # mse = c(mse, L2(theta_new, data$theta_star))
        theta_prev = theta_new
        ## plot progress.
    }
    return(theta_all)
}

                                        #  R implementation of ISGD with conditioning
isgd_RC <- function(data, gamma=1, verbose=TRUE) {
    N = nrow(data$X)
    p = ncol(data$X)
                                        # output
    psi_all = matrix(0, nrow=p, ncol=N)
    V = diag(cov(data$X))
    X = data$X
                                        # for(j in 1:p) {
                                        #   X[, j] = X[, j] / sqrt(V[j])
                                        # }
    Y = data$Y
    C = 1 / sqrt(V)

    for(i in 2:N) {
        psi = psi_all[, i-1] # current param value
        xi = X[i, ]
        yi = Y[i]
        
                                        # predi = sum(xi * (V  * psi))
                                        # xi2 = sum(xi * (V * xi))
        predi = sum(xi * psi)
        xi2 = sum(xi * (C * xi))
        gi = gamma / (i-1)
        
        implicitBound <- function(ksi) {
                                        # this returns the value  yn - h(theta_{n-1}' xn + xn^2 ξ)  -- for a GLM
                                        # this is a scalar.
            return(gi * (yi - data$glm_link(predi + xi2 * ksi)))
        }
        
                                        # 1. Define the search interval
        Bi = c(min(implicitBound(0), 0), max(implicitBound(0), 0))
                                        #Bounds2_n= c(min(implicitBound2(0), 0), max(implicitBound2(0), 0))
        
        implicit_fn <- function(u) {
            u  - implicitBound(u)
        }
        
                                        # 2. Solve implicit equation
        ksi_star = 0
        if(Bi[1] != Bi[2]) {
            ksi_star = uniroot(implicit_fn, interval=Bi)$root
        }
        ## Main update.
        psi_all[, i]  = psi + ksi_star * (C * xi)
    }
                                        # theta_all = matrix(0, nrow=p, ncol=N)
                                        # for(i in 1:nrow(psi_all)) {
                                        #   theta_all[i,] = psi_all[i, ]* sqrt(V[i])
                                        # }
                                        # return(theta_all)
    return(psi_all)
}

adagrad_R <- function(data, gamma=1, verbose=TRUE) {
    N = nrow(data$X)
    p = ncol(data$X)
                                        # output
    theta_prev = rep(0, p)
    theta_new = NA
    Cn = 0.1 * rep(1, p)
    
    
    for(i in 1:N) {
        xi = data$X[i, ]
        yi = data$Y[i]
        predi = sum(theta_prev * xi) 
        yi_hat = data$glm_link(predi)
        gradi = (yi - yi_hat) * xi
        Cn = Cn + gradi^2
        theta_new = theta_prev + gamma * Cn^(-1/2) * gradi
        theta_prev = theta_new
    }
    return(theta_new)
}

sgd_R <- function(data, gamma, theta0=NULL) {
    N = nrow(data$X)
    p = ncol(data$X)
                                        # output
    if (is.null(theta0)) {
        theta_prev = rep(0, p)
    } else {
        theta_prev = theta0
    }
    theta_new = NA
    
    for(i in 1:N) {
        xi = data$X[i, ]
        yi = data$Y[i]
        predi = sum(theta_prev * xi) 
        yi_hat = data$glm_link(predi)
        gradi = (yi - yi_hat) * xi
        theta_new = theta_prev + gamma * gradi
        theta_prev = theta_new
    }
    return(theta_new)
}

