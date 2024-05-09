## Scripts to control all simulations for Thresholding GLMs with SGD.
## Load libraries
library("progress")
## Define paths
experiment_path <- "~/repos/SGDInference"
outputs_path <- file.path(experiment_path, "outputs")
## Load functions
source("simulations.R")
source("exact_inference.R")
source("iterative_functions.R")

## Set any environment variables

## 1. Convergence to the true model
## Define settings

## Run simulations

## Save data

## 2. Oracle property
## Set hyperameters
nsim <- 500
settings <- expand.grid(family = c("binomial", "poisson"),
                        n = c(200, 500), 
                        p = c(10, 40, 100),
                        s = c(5))

make_file_name <- function(x) {
    paste0("oracle_sgd",
           "_", "family", x$family,
           "_", "n", x$n,
           "_", "p", x$p,
           "_", "s", x$s,
           ".rda")
}

settings <- settings |>
subset((n == 200 & p == 10) | (n == 200 & p == 40) | (n == 500 & p == 100))

rda_file <- file.path(outputs_path, "oracle_sgd.rda")

for (k in seq.int(nrow(settings))) {
    sett <- settings[k, ]
    sett_nm <- make_file_name(sett)
    cat("Iteration: ", k, "/", nrow(settings), sep = "")
    cat(sett_nm)
    if(file.exists(sett_nm)) {
        cat("Already run")
        next
    }
    pb <- progress_bar$new(total = nsim)
    coefs <- matrix(0, ncol = sett$p, nrow = nsim)
    ses <- matrix(0, ncol = sett$p, nrow = nsim)
    for (i in seq.int(nsim)) {
        set.seed(i)
        rr <- run_sim(model_name = model_name,
                      true_param = true_param,
                      nn = sett$n,
                      pp = sett$p,
                      ss = sett$s,
                      iB = iB,
                      B = B,
                      lvl = lvl,
                      reest = FALSE,
                      verbose = FALSE,
                      prog = FALSE)
        ## Get model
        Jhat <- with(rr, Jhat[nrow(Jhat), ])
        mm <- paste("X", (1:pp)[Jhat], sep = "")
        ## Re-estimate
        dd <- with(rr, data.frame(data$Y, data$X))
        attr(dd, "formula") <- paste0(c("y ~ -1", mm), collapse = " + ")
        names(dd) <- c("y", paste("X", 1:pp, sep = ""))
        mod <- glm(attr(dd, "formula"),
                   family = binomial(link = "logit"),
                   data = dd,
                   method = "brglmFit")
        coefs[i, Jhat] <- summary(mod)$coefficients[,1]
        X <- as.matrix(dd[, -c(1)])
        XJ <- X[, 1:ss] 
        eta <- X %*% rr$data$theta_star
        phi <- 1
        mu <- exp(-eta) / (1 + exp(-eta))
        W <- diag(as.vector(mu * (1 - mu)))
        FI <- t(XJ) %*% W %*% XJ / phi
        FIinv <- solve(FI)            
        ses[i, c(rep(TRUE, ss), rep(FALSE, pp - ss))] <- sqrt(diag(FIinv))
        pb$tick()
    }
}


## 3. Coverage of these confidence sets for one-pass SGD. ie iterate through all observations to get model, then repeat this and check coverage.

## 4. Confidence sets on the fly as an ad-hoc recommendation.


