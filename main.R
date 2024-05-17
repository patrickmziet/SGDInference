## Scripts to control all simulations for Thresholding GLMs with SGD.
## Load functions
source("simulations.R")
source("exact_inference.R")
source("iterative_functions.R")
source("support_functions.R")
## Load libraries
library("progress")

## Define paths
experiment_path <- "~/repos/SGDInference"
outputs_path <- file.path(experiment_path, "outputs")
## Set any environment variables

## 1. Convergence to the true model
## Define settings
settings <- expand.grid(family = c("binomial", "poisson"),
                        sigma_x = c("id", "toeplitz", "equicor"),
                        true_param = c("hht", "hht-switch"),
                        fixed = c(TRUE, FALSE),
                        n = c(500), 
                        p = c(10, 40, 100),
                        s = c(5, 25))

settings <- settings |>
subset(p < n) |>
subset(s < p) |>
subset(!(p == 100 & s == 25)) |>
subset(!((family == "binomial" & true_param == "hht-switch") | (family == "poisson" & true_param == "hht")))
rownames(settings) <- 1:nrow(settings)
settings

save(settings, file = file.path(outputs_path, "cns_sgd_settings.rda"))
nsim <- 20 

## Run simulations
for (k in seq.int(nrow(settings))) {
    sett <- settings[k,]
    sett_nm <- make_file_name(sett, "cns_sgd")
    cat("Setting", k, "of", nrow(settings),"\n")
    print(sett)
    if(file.exists(file.path(outputs_path, sett_nm))) {
        print("Already run")
        next
    }
    ## Define sim_vals
    n_max <- 700
    if (sett$family == "binomial") {
        n_max <- 8000 # increase number of observations for binomial example
    }
    sim_vals <- round(exp(seq(log(sett$p * 2), log(n_max), length = 10)))
    stats <- list()
    stats$sett <- sett
    stats$call <- sett_nm
    stats$sim_vals <- sim_vals
    stats$mdls <- list()
    ## Loop along sim_vals
    for (j in seq_along(sim_vals)) {
        sett$n <- sim_vals[j]
        mdls_j <- matrix(FALSE, ncol = sett$p, nrow = nsim) # models for n_j in sim_vals
        for (i in seq.int(nsim)) {
            res <- get_gamma()
            sgd_control <- res$sgd_control
            gamma_star <- sgd_control$gamma
            data <- res$data
            theta_hat <- isgd_0(sgd_control$theta0, data, gamma=gamma_star, 
                                npass=1, use_permutation = FALSE)
            psi <- get_model_psi(data, theta_hat)
            V <- psi * gamma_star * rep(1, sett$p) / sett$n
            ## Compute pivots
            pivots <- abs(theta_hat) / sqrt(V)
            thresholds <- optimise_threshold(betas = theta_hat,
                                             stderrs = sqrt(V),
                                             n = sett$n,
                                             p = sett$p,
                                             fixed = sett$fixed) 
            Jhat <- pivots > thresholds  # Estimated model
            mdls_j[i, ] <- Jhat
        }
        stats$mdls[[j]] <- mdls_j
    }
    ## Save setting
    save(stats, file = file.path(outputs_path, sett_nm))
}


## 2. Oracle property
## Set hyperameters
nsim <- 500
settings <- expand.grid(family = c("binomial", "poisson"),
                        sigma_x = c("id", "toeplitz", "equicor"),
                        true_param = c("hht", "hht-switch"),
                        fixed = c(TRUE, FALSE),
                        n = c(200, 500, 3000), 
                        p = c(10, 40, 100),
                        s = c(5, 25))

settings <- settings |>
subset(p < n) |>
subset(s < p) |>
subset(n > 2 * p) |> # ensure n is not too small relative to p
subset(!(p == 100 & s == 25)) |>
subset(!((family == "binomial" & true_param == "hht-switch") | (family == "poisson" & true_param == "hht"))) |>
subset(!(n == 3000 & (family != "binomial" | s != 25)))
rownames(settings) <- 1:nrow(settings)
settings
save(settings, file = file.path(outputs_path, "oracle_sgd_settings.rda"))

settings_df <- as.data.frame(settings)


for (k in seq.int(nrow(settings))) {
    sett <- settings[k, ]
    sett_nm <- make_file_name(sett, "oracle_sgd")
    cat("Iteration: ", k, "/", nrow(settings), "\n", sep = "")
    print(sett_nm)
    if(file.exists(file.path(outputs_path, sett_nm))) {
        print("Already run")
        next
    }
    pb <- progress_bar$new(total = nsim)
    coefs <- matrix(0, ncol = sett$p, nrow = nsim)
    ses <- matrix(0, ncol = sett$p, nrow = nsim)
    for (i in seq.int(nsim)) {
        if (i < 189) next
        res <- get_gamma(i = i,sett = sett)
        sgd_control <- res$sgd_control
        gamma_star <- sgd_control$gamma
        data <- res$data
        theta_hat <- isgd_0(sgd_control$theta0, data, gamma=gamma_star, 
                            npass=1, use_permutation = FALSE)
        psi <- get_model_psi(data, theta_hat)
        V <- psi * gamma_star * rep(1, sett$p) / sett$n
        ## Compute pivots
        pivots <- abs(theta_hat) / sqrt(V)
        thresholds <- optimise_threshold(betas = theta_hat,
                                         stderrs = sqrt(V),
                                         n = sett$n,
                                         p = sett$p,
                                         fixed = sett$fixed) 
        Jhat <- pivots > thresholds  # Estimated model
        dd <- data.frame(data$Y, data$X)
        if (sum(Jhat) > 0) {
            mdls <- paste0("(",paste0((1:sett$p)[Jhat], collapse = ","),")")
            ## Get model
            mm <- paste("X", (1:sett$p)[Jhat], sep = "")
            ## Re-estimate
            attr(dd, "formula") <- paste0(c("y ~ -1", mm), collapse = " + ")
            names(dd) <- c("y", paste("X", 1:sett$p, sep = ""))
            mod <- glm(attr(dd, "formula"),
                       family = data$model,
                       data = dd,
                       method = "brglmFit")
            coefs[i, Jhat] <- summary(mod)$coefficients[,1]
        }
        X <- as.matrix(dd[, -c(1)])
        XJ <- X[, 1:sett$s] 
        eta <- X %*% data$theta_star
        phi <- 1
        if (sett$family == "binomial") {
            mu <- exp(-eta) / (1 + exp(-eta))
            W <- diag(as.vector(mu * (1 - mu)))
        } else if (sett$family == "poisson") {
            mu <- exp(eta)
            W <- diag(as.vector(mu))
        }
        FI <- t(XJ) %*% W %*% XJ / phi
        FIinv <- solve(FI)            
        ses[i, with(sett, c(rep(TRUE, s), rep(FALSE, p - s)))] <- sqrt(diag(FIinv))
        pb$tick()
    }
    ## Save data
    stats <- list()
    stats$coefs <- coefs
    stats$theta_star <- data$theta_star
    stats$sett <- sett
    stats$ses <- ses
    stats$call <- sett_nm
    save(stats, file = file.path(outputs_path, sett_nm))
}

## 3. Coverage of these confidence sets for one-pass SGD. ie iterate through all observations to get model, then repeat this and check coverage.
nsim <- 500
B <- 200 # Bootstrap iterations
settings <- expand.grid(family = c("binomial", "poisson"),
                        sigma_x = c("id", "toeplitz", "equicor"),
                        true_param = c("hht", "hht-switch"),
                        fixed = c(TRUE, FALSE),
                        n = c(200, 500), 
                        p = c(10, 40, 100),
                        s = c(5, 25))

settings <- settings |>
subset(p < n) |>
subset(s < p) |>
subset(!(p == 100 & s == 25)) |>
subset(!((family == "binomial" & true_param == "hht-switch") | (family == "poisson" & true_param == "hht")))
rownames(settings) <- 1:nrow(settings)
settings
save(settings, file = file.path(outputs_path, "cov_sgd_settings.rda"))


## Run simulations
for (k in seq.int(nrow(settings))) {
    sett <- settings[k,]
    sett_nm <- make_file_name(sett, "cov_sgd")
    cat("Setting", k, "of", nrow(settings),"\n")
    print(sett)
    if(file.exists(file.path(outputs_path, sett_nm))) {
        print("Already run")
        next
    }
    stats <- list()
    stats$sett <- sett
    stats$call <- sett_nm
    stats$mdls <- list()
    ## Loop along sim_vals
    for (i in seq.int(nsim)) {
        mdls_i <- matrix(FALSE, ncol = sett$p, nrow = B) # models for n_i in nsim
        res <- get_gamma()
        sgd_control <- res$sgd_control
        gamma_star <- sgd_control$gamma
        data <- res$data
        for(b in seq.int(B)) {
            data_cp <- data
            boot_ref <- sample(seq.int(sett$n), sett$n, replace = TRUE)
            data_cp$Y <- data_cp$Y[boot_ref]
            data_cp$X <- data_cp$X[boot_ref, ]
            theta_hat <- isgd_0(sgd_control$theta0, data_cp, gamma=gamma_star, 
                                npass=1, use_permutation = FALSE)
            psi <- get_model_psi(data_cp, theta_hat)
            V <- psi * gamma_star * rep(1, sett$p) / sett$n
            ## Compute pivots
            pivots <- abs(theta_hat) / sqrt(V)
            thresholds <- optimise_threshold(betas = theta_hat,
                                             stderrs = sqrt(V),
                                             n = sett$n,
                                             p = sett$p,
                                             fixed = sett$fixed) 
            Jhat <- pivots > thresholds  # Estimated model
            mdls_i[b, ] <- Jhat
        }
        stats$mdls[[i]] <- mdls_i
    }
    save(stats, file = file.path(outputs_path, sett_nm))
}

## 4. Confidence sets on the fly as an ad-hoc recommendation.
## APPLY TO REAL DATA set

settings <- expand.grid(family = c("binomial", "poisson"),
                        sigma_x = c("id", "toeplitz", "equicor"),
                        true_param = c("hht", "hht-switch"),
                        fixed = c(TRUE, FALSE),
                        n = c(200, 500), 
                        p = c(10, 40, 100),
                        s = c(5, 25))


