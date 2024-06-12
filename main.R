## Scripts to control all simulations for Thresholding GLMs with SGD.
## Load functions
source("simulations.R")
source("exact_inference.R")
source("iterative_functions.R")
source("support_functions.R")
## Load libraries
library("progress")

## Set any environment variables
## Define paths
experiment_path <- "~/repos/SGDInference"
outputs_path <- file.path(experiment_path, "outputs")

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
            mdls_j[i, ] <- Jhat
        }
        stats$mdls[[j]] <- mdls_j
    }
    ## Save setting
    save(stats, file = file.path(outputs_path, sett_nm))
}

## Simulation compared to old setting
settings <- expand.grid(rho = c(0),
                        n = c(500, 2500),
                        p = c(10, 40, 100),
                        s = c(5, 25),
                        family = c("binomial", "poisson"),
                        nu = c(5),
                        link = c("logit", "log"),
                        sim_axis = c("Observations"))

settings <- settings |>
subset(p < n & p > s) |>
subset(!(family == "poisson" & link == "logit")) |>
subset(!(family == "binomial" & link == "log")) |>
subset(!(sim_axis == "Observations" & n != 2500)) |>
subset(!(p == 100 & s == 25))
rownames(settings) <- 1:nrow(settings)
settings

save(settings, file = file.path(outputs_path, "cns_comp_sgd_settings.rda"))

nsimu <- 20 ## NSIM
## Load other package
pkg_path <- "~/repos/glmodsel-sub/glmodsel" # Replace with your own file path
experiment_glmodsel_path <- "~/repos/glmodsel-sub/scripts" # Replace with your own file path
outputs_glmodsel_path <- file.path(experiment_glmodsel_path, "outputs_paper")
devtools::document(pkg_path)
source(file.path(experiment_glmodsel_path, "support_functions.R"))

for (k in seq.int(nrow(settings))) {
    sett <- settings[k, ]
    beta <- numeric(sett$p)
    beta[1:sett$s] <- 1
    current_setting <- make_glm_setting(beta,
                                        link = sett$link,
                                        family = sett$family,
                                        s = sett$s,
                                        nu = sett$nu,
                                        n = sett$n,
                                        p = sett$p,
                                        rho = sett$rho,
                                        X = sett$X,
                                        sim_axis = sett$sim_axis)
    
    sim_vals <- round(exp(seq(log(current_setting$p * 2), log(4000), length = 10)))
    ## Increase sim_vals range when more sparse
    increase_range <- with(current_setting, (p == 100 & s == 5) | (p == 40 & s == 5))
    if (isTRUE(increase_range)) {
        sim_vals <- round(exp(seq(log(current_setting$p * 2), log(200000), length = 10)))
    }

    new_rda_file <- with(current_setting,
                         output_file(sim_axis,
                                     family,
                                     link,
                                     n,
                                     p,
                                     sum(beta != 0),
                                     rho,
                                     outputs_path))

    new_rda_file <- gsub(pattern = "summary_", replacement = "summary_sgd_", x = new_rda_file)
    cat("Setting", k, "of", nrow(settings),"\n")
    print(sett)
    if (file.exists(new_rda_file)) {
        print("Setting completed")
        next
    }

    rda_file <- with(current_setting,
                     output_file(sim_axis,
                                 family,
                                 link,
                                 n,
                                 p,
                                 sum(beta != 0),
                                 rho,
                                 outputs_glmodsel_path))

    summaries <- get_object(rda_file, "summaries")
    summs_soft <- NULL
    summs_hard <- NULL
    pb <- progress_bar$new(total = length(sim_vals))
    for (j in seq_along(sim_vals)) {
        current_setting$n <- sim_vals[j]
        selmat_soft <- matrix(FALSE, nrow = current_setting$p, ncol = nsimu)
        selmat_hard <- matrix(FALSE, nrow = current_setting$p, ncol = nsimu)

        ## Generate data from gen_data just to get the objects variables
        data <- gen_data(model_name = as.character(current_setting$family),
                         true_param = "hht",
                         N = current_setting$n,
                         p = current_setting$p,
                         s = current_setting$s)
        for (l in seq.int(nsimu)) {
            set.seed(l)
            dataset <- simulate_glm(current_setting)
            data$X <- as.matrix(dataset[, names(dataset) != "y"])
            attr(data$X, "dimnames") <- NULL
            data$Y <- as.integer(dataset[, names(dataset) == "y"])
            
            init_control <- list()
            init_control$gamma.method <- "ipower"
            init_control_def <- default_init(current_setting$p)
            for(ii in names(init_control)) init_control_def[[ii]] <- init_control[[ii]]
            init_control <- init_control_def
            sgd_control <- calculate_gammaStar_wrapper(data, init_control)
            gamma_star <- sgd_control$gamma

            theta_hat <- isgd_0(sgd_control$theta0, data, gamma=gamma_star, 
                                npass=1, use_permutation = FALSE)
            psi <- get_model_psi(data, theta_hat)
            V <- psi * gamma_star * with(current_setting, rep(1, p) / n)
            ## Compute pivots
            pivots <- abs(theta_hat) / sqrt(V)
            thresholds_soft <- optimise_threshold(betas = theta_hat,
                                                  stderrs = sqrt(V),
                                                  n = current_setting$n,
                                                  p = current_setting$p,
                                                  fixed = FALSE)
            thresholds_hard <- optimise_threshold(betas = theta_hat,
                                                  stderrs = sqrt(V),
                                                  n = current_setting$n,
                                                  p = current_setting$p,
                                                  fixed = TRUE)

            selmat_soft[,l] <- pivots > thresholds_hard
            selmat_hard[,l] <- pivots > thresholds_soft
        }

        attr(selmat_soft, "setting") <- current_setting
        attr(selmat_hard, "setting") <- current_setting
        attr(selmat_soft, "call") <- call(name = "sgd", fixed = FALSE)
        attr(selmat_hard, "call") <- call(name = "sgd", fixed = TRUE)
        class(selmat_soft) <- c("evaluate", class(selmat_soft))
        class(selmat_hard) <- c("evaluate", class(selmat_hard))
        summs_soft <- rbind(summs_soft, summary(selmat_soft)$metrics)
        summs_hard <- rbind(summs_hard, summary(selmat_hard)$metrics)
        pb$tick()
    }
    summs_soft$selector_call <- "sgd_soft"
    summs_hard$selector_call <- "sgd_hard"
    summs_soft$eval_expr <- "sgd_soft"
    summs_hard$eval_expr <- "sgd_hard"
    summaries <- rbind(summaries, summs_soft, summs_hard)
    save(current_setting, summaries, file = new_rda_file)
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
subset(!(n == 3000 & (family != "binomial" | s != 25)) | (family == "binomial" & p == 100 & s == 5 & n == 3000))

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

## Oracle example comparing with MSGLM chapter
## Load other package
## Define paths
pkg_path <- "~/repos/glmodsel-sub/glmodsel" # Replace with your own file path
experiment_glmodsel_path <- "~/repos/glmodsel-sub/scripts" # Replace with your own file path
outputs_glmodsel_path <- file.path(experiment_glmodsel_path, "outputs_paper")
outputs_glmodsel_path <- file.path(experiment_glmodsel_path, "outputs_paper")
devtools::document(pkg_path)
source(file.path(experiment_glmodsel_path, "support_functions.R"))

nsim <- 500
settings <- expand.grid(rho = c(0),
                        n = c(200, 500, 3000),
                        p = c(10, 40, 100),
                        s = c(5, 25),
                        family = c("binomial", "poisson"),
                        nu = c(5),
                        link = c("logit", "log"),
                        sim_axis = c("Observations"))

settings <- settings |>
subset(p < n & p > s) |>
subset(!(family == "poisson" & link == "logit")) |>
subset(!(family == "binomial" & link == "log")) |>
subset(!(p == 100 & s == 25))
rownames(settings) <- 1:nrow(settings)
settings

save(settings, file = file.path(outputs_path, "oracle_comp_sgd_settings.rda"))

get_coefs <- function(Jhat, dd) {
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
    summary(mod)$coefficients[,1]
}

for (k in seq.int(nrow(settings))) {
    sett <- settings[k, ]
    sett_nm <- make_file_name_msglm(sett, "oracle_sgd_comp")
    cat("Iteration: ", k, "/", nrow(settings), "\n", sep = "")
    print(sett)

    beta <- numeric(sett$p)
    beta[1:sett$s] <- 1
    current_setting <- make_glm_setting(beta,
                                        link = sett$link,
                                        family = sett$family,
                                        s = sett$s,
                                        nu = sett$nu,
                                        n = sett$n,
                                        p = sett$p,
                                        rho = sett$rho,
                                        X = sett$X,
                                        sim_axis = sett$sim_axis)

    rda_file <- file.path(outputs_path, sett_nm)
    if(file.exists(rda_file)) {
        print("Already run")
        next
    }

    data <- gen_data(model_name = as.character(current_setting$family),
                     true_param = "hht",
                     N = current_setting$n,
                     p = current_setting$p,
                     s = current_setting$s)

    pb <- progress_bar$new(total = nsim)
    fits <- list()
    fits[["sgd_soft"]] <- fits[["sgd_hard"]] <- fits[["ses"]] <- matrix(0, ncol = sett$p, nrow = nsim)
    system.time(
        for (i in seq.int(nsim)) {
            set.seed(i)
            dataset <- simulate_glm(current_setting)
            data$X <- as.matrix(dataset[, names(dataset) != "y"])
            attr(data$X, "dimnames") <- NULL
            data$Y <- dataset[, names(dataset) == "y"]
            
            init_control <- list()
            init_control$gamma.method <- "ipower"
            init_control_def <- default_init(current_setting$p)
            for(ii in names(init_control)) init_control_def[[ii]] <- init_control[[ii]]
            init_control <- init_control_def
            sgd_control <- calculate_gammaStar_wrapper(data, init_control)
            gamma_star <- sgd_control$gamma

            theta_hat <- isgd_0(sgd_control$theta0, data, gamma=gamma_star, 
                                npass=1, use_permutation = FALSE)
            psi <- get_model_psi(data, theta_hat)
            V <- psi * gamma_star * with(current_setting, rep(1, p) / n)
            ## Compute pivots
            pivots <- abs(theta_hat) / sqrt(V)
            thresholds_soft <- optimise_threshold(betas = theta_hat,
                                                  stderrs = sqrt(V),
                                                  n = current_setting$n,
                                                  p = current_setting$p,
                                                  fixed = FALSE)
            thresholds_hard <- optimise_threshold(betas = theta_hat,
                                                  stderrs = sqrt(V),
                                                  n = current_setting$n,
                                                  p = current_setting$p,
                                                  fixed = TRUE)

            Jhat_soft <- pivots > thresholds_soft
            Jhat_hard <- pivots > thresholds_hard

            dd <- data.frame(data$Y, data$X)
            if (sum(Jhat_soft) > 0) {
                fits[["sgd_soft"]][i, Jhat_soft] <- get_coefs(Jhat_soft, dd)
            }
            if (sum(Jhat_hard) > 0) {
                fits[["sgd_hard"]][i, Jhat_hard] <- get_coefs(Jhat_hard, dd)
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
            fits[["ses"]][i, with(sett, c(rep(TRUE, s), rep(FALSE, p - s)))] <- sqrt(diag(FIinv))
            pb$tick()
        }
    )
    ## Save data
    stats <- list()
    stats$fits <- fits
    stats$theta_star <- data$theta_star
    stats$sett <- sett
    stats$call <- sett_nm
    save(stats, file = rda_file)
}


## 3. Coverage of these confidence sets for one-pass SGD. ie iterate through all observations to get model, then repeat this and check coverage.
nsim <- 500
tuning <- "lmin" # "ipower" or "lmin"
fn_prefix <- "cov_sgd"
if (tuning == "lmin") fn_prefix <- paste0(fn_prefix, "_", tuning)

B <- 200 # Bootstrap iterations
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
subset(!(p == 100 & s == 25)) |>
subset(!((family == "binomial" & true_param == "hht-switch") | (family == "poisson" & true_param == "hht"))) |>
subset(!(n == 3000 & (family != "binomial" | s != 25)) | (family == "binomial" & p == 100 & s == 5 & n == 3000))

rownames(settings) <- 1:nrow(settings)
settings
save(settings, file = file.path(outputs_path, paste(fn_prefix, "settings.rda", sep = "_")))

## Run simulations
for (k in seq.int(nrow(settings))) {
    sett <- settings[k,]
    sett_nm <- make_file_name(sett, fn_prefix)
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
    pb <- progress_bar$new(total = nsim)
    for (i in seq.int(nsim)) {
        mdls_i <- matrix(FALSE, ncol = sett$p, nrow = B) # models for n_i in nsim
        res <- get_gamma(i = i,sett = sett, imth = tuning)
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
        pb$tick()
    }
    save(stats, file = file.path(outputs_path, sett_nm))
}

## 4. Confidence sets on the fly as an ad-hoc recommendation.
## Diabetes data
diabetes <- read.csv(file.path(experiment_path, "diabetes.csv"))
## Main effects
## change name of Outcome to "y"
names(diabetes)[which(names(diabetes) == "Outcome")] <- "y"
mf <- model.frame(y ~ ., data = diabetes)
X <- model.matrix(mf, data = diabetes)
Y <- model.response(mf)
diabetes_main <- diabetes
attr(diabetes_main, "formula") <- formula(mf)
n <- nrow(X)
p <- ncol(X)
## Compute gamma star
data <- list()
data$X <- scale(as.matrix(X))
data$X[,1] <- 1
data$Y <- Y
data$model <- "binomial"
## Run gen_data to get link functions
zz <- gen_data(model_name = "binomial",
                 N = 100,
                 p = 10,
                 s = 5,
                 true_param = "hht")
data$glm_link <- zz$glm_link

## Copy data for loop
data_cp <- data
iB <- 1
B <- 1
sq1 <- seq(1, iB)
sq2 <- seq.int(iB + 1, n)
chnks <- split(sq2, ceiling(seq_along(sq2) / B))
names(chnks) <- NULL
chnks_cp <- chnks
chnks[[1]] <- sq1
for (j in seq.int(length(chnks_cp))) chnks[[j + 1]] <- chnks_cp[[j]]
verbose <- FALSE
prog <- TRUE
if (prog) pb <- progress_bar$new(total = length(chnks))
## Compute gammastar
init_control <- list()
init_control$gamma.method <- "ipower"
init_control_def <- default_init(p)
for(ii in names(init_control)) init_control_def[[ii]] <- init_control[[ii]]
init_control <- init_control_def
sgd_control <- calculate_gammaStar_wrapper(data, init_control)
gamma_star <- sgd_control$gamma
mdls_hard <- NA
mdls_soft <- NA
col_names <- c("1", "\text{Preg.}", "\text{Gluc.}", "\text{BP}", "\text{Skin}", "\text{Insulin}", "\text{BMI}", "\text{Ped.}", "\text{Age}")
for (i in seq.int(chnks)) {
    if (verbose) print(paste0("Batch ", i, "| n = ", max(chnks[[i]])))
    ## reference data
    ref <- seq(1, max(chnks[[i]]))
    data_cp$X <- as.matrix(data$X[ref,])
    if (i == 1 & iB == 1) data_cp$X <- t(as.matrix(data$X[ref,]))
    data_cp$Y <- data$Y[ref]
    N <- nrow(data_cp$X)
    p <- ncol(data_cp$X)
    ## Very inefficient
    theta_hat <- isgd_0(sgd_control$theta0, data_cp, gamma=gamma_star, 
                        npass=1, use_permutation = FALSE)
    psi <- get_model_psi(data_cp, theta_hat)
    V <- psi * gamma_star * rep(1, p) / N
    ## Compute pivots
    pivots <- abs(theta_hat) / sqrt(V)
    thresholds_hard <- optimise_threshold(betas = theta_hat,
                                          stderrs = sqrt(V),
                                          n = N,
                                          p = p,
                                          fixed = TRUE)
    thresholds_soft <- optimise_threshold(betas = theta_hat,
                                          stderrs = sqrt(V),
                                          n = N,
                                          p = p,
                                          fixed = FALSE)
    Jhat_hard <- pivots > thresholds_hard  # Estimated model
    Jhat_hard[1] <- TRUE
    Jhat_soft <- pivots > thresholds_soft  # Estimated model
    Jhat_soft[1] <- TRUE
    mdls_hard[i] <- paste0("$\text{Diab.} ~ ",paste0((col_names)[Jhat_hard], collapse = " + "), "$")
    mdls_soft[i] <- paste0("$\text{Diab.} ~ ",paste0((col_names)[Jhat_soft], collapse = " + "), "$")
    ## Print diagnostics as iterations go on
    ## if (verbose) print(mdls[i])
    ## if (verbose) print_ci(mdls)
    if (prog) pb$tick()
}

save(mdls_hard, mdls_soft, file = file.path(outputs_path, "conf_sets.rda"))

## Hard thresholding
print_ci(mdls_hard)
## Soft thresholding
print_ci(mdls_soft)


## Coverage of diabetes example
res <- glm()
## Pick Model: Diab. ∼ 1 + Gluc. + BMI

## Set nsim
nsim <- 500
B <- 200
for (i in seq.int(nsim)) {
    set.seed(i)
    ## Generate new response variable

    ## Get confidence set on the fly

    ## Get confidence set via bootstrap glmodsel
    
}


