run_sim <- function(model_name, true_param, nn, pp, ss, B, lvl, verbose = TRUE) {
    ## Generate data
    data <- gen_data(model_name = model_name,
                     N = nn,
                     p = pp,
                     s = ss,
                     true_param = true_param)
    sq <- seq.int(nn)
    chnks <- split(sq, ceiling(seq_along(sq) / B))

    init_control <- list()
    init_control$gamma.method <- "ipower"
    init_control_def <- default_init(pp)
    for(n in names(init_control)) init_control_def[[n]] <- init_control[[n]]
    init_control <- init_control_def

    pivots <- matrix(0, nrow = length(chnks), ncol = pp)
    Jhat <- matrix(FALSE, nrow = length(chnks), ncol = pp)
    mdls <- NA
    data_cp <- data

    for (i in seq.int(chnks)) {
        if (verbose) print(paste0("Batch ", i, "\n"))
        ## reference data
        ref <- seq(1, max(chnks[[i]]))
        data_cp$X <- data$X[ref,]
        data_cp$Y <- data$Y[ref]
        N <- nrow(data_cp$X)
        p <- ncol(data_cp$X)
        ## Compute gammastar
        sgd_control <- calculate_gammaStar_wrapper(data_cp, init_control=init_control)
        gamma_star <- sgd_control$gamma
        ## Very inefficient
        theta_hat <- isgd_0(sgd_control$theta0, data_cp, gamma=gamma_star, 
                            npass=1, use_permutation = FALSE)
        psi <- get_model_psi(data_cp, theta_hat)

        V <- psi * gamma_star * rep(1, p) / N
        ## Compute pivots
        pivots[i,] <- abs(theta_hat) / sqrt(V)
        thresholds <- optimise_threshold(betas = theta_hat,
                                         stderrs = sqrt(V),
                                         n = N,
                                         p = p,
                                         fixed=TRUE) 
        Jhat[i, ] <- pivots[i, ] > thresholds  # Estimated model
        mdls[i] <- paste0("(",paste0((1:p)[Jhat[i, ]], collapse = ","),")")
        ## Print diagnostics as iterations go on
        if (verbose) print(mdls[i])
        if (verbose) print_ci(mdls)
    }
    return(list(mdls = mdls, pivots = pivots, theta_hat = theta_hat, Jhat = Jhat))
}

print_ci <- function(mm, lvl=0.95, verbose = TRUE) {
    prps <- round(table(mm) / length(mm), 4)
    result_df <- data.frame(Value = names(prps), Proportion = as.numeric(prps))
    result_df <- result_df[order(result_df$Proportion, decreasing = TRUE), ]
    result_df$cm_prop <- cumsum(result_df$Proportion)
    if (verbose) print(result_df, row.names = FALSE)
    ci_ref <- which(result_df$cm_prop <= lvl)
    if (length(ci_ref) != 0) {
        while (result_df$cm_prop[tail(ci_ref, 1)] < lvl) {
            ci_ref <- c(ci_ref, tail(ci_ref, 1) + 1)
        }
    } else {
        ci_ref[1] <- 1
    }
    ci <- result_df$Value[ci_ref]
    if (verbose) print(paste0("{", paste0(ci, collapse = ", "), "}"))
}

g <- function(n, u, p=NULL, type="naive", fixed=TRUE) {
    ## Threshold function
    switch(type,
           "naive"     = n^u,
           "enhanced1" = n^u + exp(1 / (1 - 2 * u) - 1) - 2,
           "enhanced2" = n^u - (1 - 4 * u) / (1 - 2 * u),
           "dBIC"      = sqrt((n^(1/n) - 1) * (n - p)),
           stop("unknown type"))
}

optimise_threshold <- function(betas, stderrs, n, p, fixed=TRUE, w_type="bicd", g_type="enhanced2") {
    ## Compute w's
    w <- function(..., type = "bicd") {
        dots <- list(...)
        with(dots,
             switch(type,
                    "bicd" = do.call("w_bicd", dots),
                    stop("unknown type")))
    }

    w_args <- function(type = "bicd") {
        switch(type,
               "bicd" = c("beta", "stderr", "n", "p"),
               stop("unknown type"))
    }

    w_bicd <- function(...) {
        with(list(...), {
            tsq <- (beta / stderr)^2
            dbic <- n * log(tsq / (n - p) + 1) - log(n)
            if (isTRUE(fixed)) dbic < 0 else pnorm(-dbic)
        })
    }

    res <- list(beta = betas,
                stderr = stderrs,
                n = n,
                p = p)
    w_fun <- function(...) w(..., type = w_type)
    w_args <- list()
    w_args$beta <- res$beta
    w_args$stderr <- res$stderr
    w_args$n <- res$n
    w_args$p <- res$p

    ws <- do.call("w_fun", w_args)

    ## Objects to use in objective
    beta_obj <- betas
    stderr_obj <- stderrs
    ws_obj <- ws
    n_obj <- n
    p_obj <- p
    ## Compute objective
    u_min <- obj_min <- numeric(p_obj)
    ## Function for g's
    g_fun <- function(n, u) g(n = n, u = u, type = g_type)
    ## Objective function
    type1 <- function(beta, stderr, g, w, n, p) {
        2 * pnorm(-g)
    }

    type2 <- function(beta, stderr, g, w, n, p, empirical = FALSE) {
        z <- beta / stderr
        pnorm(g - z) - pnorm(-g - z)
    }

    base_objective_empirical <- function(beta, stderr, g, w, n, p) {
        ## convex comb of Type I and Type II error in terms of u (i.e gamma in notes)
        w * type1(beta, stderr, g, w, n, p) +
            (1 - w) * type2(beta, stderr, g, w, n, p, empirical = TRUE)
    }
    objective <- function(beta, stderr, g, w, n, p) {
        base_objective_empirical(beta, stderr, g, w, n, p)
    }
    for (j in 1:p) {
        gj <- function(u) g_fun(n = n_obj, u = u)
        result <- optimize(function(u0) {
            objective(beta_obj[j], stderr_obj[j], gj(u0), ws_obj[j], n_obj, p_obj)
        }, lower = 0, upper = 0.5)
        u_min[j] <- result$minimum
        obj_min[j] <- result$objective
    }
    g_fun(n, u_min)
}
